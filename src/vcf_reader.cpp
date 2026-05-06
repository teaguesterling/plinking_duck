#include "vcf_reader.hpp"
#include "plink_common.hpp"
#include "vcf_genotype_parse.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/string_util.hpp"

#include <pgenlib_misc.h>
#include <pgenlib_ffi_support.h>

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>

namespace duckdb {

// ---------------------------------------------------------------------------
// Buffered VFS line reader
// ---------------------------------------------------------------------------

static constexpr idx_t VCF_READ_BUF_SIZE = 65536;

struct BufferedLineReader {
	FileHandle &handle;
	std::unique_ptr<char[]> buf;
	idx_t buf_len = 0;
	idx_t buf_pos = 0;

	explicit BufferedLineReader(FileHandle &handle_p)
	    : handle(handle_p), buf(new char[VCF_READ_BUF_SIZE]) {
	}

	bool ReadLine(string &line) {
		line.clear();
		bool read_any = false;
		while (true) {
			if (buf_pos >= buf_len) {
				buf_len = handle.Read(buf.get(), VCF_READ_BUF_SIZE);
				buf_pos = 0;
				if (buf_len == 0) {
					return read_any;
				}
			}
			auto scan_start = buf_pos;
			while (buf_pos < buf_len) {
				char c = buf[buf_pos];
				if (c == '\n') {
					line.append(buf.get() + scan_start, buf_pos - scan_start);
					buf_pos++;
					if (!line.empty() && line.back() == '\r') {
						line.pop_back();
					}
					return true;
				}
				buf_pos++;
			}
			line.append(buf.get() + scan_start, buf_pos - scan_start);
			read_any = true;
		}
	}
};

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct VcfBindData : public TableFunctionData {
	string file_path;
	vector<string> sample_names;
	uint32_t sample_ct;

	GenotypeMode genotype_mode;
	bool include_phased;

	// Region filter
	string filter_chrom;
	int64_t filter_start = 0;
	int64_t filter_end = INT64_MAX;
	bool has_region_filter = false;

	// Quality filter
	int32_t min_gq = -1;
	int32_t min_dp = -1;
	int32_t max_dp = -1;

	// Half-call mode
	VcfHalfCallMode halfcall_mode = kHalfCallMissing;

	// Columns mode layout
	vector<string> genotype_column_names;
	idx_t columns_mode_first_geno_col = 0;
	idx_t columns_mode_geno_col_count = 0;

	// Column identifiers for projection pushdown
	static constexpr idx_t CHROM_COL = 0;
	static constexpr idx_t POS_COL = 1;
	static constexpr idx_t ID_COL = 2;
	static constexpr idx_t REF_COL = 3;
	static constexpr idx_t ALT_COL = 4;
	static constexpr idx_t GENOTYPES_COL = 5;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct VcfGlobalState : public GlobalTableFunctionState {
	mutex lock;
	unique_ptr<FileHandle> file_handle;
	unique_ptr<BufferedLineReader> reader;
	bool done = false;
	uint64_t multiallelic_skipped = 0;
	string line_buf;
	bool warned_multiallelic = false;

	// Projection pushdown
	vector<column_t> column_ids;
	bool need_genotypes = false;

	idx_t MaxThreads() const override {
		return 1; // Sequential text file reading
	}
};

// ---------------------------------------------------------------------------
// Local state
// ---------------------------------------------------------------------------

struct VcfLocalState : public LocalTableFunctionState {
	AlignedBuffer genovec_buf;
	AlignedBuffer phasepresent_buf;
	AlignedBuffer phaseinfo_buf;
	vector<int8_t> genotype_bytes;
	vector<int8_t> phased_pairs;
	bool initialized = false;
};

// ---------------------------------------------------------------------------
// Region parsing
// ---------------------------------------------------------------------------

static void ParseVcfRegion(const string &region_str, string &chrom, int64_t &start, int64_t &end) {
	auto colon_pos = region_str.find(':');
	if (colon_pos == string::npos) {
		chrom = region_str;
		start = 0;
		end = INT64_MAX;
		return;
	}

	chrom = region_str.substr(0, colon_pos);
	if (chrom.empty()) {
		throw InvalidInputException("read_plink_vcf: invalid region format '%s' (empty chromosome)", region_str);
	}

	auto range_str = region_str.substr(colon_pos + 1);
	auto dash_pos = range_str.find('-');
	if (dash_pos == string::npos) {
		throw InvalidInputException("read_plink_vcf: invalid region format '%s' (expected chr:start-end)", region_str);
	}

	auto start_str = range_str.substr(0, dash_pos);
	auto end_str = range_str.substr(dash_pos + 1);

	if (start_str.empty()) {
		throw InvalidInputException("read_plink_vcf: invalid region format '%s' (empty start position)", region_str);
	}

	char *parse_end;
	errno = 0;
	start = std::strtol(start_str.c_str(), &parse_end, 10);
	if (parse_end == start_str.c_str() || *parse_end != '\0' || errno != 0 || start < 0) {
		throw InvalidInputException("read_plink_vcf: invalid region start '%s' in '%s'", start_str, region_str);
	}

	if (end_str.empty()) {
		end = INT64_MAX;
	} else {
		errno = 0;
		end = std::strtol(end_str.c_str(), &parse_end, 10);
		if (parse_end == end_str.c_str() || *parse_end != '\0' || errno != 0 || end < 0) {
			throw InvalidInputException("read_plink_vcf: invalid region end '%s' in '%s'", end_str, region_str);
		}
	}

	if (start > end) {
		throw InvalidInputException("read_plink_vcf: region start (%lld) > end (%lld) in '%s'",
		                            static_cast<long long>(start), static_cast<long long>(end), region_str);
	}
}

// ---------------------------------------------------------------------------
// Halfcall mode parsing
// ---------------------------------------------------------------------------

static VcfHalfCallMode ParseHalfCallMode(const string &mode_str) {
	auto mode = StringUtil::Lower(mode_str);
	if (mode == "missing") {
		return kHalfCallMissing;
	} else if (mode == "reference") {
		return kHalfCallReference;
	} else if (mode == "haploid") {
		return kHalfCallHaploid;
	} else if (mode == "error") {
		return kHalfCallError;
	} else {
		throw InvalidInputException(
		    "read_plink_vcf: invalid halfcall value '%s' (expected 'missing', 'reference', 'haploid', or 'error')",
		    mode_str);
	}
}

// ---------------------------------------------------------------------------
// Bind
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> VcfBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<VcfBindData>();
	bind_data->file_path = input.inputs[0].GetValue<string>();

	// --- Parse the VCF header ---
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(bind_data->file_path, FileFlags::FILE_FLAGS_READ | FileCompressionType::AUTO_DETECT);
	BufferedLineReader header_reader(*handle);

	string line;
	bool found_header = false;
	while (header_reader.ReadLine(line)) {
		if (line.empty()) {
			continue;
		}
		// Skip ## metadata lines
		if (line.size() >= 2 && line[0] == '#' && line[1] == '#') {
			continue;
		}
		// Must be #CHROM line
		if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
			found_header = true;
			break;
		}
		// Data line before header
		throw InvalidInputException("read_plink_vcf: file '%s' has no #CHROM header line", bind_data->file_path);
	}

	if (!found_header) {
		throw InvalidInputException("read_plink_vcf: file '%s' is empty or has no #CHROM header line",
		                            bind_data->file_path);
	}

	// Parse #CHROM header: split on tabs, extract sample names from columns 9+
	vector<string> header_fields;
	{
		size_t start = 0;
		size_t pos = line.find('\t');
		while (pos != string::npos) {
			header_fields.push_back(line.substr(start, pos - start));
			start = pos + 1;
			pos = line.find('\t', start);
		}
		header_fields.push_back(line.substr(start));
	}

	// VCF fixed columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT (0-8)
	if (header_fields.size() < 10) {
		throw InvalidInputException("read_plink_vcf: file '%s' has no sample columns (need at least 10 tab-separated "
		                            "fields in #CHROM line, got %llu)",
		                            bind_data->file_path, static_cast<unsigned long long>(header_fields.size()));
	}

	for (size_t i = 9; i < header_fields.size(); i++) {
		bind_data->sample_names.push_back(header_fields[i]);
	}
	bind_data->sample_ct = static_cast<uint32_t>(bind_data->sample_names.size());

	// --- Parse named parameters ---
	string genotypes_str = "auto";
	bind_data->include_phased = false;

	for (auto &kv : input.named_parameters) {
		auto key = StringUtil::Lower(kv.first);
		if (key == "genotypes") {
			genotypes_str = kv.second.GetValue<string>();
		} else if (key == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (key == "region") {
			auto region_str = kv.second.GetValue<string>();
			ParseVcfRegion(region_str, bind_data->filter_chrom, bind_data->filter_start, bind_data->filter_end);
			bind_data->has_region_filter = true;
		} else if (key == "min_gq") {
			bind_data->min_gq = kv.second.GetValue<int32_t>();
		} else if (key == "min_dp") {
			bind_data->min_dp = kv.second.GetValue<int32_t>();
		} else if (key == "max_dp") {
			bind_data->max_dp = kv.second.GetValue<int32_t>();
		} else if (key == "halfcall") {
			bind_data->halfcall_mode = ParseHalfCallMode(kv.second.GetValue<string>());
		}
	}

	bind_data->genotype_mode = ResolveGenotypeMode(genotypes_str, bind_data->sample_ct, "read_plink_vcf");

	// Only ARRAY, LIST, and COLUMNS modes are supported — STRUCT/COUNTS/STATS require PgrGetCounts
	if (bind_data->genotype_mode != GenotypeMode::ARRAY &&
	    bind_data->genotype_mode != GenotypeMode::LIST &&
	    bind_data->genotype_mode != GenotypeMode::COLUMNS) {
		throw InvalidInputException("read_plink_vcf: genotypes := '%s' is not supported (use 'array', 'list', or 'columns')",
		                            genotypes_str);
	}

	if (bind_data->include_phased && IsAggregateGenotypeMode(bind_data->genotype_mode)) {
		throw InvalidInputException("read_plink_vcf: phased := true is incompatible with genotypes := '%s'",
		                            genotypes_str);
	}

	// --- Build return schema ---
	names.push_back("chrom");
	return_types.push_back(LogicalType::VARCHAR);

	names.push_back("pos");
	return_types.push_back(LogicalType::INTEGER);

	names.push_back("id");
	return_types.push_back(LogicalType::VARCHAR);

	names.push_back("ref");
	return_types.push_back(LogicalType::VARCHAR);

	names.push_back("alt");
	return_types.push_back(LogicalType::VARCHAR);

	if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
		// One column per sample
		bind_data->columns_mode_first_geno_col = 5;
		bind_data->columns_mode_geno_col_count = bind_data->sample_ct;

		for (uint32_t i = 0; i < bind_data->sample_ct; i++) {
			names.push_back(bind_data->sample_names[i]);
			if (bind_data->include_phased) {
				return_types.push_back(LogicalType::ARRAY(LogicalType::TINYINT, 2));
			} else {
				return_types.push_back(LogicalType::TINYINT);
			}
		}
		bind_data->genotype_column_names = bind_data->sample_names;
	} else {
		// Single genotypes column
		LogicalType geno_type;
		if (bind_data->include_phased) {
			auto pair_type = LogicalType::ARRAY(LogicalType::TINYINT, 2);
			if (bind_data->genotype_mode == GenotypeMode::ARRAY) {
				geno_type = LogicalType::ARRAY(pair_type, bind_data->sample_ct);
			} else {
				geno_type = LogicalType::LIST(pair_type);
			}
		} else {
			if (bind_data->genotype_mode == GenotypeMode::ARRAY) {
				geno_type = LogicalType::ARRAY(LogicalType::TINYINT, bind_data->sample_ct);
			} else {
				geno_type = LogicalType::LIST(LogicalType::TINYINT);
			}
		}
		names.push_back("genotypes");
		return_types.push_back(geno_type);
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> VcfInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<VcfBindData>();
	auto state = make_uniq<VcfGlobalState>();

	state->column_ids = input.column_ids;

	// Determine if genotypes need to be parsed (projection pushdown)
	state->need_genotypes = false;
	for (auto col_id : input.column_ids) {
		if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		if (bind_data.genotype_mode == GenotypeMode::COLUMNS) {
			if (col_id >= bind_data.columns_mode_first_geno_col &&
			    col_id < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
				state->need_genotypes = true;
				break;
			}
		} else if (col_id == VcfBindData::GENOTYPES_COL) {
			state->need_genotypes = true;
			break;
		}
	}

	auto &fs = FileSystem::GetFileSystem(context);
	state->file_handle = fs.OpenFile(bind_data.file_path, FileFlags::FILE_FLAGS_READ | FileCompressionType::AUTO_DETECT);
	state->reader = make_uniq<BufferedLineReader>(*state->file_handle);

	string line;
	while (state->reader->ReadLine(line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '#') {
			continue;
		}
		// First data line — store for scan to pick up
		state->line_buf = line;
		break;
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> VcfInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                        GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<VcfBindData>();
	auto state = make_uniq<VcfLocalState>();

	auto sample_ct = bind_data.sample_ct;

	// Allocate genovec buffer
	auto genovec_word_ct = plink2::NypCtToAlignedWordCt(sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));

	if (bind_data.include_phased) {
		auto bitvec_word_ct = plink2::BitCtToAlignedWordCt(sample_ct);
		state->phasepresent_buf.Allocate(bitvec_word_ct * sizeof(uintptr_t));
		state->phaseinfo_buf.Allocate(bitvec_word_ct * sizeof(uintptr_t));
	}

	state->genotype_bytes.resize(sample_ct);
	if (bind_data.include_phased) {
		state->phased_pairs.resize(sample_ct * 2);
	}

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// FORMAT field parsing: find GT position and quality field skips
// ---------------------------------------------------------------------------

struct FormatParseResult {
	int32_t gt_pos = -1;
	int32_t gq_pos = -1;
	int32_t dp_pos = -1;
};

static FormatParseResult ParseFormatField(const char *format_start, size_t format_len) {
	FormatParseResult result;
	int32_t field_idx = 0;
	const char *p = format_start;
	const char *end = format_start + format_len;
	const char *field_start = p;

	while (true) {
		if (p == end || *p == ':') {
			size_t field_len = static_cast<size_t>(p - field_start);
			if (field_len == 2 && field_start[0] == 'G' && field_start[1] == 'T') {
				result.gt_pos = field_idx;
			} else if (field_len == 2 && field_start[0] == 'G' && field_start[1] == 'Q') {
				result.gq_pos = field_idx;
			} else if (field_len == 2 && field_start[0] == 'D' && field_start[1] == 'P') {
				result.dp_pos = field_idx;
			}
			field_idx++;
			if (p == end) {
				break;
			}
			field_start = p + 1;
		}
		p++;
	}

	return result;
}

// ---------------------------------------------------------------------------
// Build VcfParseContext for a given line's FORMAT field
// ---------------------------------------------------------------------------

static void BuildParseContext(VcfParseContext &ctx, const VcfBindData &bind_data, const FormatParseResult &format_info) {
	ctx.sample_ct = bind_data.sample_ct;
	ctx.halfcall_mode = bind_data.halfcall_mode;
	ctx.qual_field_ct = 0;
	ctx.qual_field_skips[0] = 0;
	ctx.qual_field_skips[1] = 0;
	ctx.qual_line_mins[0] = INT32_MIN;
	ctx.qual_line_maxs[0] = INT32_MAX;
	ctx.qual_line_mins[1] = INT32_MIN;
	ctx.qual_line_maxs[1] = INT32_MAX;

	if (format_info.gt_pos != 0) {
		return;
	}

	// Collect requested quality fields with their FORMAT positions, then sort
	// by position so skip deltas are always non-negative.
	struct QualField {
		int32_t format_pos;
		int32_t min_val;
		int32_t max_val;
	};
	QualField fields[2];
	uint32_t n_fields = 0;

	if (bind_data.min_gq >= 0 && format_info.gq_pos > 0) {
		fields[n_fields++] = {format_info.gq_pos, bind_data.min_gq, INT32_MAX};
	}
	if ((bind_data.min_dp >= 0 || bind_data.max_dp >= 0) && format_info.dp_pos > 0 && n_fields < 2) {
		fields[n_fields++] = {format_info.dp_pos,
		                      bind_data.min_dp >= 0 ? bind_data.min_dp : INT32_MIN,
		                      bind_data.max_dp >= 0 ? bind_data.max_dp : INT32_MAX};
	}

	if (n_fields == 2 && fields[1].format_pos < fields[0].format_pos) {
		std::swap(fields[0], fields[1]);
	}

	int32_t last_pos = 0;
	for (uint32_t i = 0; i < n_fields; i++) {
		ctx.qual_field_skips[i] = static_cast<uint32_t>(fields[i].format_pos - last_pos);
		ctx.qual_line_mins[i] = fields[i].min_val;
		ctx.qual_line_maxs[i] = fields[i].max_val;
		last_pos = fields[i].format_pos;
	}
	ctx.qual_field_ct = n_fields;
}

// ---------------------------------------------------------------------------
// Tab-split helper: find tab positions in line (avoids allocation)
// ---------------------------------------------------------------------------

// Finds the start of each tab-delimited field and returns field-start offsets.
static void FindTabFields(const string &line, vector<size_t> &tab_positions) {
	tab_positions.clear();
	tab_positions.push_back(0); // first field starts at 0
	for (size_t i = 0; i < line.size(); i++) {
		if (line[i] == '\t') {
			tab_positions.push_back(i + 1);
		}
	}
}

// Get a field from the line given tab positions (returns a string)
static inline string GetField(const string &line, const vector<size_t> &tab_positions, size_t field_idx) {
	size_t start = tab_positions[field_idx];
	size_t end;
	if (field_idx + 1 < tab_positions.size()) {
		end = tab_positions[field_idx + 1] - 1; // -1 for the tab character
	} else {
		end = line.size();
	}
	return line.substr(start, end - start);
}

// ---------------------------------------------------------------------------
// POS parsing helper with validation
// ---------------------------------------------------------------------------

static int32_t ParseVcfPos(const string &line, const vector<size_t> &tab_positions) {
	size_t start = tab_positions[1];
	size_t end = (tab_positions.size() > 2) ? tab_positions[2] - 1 : line.size();
	char *parse_end;
	errno = 0;
	long val = std::strtol(line.c_str() + start, &parse_end, 10);
	if (parse_end == line.c_str() + start || errno != 0 || val < 0 || val > INT32_MAX) {
		throw IOException("read_plink_vcf: invalid POS value '%s'",
		                  line.substr(start, end - start));
	}
	return static_cast<int32_t>(val);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void VcfScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<VcfBindData>();
	auto &gstate = data_p.global_state->Cast<VcfGlobalState>();
	auto &lstate = data_p.local_state->Cast<VcfLocalState>();

	auto &column_ids = gstate.column_ids;
	auto sample_ct = bind_data.sample_ct;

	idx_t row_count = 0;
	vector<size_t> tab_positions;
	tab_positions.reserve(9 + sample_ct + 1);

	VcfParseContext parse_ctx;

	while (row_count < STANDARD_VECTOR_SIZE) {
		string line;
		{
			lock_guard<mutex> guard(gstate.lock);
			if (gstate.done) {
				break;
			}

			if (!gstate.line_buf.empty()) {
				line = std::move(gstate.line_buf);
				gstate.line_buf.clear();
			} else {
				if (!gstate.reader->ReadLine(line)) {
					gstate.done = true;
					break;
				}
			}
		}

		if (line.empty()) {
			continue;
		}

		FindTabFields(line, tab_positions);

		size_t min_fields = 9 + sample_ct;
		if (tab_positions.size() < min_fields) {
			continue;
		}

		// Region filter — uses raw line positions to avoid string allocation
		if (bind_data.has_region_filter) {
			auto chrom = GetField(line, tab_positions, 0);
			if (chrom != bind_data.filter_chrom) {
				continue;
			}
			int32_t pos_val = ParseVcfPos(line, tab_positions);
			if (pos_val < bind_data.filter_start || pos_val > bind_data.filter_end) {
				continue;
			}
		}

		// Multiallelic check: skip if ALT contains a comma
		{
			size_t alt_start = tab_positions[4];
			size_t alt_end = (tab_positions.size() > 5) ? tab_positions[5] - 1 : line.size();
			if (line.find(',', alt_start) < alt_end) {
				lock_guard<mutex> guard(gstate.lock);
				gstate.multiallelic_skipped++;
				continue;
			}
		}

		// Parse genotypes if needed
		if (gstate.need_genotypes) {
			auto format_field = GetField(line, tab_positions, 8);
			auto format_info = ParseFormatField(format_field.data(), format_field.size());

			if (format_info.gt_pos < 0) {
				std::fill(lstate.genotype_bytes.begin(), lstate.genotype_bytes.end(), static_cast<int8_t>(-9));
				if (bind_data.include_phased) {
					std::fill(lstate.phased_pairs.begin(), lstate.phased_pairs.end(), static_cast<int8_t>(-9));
				}
			} else if (format_info.gt_pos != 0) {
				throw IOException("read_plink_vcf: GT must be the first FORMAT subfield, got '%s'",
				                  format_field);
			} else {
				BuildParseContext(parse_ctx, bind_data, format_info);

				const char *sample_data_ptr = line.c_str() + tab_positions[9];
				auto *genovec = lstate.genovec_buf.As<uintptr_t>();

				VcfGenoParseResult result;
				if (bind_data.include_phased) {
					auto *phasepresent = lstate.phasepresent_buf.As<uintptr_t>();
					auto *phaseinfo = lstate.phaseinfo_buf.As<uintptr_t>();
					auto bitvec_word_ct = plink2::BitCtToAlignedWordCt(sample_ct);
					memset(phasepresent, 0, bitvec_word_ct * sizeof(uintptr_t));
					memset(phaseinfo, 0, bitvec_word_ct * sizeof(uintptr_t));
					result = ParsePhasedBiallelicGT(parse_ctx, sample_data_ptr, genovec, phasepresent, phaseinfo);
				} else {
					result = ParseUnphasedBiallelicGT(parse_ctx, sample_data_ptr, genovec);
				}

				if (result != VcfGenoParseResult::OK) {
					auto chrom = GetField(line, tab_positions, 0);
					auto pos = GetField(line, tab_positions, 1);
					throw IOException("read_plink_vcf: failed to parse genotypes at %s:%s (error: %d)",
					                  chrom, pos, static_cast<int>(result));
				}

				plink2::GenoarrToBytesMinus9(genovec, sample_ct, lstate.genotype_bytes.data());

				if (bind_data.include_phased) {
					UnpackPhasedGenotypes(lstate.genotype_bytes.data(), lstate.phasepresent_buf.As<uintptr_t>(),
					                      lstate.phaseinfo_buf.As<uintptr_t>(), sample_ct, lstate.phased_pairs.data());
				}
			}
		}

		// Fill output vectors — fixed columns first (shared), then genotype columns
		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];

			if (file_col == VcfBindData::CHROM_COL) {
				FlatVector::GetData<string_t>(vec)[row_count] =
				    StringVector::AddString(vec, GetField(line, tab_positions, 0));
			} else if (file_col == VcfBindData::POS_COL) {
				FlatVector::GetData<int32_t>(vec)[row_count] = ParseVcfPos(line, tab_positions);
			} else if (file_col == VcfBindData::ID_COL) {
				auto id_field = GetField(line, tab_positions, 2);
				if (id_field == ".") {
					FlatVector::SetNull(vec, row_count, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row_count] =
					    StringVector::AddString(vec, id_field);
				}
			} else if (file_col == VcfBindData::REF_COL) {
				FlatVector::GetData<string_t>(vec)[row_count] =
				    StringVector::AddString(vec, GetField(line, tab_positions, 3));
			} else if (file_col == VcfBindData::ALT_COL) {
				auto alt_field = GetField(line, tab_positions, 4);
				if (alt_field == ".") {
					FlatVector::SetNull(vec, row_count, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row_count] =
					    StringVector::AddString(vec, alt_field);
				}
			} else if (bind_data.genotype_mode == GenotypeMode::COLUMNS &&
			           file_col >= bind_data.columns_mode_first_geno_col &&
			           file_col < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
				auto sample_idx = static_cast<uint32_t>(file_col - bind_data.columns_mode_first_geno_col);
				if (bind_data.include_phased) {
					auto &child = ArrayVector::GetEntry(vec);
					auto *child_data = FlatVector::GetData<int8_t>(child);
					idx_t base = row_count * 2;
					int8_t a1 = lstate.phased_pairs[sample_idx * 2];
					int8_t a2 = lstate.phased_pairs[sample_idx * 2 + 1];
					if (a1 == -9) {
						FlatVector::SetNull(vec, row_count, true);
						child_data[base] = 0;
						child_data[base + 1] = 0;
					} else {
						child_data[base] = a1;
						child_data[base + 1] = a2;
					}
				} else {
					int8_t geno = lstate.genotype_bytes[sample_idx];
					if (geno == -9) {
						FlatVector::SetNull(vec, row_count, true);
					} else {
						FlatVector::GetData<int8_t>(vec)[row_count] = geno;
					}
				}
			} else if (file_col == VcfBindData::GENOTYPES_COL) {
				FillGenotypeVector(vec, row_count, bind_data.genotype_mode, sample_ct,
				                   lstate.genotype_bytes.data(),
				                   bind_data.include_phased ? lstate.phased_pairs.data() : nullptr,
				                   bind_data.include_phased);
			}
		}

		row_count++;
	}

	// Emit warning for skipped multiallelic variants (once per query)
	if (row_count == 0 && gstate.multiallelic_skipped > 0 && !gstate.warned_multiallelic) {
		gstate.warned_multiallelic = true;
		Printer::Print(StringUtil::Format(
		    "read_plink_vcf: skipped %llu multiallelic variant(s)",
		    static_cast<unsigned long long>(gstate.multiallelic_skipped)));
	}

	output.SetCardinality(row_count);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkVcfReader(ExtensionLoader &loader) {
	TableFunction func("read_plink_vcf", {LogicalType::VARCHAR}, VcfScan, VcfBind, VcfInitGlobal, VcfInitLocal);
	func.named_parameters["genotypes"] = LogicalType::VARCHAR;
	func.named_parameters["phased"] = LogicalType::BOOLEAN;
	func.named_parameters["region"] = LogicalType::VARCHAR;
	func.named_parameters["min_gq"] = LogicalType::INTEGER;
	func.named_parameters["min_dp"] = LogicalType::INTEGER;
	func.named_parameters["max_dp"] = LogicalType::INTEGER;
	func.named_parameters["halfcall"] = LogicalType::VARCHAR;
	func.projection_pushdown = true;

	loader.RegisterFunction(func);
}

} // namespace duckdb
