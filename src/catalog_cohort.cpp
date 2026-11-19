#include "catalog_cohort.hpp"

#include "duckdb/common/string_util.hpp"
#include "named_parameter_compat.hpp"

#include <algorithm>

namespace duckdb {

string ResolveAnalysisCatalog(ClientContext &context, const TableFunctionBindInput &input,
                              const vector<string> &rejected, const string &func_name) {
	auto &fs = FileSystem::GetFileSystem(context);
	enum class CatalogMode { AUTO, OFF, REQUIRE, PATH };
	CatalogMode mode = CatalogMode::AUTO;
	string catalog_path;
	const Value *cat_val = FindNamedParameter(input.named_parameters, "catalog");
	if (cat_val && !cat_val->IsNull()) {
		if (cat_val->type().id() == LogicalTypeId::BOOLEAN) {
			mode = cat_val->GetValue<bool>() ? CatalogMode::REQUIRE : CatalogMode::OFF;
		} else {
			mode = CatalogMode::PATH;
			catalog_path = cat_val->ToString();
			if (!input.inputs.empty()) {
				throw InvalidInputException("%s: catalog := '<path>' cannot be combined with a positional path (the "
				                            "catalog names the inputs)",
				                            func_name);
			}
		}
	}
	string arg;
	if (!input.inputs.empty() && !input.inputs[0].IsNull()) {
		arg = input.inputs[0].GetValue<string>();
	}
	if (mode == CatalogMode::AUTO || mode == CatalogMode::REQUIRE) {
		if (StringUtil::EndsWith(arg, ".pfile.json")) {
			catalog_path = arg;
		} else if (!arg.empty() && !StringUtil::EndsWith(arg, ".pgen") && arg.find("://") == string::npos &&
		           !StringUtil::StartsWith(arg, "pathmacro:") && arg.find_first_of("*?[") == string::npos) {
			catalog_path = ResolveExistingPath(context, fs, arg + ".pfile.json");
		}
		if (mode == CatalogMode::REQUIRE && catalog_path.empty()) {
			throw InvalidInputException("%s: catalog := true, but there is no catalog at '%s.pfile.json' (discovery "
			                            "looks only at local paths; name a remote catalog with catalog := '<path>')",
			                            func_name, arg.empty() ? string("(no positional path)") : arg);
		}
	}
	if (catalog_path.empty()) {
		if (input.inputs.empty()) {
			throw InvalidInputException("%s: needs a .pgen path, a catalog path, or catalog := '<path>'", func_name);
		}
		return catalog_path;
	}
	for (auto &key : rejected) {
		const Value *val = FindNamedParameter(input.named_parameters, key.c_str());
		if (!val || val->IsNull()) {
			continue;
		}
		if (key == "pvar" || key == "psam") {
			throw InvalidInputException(
			    "%s: a catalog (%s) cannot be combined with %s := (the catalog names every "
			    "cell's companions). To read the path as a .pgen instead, pass catalog := false.",
			    func_name, catalog_path, key);
		}
		throw InvalidInputException("%s: %s := %s is not supported with a catalog (%s)", func_name, key,
		                            val->ToString(), catalog_path);
	}
	return catalog_path;
}

const CatalogCohort::Segment &CatalogCohort::SegmentAt(uint32_t pos) const {
	auto it = std::upper_bound(segments.begin(), segments.end(), pos,
	                           [](uint32_t p, const Segment &s) { return p < s.scan_start; });
	D_ASSERT(it != segments.begin());
	return *(it - 1);
}

uint64_t CatalogCohort::CountMultiallelic() const {
	uint64_t ct = 0;
	for (auto &seg : segments) {
		ct += CountMultiallelicInRange(blocks[seg.block].variants, seg.start, seg.end);
	}
	return ct;
}

//! The variant and sample counts in a .pgen header.
static void ReadPgenCounts(const string &path, uint32_t &variant_ct, uint32_t &sample_ct, const string &func_name) {
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);
	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;
	plink2::PglErr err = plink2::PgfiInitPhase1(path.c_str(), nullptr, UINT32_MAX, UINT32_MAX, &header_ctrl, &pgfi,
	                                            &pgfi_alloc_cacheline_ct, errstr_buf);
	variant_ct = pgfi.raw_variant_ct;
	sample_ct = pgfi.raw_sample_ct;
	plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
	plink2::CleanupPgfi(&pgfi, &cleanup_err);
	if (err != plink2::kPglRetSuccess) {
		throw IOException("%s: failed to open '%s': %s", func_name, path, errstr_buf);
	}
}

unique_ptr<CatalogCohort> LoadCatalogCohort(ClientContext &context, const string &catalog_path,
                                            const string &func_name) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto grid = LoadCatalogGrid(context, fs, catalog_path, func_name);
	auto cohort = make_uniq<CatalogCohort>();
	cohort->path = grid.path;
	const idx_t K = grid.sample_block_ct;
	const idx_t V = grid.cells.size() / K;

	// A cell whose catalog entry omits a companion (a pvar_parquet catalog omits `pvar`) finds it beside
	// its .pgen, parquet companion first, as read_pfile does.
	auto companion = [&](const MatrixCell &cell, const string &declared, const vector<string> &exts, const char *what,
	                     const char *axis, idx_t block_idx) {
		if (!declared.empty()) {
			return declared;
		}
		auto found = FindCompanionFileWithParquet(context, fs, cell.prefix + ".pgen", exts);
		if (found.empty()) {
			throw InvalidInputException("%s: catalog '%s' names no %s for %s block %llu, and none is beside '%s.pgen'",
			                            func_name, grid.path, what, axis, static_cast<unsigned long long>(block_idx),
			                            cell.prefix);
		}
		return found;
	};

	// The samples: sample block sb's .psam is the one its cells name (row 0 holds one cell per block).
	vector<uint32_t> block_sample_ct(K);
	bool all_sexes = true;
	bool any_fids = false;
	vector<SampleInfo> parts(K);
	vector<string> psam_paths(K);
	for (idx_t sb = 0; sb < K; sb++) {
		auto psam = companion(grid.cells[sb], grid.cells[sb].psam, {".psam", ".fam"}, ".psam", "sample", sb);
		psam_paths[sb] = psam;
		parts[sb] = LoadSampleMetadata(context, psam);
		block_sample_ct[sb] = static_cast<uint32_t>(parts[sb].sample_ct);
		all_sexes = all_sexes && parts[sb].sexes.size() == parts[sb].sample_ct;
		any_fids = any_fids || !parts[sb].fids.empty();
	}
	auto &out = cohort->samples;
	out.sample_ct = 0;
	for (auto &part : parts) {
		out.iids.insert(out.iids.end(), part.iids.begin(), part.iids.end());
		if (any_fids) {
			if (part.fids.empty()) {
				out.fids.resize(out.fids.size() + part.sample_ct);
			} else {
				out.fids.insert(out.fids.end(), part.fids.begin(), part.fids.end());
			}
		}
		if (all_sexes) {
			out.sexes.insert(out.sexes.end(), part.sexes.begin(), part.sexes.end());
		}
		out.sample_ct += part.sample_ct;
	}
	cohort->sample_ct = static_cast<uint32_t>(out.sample_ct);

	// The variant blocks, in catalog order. Check every cell's .pgen against its .pvar and .psam.
	for (idx_t vb = 0; vb < V; vb++) {
		CatalogCohort::Block block;
		auto pvar = companion(grid.cells[vb * K], grid.cells[vb * K].pvar, {".pvar", ".bim"}, ".pvar", "variant", vb);
		block.variants = LoadVariantMetadata(context, pvar, func_name);
		block.variant_ct = static_cast<uint32_t>(block.variants.variant_ct);
		uint32_t offset = 0;
		for (idx_t sb = 0; sb < K; sb++) {
			PgenSampleShard shard;
			shard.pgen_path = grid.cells[vb * K + sb].prefix + ".pgen";
			PgenLocalizeGuard guard;
			LocalizePgenIfRequested(context, shard.pgen_path, guard);
			cohort->localize_guards.push_back(std::move(guard));
			const bool use_vfs = PgenIoUseVfs(context, shard.pgen_path);
			cohort->use_vfs = cohort->use_vfs || use_vfs;
			uint32_t pgen_variant_ct = 0;
			{
				PgenVfsScope pgen_vfs_scope(context, use_vfs, PgenVfsRole::HEADER);
				ReadPgenCounts(shard.pgen_path, pgen_variant_ct, shard.raw_sample_ct, func_name);
			}
			if (pgen_variant_ct != block.variant_ct) {
				throw InvalidInputException("%s: variant count mismatch: .pgen '%s' has %u variants, .pvar/.bim '%s' "
				                            "has %u variants",
				                            func_name, shard.pgen_path, pgen_variant_ct, pvar, block.variant_ct);
			}
			if (shard.raw_sample_ct != block_sample_ct[sb]) {
				throw InvalidInputException("%s: sample count mismatch: .pgen '%s' has %u samples, .psam/.fam '%s' has "
				                            "%u samples",
				                            func_name, shard.pgen_path, shard.raw_sample_ct, psam_paths[sb],
				                            block_sample_ct[sb]);
			}
			shard.sample_offset = offset;
			offset += shard.raw_sample_ct;
			block.shards.push_back(std::move(shard));
		}
		cohort->variant_ct += block.variant_ct;
		cohort->blocks.push_back(std::move(block));
	}
	return cohort;
}

ResolvedSelection ResolveCohortSelection(const Selection &selection, CatalogCohort &cohort, const string &func_name) {
	// The sample axis (and the rejection of an arbitrary `variants` set) is the shared handler's. The
	// variant axis is resolved here, block by block.
	Selection sample_axis = selection;
	sample_axis.has_region = false;
	sample_axis.has_variant_range = false;
	auto &first = cohort.blocks[0];
	ResolvedSelection resolved =
	    ResolveSelection(sample_axis, first.variants, &cohort.samples, cohort.sample_ct, first.variant_ct, func_name);
	if (selection.has_variant_range && selection.has_region) {
		throw InvalidInputException("%s: specify either `region` or `variant_range`, not both", func_name);
	}

	cohort.segments.clear();
	cohort.scan_ct = 0;
	uint32_t first_variant = 0; // cohort-wide ordinal of the block's first variant
	for (idx_t b = 0; b < cohort.blocks.size(); b++) {
		auto &block = cohort.blocks[b];
		uint32_t start = 0, end = block.variant_ct;
		if (selection.has_variant_range) {
			// [lo, hi) over cohort-wide ordinals, clamped to the cohort, as over one fileset.
			auto lo = std::min<int64_t>(selection.variant_range_lo, cohort.variant_ct);
			auto hi = std::min<int64_t>(selection.variant_range_hi, cohort.variant_ct);
			lo = std::max<int64_t>(lo, first_variant);
			hi = std::min<int64_t>(hi, static_cast<int64_t>(first_variant) + block.variant_ct);
			start = static_cast<uint32_t>(std::max<int64_t>(lo - first_variant, 0));
			end = static_cast<uint32_t>(std::max<int64_t>(hi - first_variant, start));
		} else if (selection.has_region) {
			auto range = ParseRegion(selection.region.GetValue<string>(), block.variants, func_name);
			start = range.start_idx;
			end = range.end_idx;
		}
		if (end > start) {
			CatalogCohort::Segment seg;
			seg.block = b;
			seg.start = start;
			seg.end = end;
			seg.scan_start = cohort.scan_ct;
			cohort.segments.push_back(seg);
			cohort.scan_ct += end - start;
		}
		first_variant += block.variant_ct;
	}
	return resolved;
}

ConcatPgenReader &CatalogCohortReader::Open(ClientContext &context, const CatalogCohort &cohort, idx_t block_idx,
                                            bool has_subset, const vector<uint32_t> &subset_sorted) {
	if (!reader || block != block_idx) {
		// A fresh reader replaces the previous block's, which closes its files.
		reader = make_uniq<ConcatPgenReader>();
		auto &b = cohort.blocks[block_idx];
		reader->Open(context, cohort.use_vfs, b.shards, b.variant_ct, has_subset, subset_sorted);
		block = block_idx;
	}
	return *reader;
}

} // namespace duckdb
