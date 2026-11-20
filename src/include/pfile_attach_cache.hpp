#pragma once

#include "duckdb.hpp"
#include "duckdb/storage/object_cache.hpp"
#include "plink_common.hpp"
#include "pfile_sample_metadata.hpp"

#include <limits>
#include <mutex>

namespace duckdb {

//! The metadata cache of one ATTACH ... (TYPE pfile_catalog) (design Part 2). It lives in DuckDB's
//! ObjectCache under a key unique to that ATTACH; the views pass the key to read_pfile as
//! attach_cache := '<key>', and DETACH deletes the entry. The ATTACH is a snapshot: nothing here is
//! re-validated against the files. Every access takes `lock`, and a Get copies out, so two queries
//! binding at once never share state. A missing entry (evicted under memory pressure) means no cache.
class PfileAttachCache : public ObjectCacheEntry {
public:
	//! What discovery resolved for one prefix (LoadPfileSource, no overrides).
	struct CachedCompanions {
		string pgen;
		string pvar;
		string psam;
		string eff_prefix;
	};

	static string ObjectType() {
		return "plinking_attach_cache";
	}
	string GetObjectType() override {
		return ObjectType();
	}
	optional_idx GetEstimatedCacheMemory() const override {
		std::lock_guard<std::mutex> guard(lock);
		return optional_idx(bytes);
	}

	bool GetGrid(const string &catalog_path, CatalogGrid &out) const {
		std::lock_guard<std::mutex> guard(lock);
		auto it = grids.find(catalog_path);
		if (it == grids.end()) {
			return false;
		}
		out = it->second;
		return true;
	}
	void PutGrid(const string &catalog_path, const CatalogGrid &g) {
		std::lock_guard<std::mutex> guard(lock);
		if (grids.emplace(catalog_path, g).second) {
			bytes += 256 * (g.cells.size() + 1);
			grown = true;
		}
	}

	bool GetCompanions(const string &prefix, CachedCompanions &out) const {
		std::lock_guard<std::mutex> guard(lock);
		auto it = companions.find(prefix);
		if (it == companions.end()) {
			return false;
		}
		out = it->second;
		return true;
	}
	void PutCompanions(const string &prefix, const CachedCompanions &c) {
		std::lock_guard<std::mutex> guard(lock);
		if (companions.emplace(prefix, c).second) {
			bytes += prefix.size() + c.pgen.size() + c.pvar.size() + c.psam.size() + c.eff_prefix.size();
			grown = true;
		}
	}

	//! `kind` names the loader ("full", "source", "info", "count"), so one path can hold several
	//! loads. `meta` is null for the loaders that produce only a SampleInfo.
	//! Cache tier 2, step A (docs/superpowers/specs/2026-11-20-metadata-cache-tiers-design.md; Teague,
	//! 2026-11-20): an entry keeps the IIDs, FIDs and SEX codes in a compact form, and a parquet
	//! source's column header. It never keeps a text .psam's rows: a load that needs them (a text
	//! "full" / "source" metadata load) is not cached, and the query re-parses the file. A Get rebuilds
	//! the SampleInfo (iid_to_idx is left for EnsureIidMap to build when a caller needs it).
	bool GetSamples(const string &kind, const string &path, SampleInfo &info, PfileSampleMetadata *meta) const {
		std::lock_guard<std::mutex> guard(lock);
		auto it = samples.find(kind + '\n' + path);
		if (it == samples.end()) {
			return false;
		}
		const auto &e = it->second;
		if (meta && !e.has_meta) {
			return false;
		}
		info = SampleInfo();
		e.iids.Expand(info.iids);
		e.fids.Expand(info.fids);
		info.sexes = e.sexes;
		info.sample_ct = e.sample_ct;
		if (meta) {
			*meta = e.meta;
		}
		return true;
	}
	void PutSamples(const string &kind, const string &path, const SampleInfo &info, const PfileSampleMetadata *meta) {
		if (meta && (!meta->from_parquet || !meta->rows.empty())) {
			return; // a text .psam's rows are never cached (re-parsed per query)
		}
		CachedSamples entry;
		if (!entry.iids.Assign(info.iids) || !entry.fids.Assign(info.fids)) {
			return; // over 4 GiB of IID text: not cached
		}
		entry.sexes = info.sexes;
		entry.sample_ct = info.sample_ct;
		idx_t sz = entry.iids.Bytes() + entry.fids.Bytes() + entry.sexes.size() + 64;
		if (meta) {
			entry.has_meta = true;
			entry.meta = *meta;
			for (auto &n : entry.meta.header.column_names) {
				sz += n.size() + 64;
			}
			for (auto &p : entry.meta.source_paths) {
				sz += p.size() + 32;
			}
		}
		std::lock_guard<std::mutex> guard(lock);
		if (samples.emplace(kind + '\n' + path, std::move(entry)).second) {
			bytes += sz;
			grown = true;
		}
	}

	//! True once after an insert: ObjectCache reads GetEstimatedCacheMemory only at Put, so the
	//! owner re-Puts the entry (RefreshAttachCacheWeight) to have the new size counted.
	bool TakeGrown() {
		std::lock_guard<std::mutex> guard(lock);
		bool g = grown;
		grown = false;
		return g;
	}

private:
	//! Strings as one buffer plus offsets (4 + length bytes each), in place of a vector<string>
	//! (32 bytes or more each).
	struct CompactStrings {
		string buf;
		vector<uint32_t> off; //! size n + 1 (empty when n == 0)
		bool Assign(const vector<string> &v) {
			buf.clear();
			off.clear();
			if (v.empty()) {
				return true;
			}
			idx_t total = 0;
			for (auto &x : v) {
				total += x.size();
			}
			if (total > std::numeric_limits<uint32_t>::max()) {
				return false;
			}
			buf.reserve(total);
			off.reserve(v.size() + 1);
			off.push_back(0);
			for (auto &x : v) {
				buf += x;
				off.push_back(static_cast<uint32_t>(buf.size()));
			}
			return true;
		}
		void Expand(vector<string> &out) const {
			out.clear();
			if (off.empty()) {
				return;
			}
			out.reserve(off.size() - 1);
			for (idx_t i = 0; i + 1 < off.size(); i++) {
				out.emplace_back(buf, off[i], off[i + 1] - off[i]);
			}
		}
		idx_t Bytes() const {
			return buf.size() + off.size() * sizeof(uint32_t);
		}
	};
	struct CachedSamples {
		CompactStrings iids;
		CompactStrings fids;
		vector<uint8_t> sexes;
		idx_t sample_ct = 0;
		bool has_meta = false;
		PfileSampleMetadata meta; //! parquet sources only (no rows)
	};
	mutable std::mutex lock;
	idx_t bytes = 0;
	bool grown = false;
	unordered_map<string, CatalogGrid> grids;
	unordered_map<string, CachedCompanions> companions;
	unordered_map<string, CachedSamples> samples;
};

} // namespace duckdb
