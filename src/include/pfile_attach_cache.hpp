#pragma once

#include "duckdb.hpp"
#include "duckdb/storage/object_cache.hpp"
#include "plink_common.hpp"
#include "pfile_sample_metadata.hpp"

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
		}
	}

	//! `kind` names the loader ("full", "source", "info", "count"), so one path can hold several
	//! loads. `meta` is null for the loaders that produce only a SampleInfo.
	bool GetSamples(const string &kind, const string &path, SampleInfo &info, PfileSampleMetadata *meta) const {
		std::lock_guard<std::mutex> guard(lock);
		auto it = samples.find(kind + '\n' + path);
		if (it == samples.end()) {
			return false;
		}
		info = it->second.info;
		if (meta) {
			*meta = it->second.meta;
		}
		return true;
	}
	void PutSamples(const string &kind, const string &path, const SampleInfo &info, const PfileSampleMetadata *meta) {
		std::lock_guard<std::mutex> guard(lock);
		CachedSamples entry;
		entry.info = info;
		if (meta) {
			entry.meta = *meta;
		}
		idx_t sz = 64 * (info.iids.size() + info.fids.size() + 1) + info.sexes.size();
		for (auto &row : entry.meta.rows) {
			for (auto &f : row) {
				sz += f.size() + 32;
			}
		}
		if (samples.emplace(kind + '\n' + path, std::move(entry)).second) {
			bytes += sz;
		}
	}

private:
	struct CachedSamples {
		SampleInfo info;
		PfileSampleMetadata meta;
	};
	mutable std::mutex lock;
	idx_t bytes = 0;
	unordered_map<string, CatalogGrid> grids;
	unordered_map<string, CachedCompanions> companions;
	unordered_map<string, CachedSamples> samples;
};

} // namespace duckdb
