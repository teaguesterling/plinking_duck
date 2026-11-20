#pragma once

#include "duckdb.hpp"
#include "duckdb/storage/object_cache.hpp"
#include "plink_common.hpp"

#include <atomic>

namespace duckdb {

//! A decoded .pgen header: pgenlib's PgenFileInfo after PgfiInitPhase1 + 2, and the allocation its
//! per-variant arrays (vrtypes, var_fpos: 9 B per variant) live in. Its shared_ff is closed, so a
//! PgrInit on it only reads it (pgenlib mode 2: each reader opens its own FILE*). A PgenReader copies
//! the struct but points into `alloc`, so a reader must hold the shared_ptr while it is open.
//! Cache tier 1 (docs/superpowers/specs/2026-11-20-metadata-cache-tiers-design.md).
class PgenHeader : public ObjectCacheEntry {
public:
	plink2::PgenFileInfo pgfi;
	AlignedBuffer alloc;
	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;
	idx_t bytes = 0;

	PgenHeader() {
		plink2::PreinitPgfi(&pgfi);
	}
	~PgenHeader() override {
		plink2::PglErr ce = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &ce);
	}
	static string ObjectType() {
		return "plinking_pgen_header";
	}
	string GetObjectType() override {
		return ObjectType();
	}
	optional_idx GetEstimatedCacheMemory() const override {
		return optional_idx(bytes);
	}
};

//! The decoded header of `pgen_path`. `variant_ct` / `sample_ct` are the counts the caller expects
//! (UINT32_MAX = unknown), as PgfiInitPhase1 takes them. `what` names the file in errors. The caller
//! must already hold PgenVfsScope(use_vfs) (a nested scope would clear the opener when it closes).
shared_ptr<PgenHeader> AcquirePgenHeader(ClientContext &context, bool use_vfs, const string &pgen_path,
                                         uint32_t variant_ct, uint32_t sample_ct, const string &what);

//! Open `pgr` on `header` (PgrInit, which opens its own FILE* through plinking_pgen_fopen; the caller
//! holds any PgenVfsScope). `pgr_alloc` receives the reader's buffer. Throws IOException naming `what`.
void InitPgenReaderOn(const PgenHeader &header, const string &pgen_path, plink2::PgenReader &pgr,
                      AlignedBuffer &pgr_alloc, const string &what);

//! Test-only: header decodes (PgfiInitPhase2 runs) since the last reset.
extern std::atomic<int64_t> g_pgen_header_decodes;

} // namespace duckdb
