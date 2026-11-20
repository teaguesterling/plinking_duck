#include "pgen_header_cache.hpp"
#include "pgen_vfs_opener.hpp"

namespace duckdb {

std::atomic<int64_t> g_pgen_header_decodes {0};

static shared_ptr<PgenHeader> DecodePgenHeader(ClientContext &context, bool use_vfs, const string &pgen_path,
                                               uint32_t variant_ct, uint32_t sample_ct, const string &what) {
	// No PgenVfsScope here: the caller holds one, and a nested scope's destructor clears the global
	// opener (pgen_vfs_opener.cpp), so the caller's next open (PgrInit) would bypass the VFS.
	auto h = make_shared_ptr<PgenHeader>();
	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;
	plink2::PglErr err = plink2::PgfiInitPhase1(pgen_path.c_str(), nullptr, variant_ct, sample_ct, &header_ctrl,
	                                            &h->pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		throw IOException("%s (phase 1) for '%s': %s", what, pgen_path, errstr_buf);
	}
	if (pgfi_alloc_cacheline_ct > 0) {
		h->alloc.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}
	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, h->pgfi.raw_variant_ct, &h->max_vrec_width, &h->pgfi,
	                             h->alloc.As<unsigned char>(), &h->pgr_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		throw IOException("%s (phase 2) for '%s': %s", what, pgen_path, errstr_buf);
	}
	// Every reader opens its own handle, so no PgrInit writes to the shared struct.
	if (h->pgfi.shared_ff) {
		fclose(h->pgfi.shared_ff);
		h->pgfi.shared_ff = nullptr;
	}
	h->bytes = pgfi_alloc_cacheline_ct * plink2::kCacheline + sizeof(PgenHeader);
	g_pgen_header_decodes.fetch_add(1, std::memory_order_relaxed);
	return h;
}

shared_ptr<PgenHeader> AcquirePgenHeader(ClientContext &context, bool use_vfs, const string &pgen_path,
                                         uint32_t variant_ct, uint32_t sample_ct, const string &what) {
	return DecodePgenHeader(context, use_vfs, pgen_path, variant_ct, sample_ct, what);
}

void InitPgenReaderOn(const PgenHeader &header, const string &pgen_path, plink2::PgenReader &pgr,
                      AlignedBuffer &pgr_alloc, const string &what) {
	if (header.pgr_alloc_cacheline_ct > 0) {
		pgr_alloc.Allocate(header.pgr_alloc_cacheline_ct * plink2::kCacheline);
	} else {
		pgr_alloc.Reset();
	}
	plink2::PreinitPgr(&pgr);
	// PgrInit takes a non-const pointer; with shared_ff closed it only reads the struct (it copies it
	// into pgr.fi and opens pgen_path itself).
	auto err = plink2::PgrInit(pgen_path.c_str(), header.max_vrec_width,
	                           const_cast<plink2::PgenFileInfo *>(&header.pgfi), &pgr, pgr_alloc.As<unsigned char>());
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr ce = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&pgr, &ce);
		throw IOException("%s: PgrInit failed for '%s'", what, pgen_path);
	}
}

} // namespace duckdb
