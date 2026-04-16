#pragma once

// Bind-phase timing probes, gated on the PLINKING_BIND_PROFILE env var.
// Header-only to avoid touching the build. Emits to stderr so output is
// captured even when an exception unwinds the bind call.

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace duckdb {

inline bool BindProfileEnabled() {
	static int cached = -1;
	if (cached == -1) {
		const char *v = std::getenv("PLINKING_BIND_PROFILE");
		cached = (v && v[0] && std::strcmp(v, "0") != 0) ? 1 : 0;
	}
	return cached == 1;
}

struct BindPhaseTimer {
	using Clock = std::chrono::steady_clock;
	std::string label;
	Clock::time_point start;
	bool active;

	explicit BindPhaseTimer(std::string lbl) : label(std::move(lbl)), active(BindProfileEnabled()) {
		if (active) {
			start = Clock::now();
			std::fprintf(stderr, "[BIND_PROFILE] ENTER %s\n", label.c_str());
			std::fflush(stderr);
		}
	}

	void Note(const char *fmt, ...) {
		if (!active) {
			return;
		}
		auto now = Clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
		std::fprintf(stderr, "[BIND_PROFILE]   %s @ %lldms: ", label.c_str(), static_cast<long long>(ms));
		va_list args;
		va_start(args, fmt);
		std::vfprintf(stderr, fmt, args);
		va_end(args);
		std::fprintf(stderr, "\n");
		std::fflush(stderr);
	}

	~BindPhaseTimer() {
		if (active) {
			auto now = Clock::now();
			auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
			std::fprintf(stderr, "[BIND_PROFILE] LEAVE %s: %lldms\n", label.c_str(), static_cast<long long>(ms));
			std::fflush(stderr);
		}
	}
};

} // namespace duckdb
