// SPDX-License-Identifier: GPL-3.0-or-later
//
// Shared fixtures for the split preamble-SIC core QA shards.
//
// The original single-TU qa_uwb_preamble_sic_core.cc (3,572 lines / 33 cases)
// was split into three executables to break its ~17 s compile and ~10 s run:
//
//   qa_uwb_preamble_sic_core_basics.cc  integer/fractional/CFO core cases
//   qa_uwb_preamble_sic_core_model.cc   S1 eight-symbol model guard cases
//   qa_uwb_preamble_sic_core_perf.cc    E0/R01 perf and cache-equivalence
//
// Every case and assertion was moved 1:1.

#ifndef OPENUWB_QA_PREAMBLE_SIC_TEST_HELPERS_H
#define OPENUWB_QA_PREAMBLE_SIC_TEST_HELPERS_H

#include <boost/test/unit_test.hpp>

#include <gnuradio/uwb/uwb_preamble_sic_core.h>

#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace {

using cplx = std::complex<float>;
using dcp = std::complex<double>;
namespace sic = openuwb::sic;

using sic::SicConfig;
using sic::SicMode;
using sic::SicModel;
using sic::SicResult;
using sic::SicScratch;
using sic::SicState;
using sic::SicSupportMode;
using sic::SicTemplates;
using sic::SicTruth;

constexpr double kPi = 3.14159265358979323846;
constexpr size_t kSps = 1016;
constexpr size_t kPre = 2000;
constexpr double kFs = 998.4e6;

// ---------------------------------------------------------------------------
// Synthetic model
// ---------------------------------------------------------------------------
double sinc_pi(double t)
{
    if (std::abs(t) < 1e-12)
        return 1.0;
    const double p = kPi * t;
    return std::sin(p) / p;
}

// Generator-side fractional placement model (Blackman-windowed sinc, H=24):
// deliberately different from the core's Hamming-windowed interpolator.
cplx gen_interp(const std::vector<cplx>& p, double x)
{
    constexpr int H = 24;
    const int64_t i0 = static_cast<int64_t>(std::floor(x));
    const double frac = x - std::floor(x);
    double ar = 0.0, ai = 0.0, wt = 0.0;
    for (int j = -H; j <= H; ++j) {
        const int64_t idx = i0 + j;
        if (idx < 0 || idx >= static_cast<int64_t>(p.size()))
            continue;
        const double t = static_cast<double>(j) - frac;
        const double w =
            0.42 + 0.5 * std::cos(kPi * t / H) + 0.08 * std::cos(2.0 * kPi * t / H);
        const double tap = sinc_pi(t) * w;
        ar += tap * p[static_cast<size_t>(idx)].real();
        ai += tap * p[static_cast<size_t>(idx)].imag();
        wt += tap;
    }
    if (std::abs(wt) > 1e-9) {
        ar /= wt;
        ai /= wt;
    }
    return cplx(static_cast<float>(ar), static_cast<float>(ai));
}

struct Model {
    size_t reps = 0;
    std::vector<cplx> preamble;
    std::vector<cplx> sync;
};

Model build_model(size_t code, size_t reps)
{
    Model m;
    m.reps = reps;
    const size_t cps = kSps / 2; // 2 samples/chip
    // Gaussian shaping pulse (length 9).
    const size_t T = 9;
    std::vector<double> pulse(T);
    double ps = 0.0;
    for (size_t i = 0; i < T; ++i) {
        const double d = static_cast<double>(i) - 4.0;
        pulse[i] = std::exp(-(d * d) / (2.0 * 1.2 * 1.2));
        ps += pulse[i];
    }
    for (double& v : pulse)
        v /= ps;
    // Deterministic +/-1 code.
    uint64_t s = 0x9E3779B97F4A7C15ull ^ (static_cast<uint64_t>(code) * 2654435761ull);
    std::vector<int> chips(cps);
    for (size_t i = 0; i < cps; ++i) {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        chips[i] = (s & 1u) ? 1 : -1;
    }
    // Upsampled chip stream, repeated per SYNC (same code each repetition).
    const size_t total_chips = reps * cps;
    std::vector<double> up(total_chips * 2, 0.0);
    for (size_t k = 0; k < total_chips; ++k)
        up[2 * k] = static_cast<double>(chips[k % cps]);
    const size_t plen = up.size() + T - 1;
    m.preamble.assign(plen, cplx(0.f, 0.f));
    for (size_t i = 0; i < up.size(); ++i) {
        if (up[i] == 0.0)
            continue;
        for (size_t j = 0; j < T; ++j)
            m.preamble[i + j] += cplx(static_cast<float>(up[i] * pulse[j]), 0.f);
    }
    m.sync.assign(m.preamble.begin(), m.preamble.begin() + kSps);
    return m;
}

struct Path {
    double tau;
    cplx g;
};

// Add one delayed + CFO-rotated copy of `pre` into acc (scaled by g).
void add_path(std::vector<cplx>& acc,
              const std::vector<cplx>& pre,
              int64_t base,
              double tau,
              cplx g,
              double cfo_hz)
{
    const int64_t span0 = base + static_cast<int64_t>(std::floor(tau)) - 32;
    const int64_t span1 = base + static_cast<int64_t>(std::ceil(tau)) +
                          static_cast<int64_t>(pre.size()) + 32;
    const int64_t lo = std::max<int64_t>(0, span0);
    const int64_t hi = std::min<int64_t>(static_cast<int64_t>(acc.size()), span1);
    const double w = 2.0 * kPi * cfo_hz / kFs;
    for (int64_t i = lo; i < hi; ++i) {
        const cplx v =
            gen_interp(pre, static_cast<double>(i) - static_cast<double>(base) - tau);
        const double ph = w * static_cast<double>(i);
        acc[static_cast<size_t>(i)] +=
            g * v *
            cplx(static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
    }
}

struct Signal {
    std::vector<cplx> y;
    std::vector<cplx> i;
};

Signal build_signal(const Model& m,
                    size_t n,
                    int64_t comm_start,
                    const std::vector<Path>& paths,
                    double cfo_hz,
                    double noise_std,
                    uint64_t seed)
{
    Signal s;
    s.i.assign(n, cplx(0.f, 0.f));
    for (const Path& p : paths)
        add_path(s.i, m.preamble, comm_start, p.tau, p.g, cfo_hz);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, noise_std > 0.0 ? noise_std : 1.0);
    s.y.assign(n, cplx(0.f, 0.f));
    for (size_t k = 0; k < n; ++k) {
        cplx w(0.f, 0.f);
        if (noise_std > 0.0)
            w = cplx(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
        s.y[k] = s.i[k] + w;
    }
    return s;
}

double nmse(const std::vector<cplx>& a, const std::vector<cplx>& b, size_t lo, size_t hi)
{
    double num = 0.0, den = 0.0;
    for (size_t k = lo; k < hi && k < a.size() && k < b.size(); ++k) {
        num += std::norm(a[k] - b[k]);
        den += std::norm(a[k]);
    }
    return (den > 0.0) ? num / den : 0.0;
}

double suppression_db(const std::vector<cplx>& a,
                      const std::vector<cplx>& b,
                      size_t lo,
                      size_t hi)
{
    const double e = nmse(a, b, lo, hi);
    return (e > 0.0) ? 10.0 * std::log10(1.0 / e) : 300.0;
}

// Realistic small AWGN keeps the holdout residual correlation well-defined.
constexpr double kQaNoise = 1e-3;

// ---------------------------------------------------------------------------
// Config helpers
// ---------------------------------------------------------------------------
struct Fixture {
    Model model;
    size_t n = 0;
    int64_t comm_start = 0;
    SicConfig cfg;
    SicScratch scratch;
    std::vector<cplx> out;
    std::vector<cplx> ihat;

    SicTemplates templates() const
    {
        SicTemplates t;
        t.sync = model.sync.data();
        t.sync_len = model.sync.size();
        t.preamble = model.preamble.data();
        t.preamble_len = model.preamble.size();
        return t;
    }
};

// Build a fixture: model of `reps` SYNC at code, window n, comm at comm_start,
// ROI around the comm start, and prepare the core.
Fixture make_fixture(size_t code,
                     size_t reps,
                     int64_t comm_start,
                     size_t roi_pad_lo,
                     size_t roi_pad_hi,
                     size_t k_paths,
                     size_t fit_reps,
                     size_t validation_reps,
                     size_t post_tail,
                     bool roi_wide = false)
{
    Fixture f;
    f.model = build_model(code, reps);
    f.n = kPre + reps * kSps + 3000;
    f.comm_start = comm_start;
    f.cfg.mode = SicMode::Estimate;
    f.cfg.comm_code = code;
    f.cfg.sync_reps = reps;
    f.cfg.samples_per_symbol = kSps;
    f.cfg.sample_rate = kFs;
    f.cfg.k_paths = k_paths;
    f.cfg.fit_reps = fit_reps;
    f.cfg.validation_reps = validation_reps;
    if (roi_wide) {
        f.cfg.roi_lo = 0;
        f.cfg.roi_hi = static_cast<int64_t>(f.n);
    } else {
        f.cfg.roi_lo = comm_start - static_cast<int64_t>(roi_pad_lo);
        f.cfg.roi_hi = comm_start + static_cast<int64_t>(reps * kSps) +
                       static_cast<int64_t>(roi_pad_hi);
    }
    f.cfg.coarse_stride = 14;
    f.cfg.path_min_spacing = 0;
    f.cfg.post_tail = post_tail;
    f.cfg.enable_fractional = true;
    f.cfg.enable_cfo = true;
    f.cfg.detect_threshold = 0.5;
    f.cfg.validate_threshold = 0.25;
    f.out.assign(f.n, cplx(0.f, 0.f));
    f.ihat.assign(f.n, cplx(0.f, 0.f));
    std::string err;
    const bool ok = sic::prepare_sic_core(f.cfg,
                                          f.model.sync.data(),
                                          f.model.sync.size(),
                                          f.model.preamble.data(),
                                          f.model.preamble.size(),
                                          f.scratch,
                                          &err);
    BOOST_CHECK_MESSAGE(ok, "prepare_sic_core failed: " << err);
    return f;
}

// Run one window and return the result (out/ihat live in the fixture).
SicResult run(Fixture& f, const Signal& s, const SicTruth* truth = nullptr)
{
    SicResult r;
    r.out = f.out.data();
    r.i_hat = f.ihat.data();
    const bool ok =
        sic::sic_one(s.y.data(), s.y.size(), f.templates(), f.cfg, f.scratch, r, truth);
    BOOST_CHECK_MESSAGE(
        ok, "sic_one returned false (state=" << static_cast<int>(r.state) << ")");
    return r;
}

const char* state_name(SicState s)
{
    switch (s) {
    case SicState::Bypassed:
        return "Bypassed";
    case SicState::Committed:
        return "Committed";
    case SicState::Failed:
        return "Failed";
    case SicState::Ambiguous:
        return "Ambiguous";
    }
    return "?";
}

// ---------------------------------------------------------------------------
// Minimal NPY / NPZ (zip, STORED) reader for the golden cross-check.
// ---------------------------------------------------------------------------
bool read_file(const std::string& path, std::vector<uint8_t>& bytes)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const std::streamoff sz = f.tellg();
    if (sz <= 0)
        return false;
    f.seekg(0, std::ios::beg);
    bytes.resize(static_cast<size_t>(sz));
    f.read(reinterpret_cast<char*>(bytes.data()), sz);
    return static_cast<bool>(f);
}

uint16_t le16(const uint8_t* p) { return static_cast<uint16_t>(p[0] | (p[1] << 8)); }
uint32_t le32(const uint8_t* p)
{
    return static_cast<uint32_t>(p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24));
}

// Extract one entry (STORED only) from an npz buffer by name.
bool npz_get(const std::vector<uint8_t>& z,
             const std::string& name,
             std::vector<uint8_t>& out)
{
    if (z.size() < 22)
        return false;
    // locate End Of Central Directory
    size_t eocd = std::string::npos;
    for (size_t i = z.size() - 22 + 1; i-- > 0;) {
        if (le32(&z[i]) == 0x06054b50u) {
            eocd = i;
            break;
        }
        if (z.size() - i > 66000)
            break;
    }
    if (eocd == std::string::npos)
        return false;
    const uint16_t nent = le16(&z[eocd + 10]);
    size_t cd = le32(&z[eocd + 16]);
    for (uint16_t e = 0; e < nent; ++e) {
        if (cd + 46 > z.size() || le32(&z[cd]) != 0x02014b50u)
            return false;
        const uint16_t method = le16(&z[cd + 10]);
        const uint32_t csize = le32(&z[cd + 20]);
        const uint16_t fnlen = le16(&z[cd + 28]);
        const uint16_t extlen = le16(&z[cd + 30]);
        const uint16_t cmlen = le16(&z[cd + 32]);
        const uint32_t lho = le32(&z[cd + 42]);
        const std::string fname(reinterpret_cast<const char*>(&z[cd + 46]), fnlen);
        if (fname == name) {
            if (method != 0)
                return false; // deflate: cannot inflate without zlib
            if (lho + 30 > z.size() || le32(&z[lho]) != 0x04034b50u)
                return false;
            const uint16_t lfn = le16(&z[lho + 26]);
            const uint16_t lex = le16(&z[lho + 28]);
            const size_t d0 = lho + 30 + lfn + lex;
            if (d0 + csize > z.size())
                return false;
            out.assign(z.begin() + d0, z.begin() + d0 + csize);
            return true;
        }
        cd += 46u + fnlen + extlen + cmlen;
    }
    return false;
}

// Parse a .npy buffer (complex64 / complex128, C order) into cplx.
bool npy_complex(const std::vector<uint8_t>& np, std::vector<cplx>& out)
{
    if (np.size() < 10 || np[0] != 0x93 || np[1] != 'N' || np[2] != 'U' || np[3] != 'M' ||
        np[4] != 'P' || np[5] != 'Y')
        return false;
    const uint8_t major = np[6];
    size_t hlen = 0, hoff = 0;
    if (major <= 1) {
        hlen = le16(&np[8]);
        hoff = 10;
    } else {
        hlen = le32(&np[8]);
        hoff = 12;
    }
    if (hoff + hlen > np.size())
        return false;
    const std::string h(reinterpret_cast<const char*>(&np[hoff]), hlen);
    if (h.find("'fortran_order': True") != std::string::npos)
        return false;
    const size_t sp = h.find("'shape': (");
    if (sp == std::string::npos)
        return false;
    const size_t ne = h.find(')', sp);
    if (ne == std::string::npos)
        return false;
    size_t count = 1;
    {
        std::string sh = h.substr(sp + 10, ne - sp - 10);
        std::string cur;
        for (char c : sh + ",") {
            if (c == ',' || c == ' ') {
                if (!cur.empty()) {
                    count *= static_cast<size_t>(std::stoul(cur));
                    cur.clear();
                }
            } else {
                cur.push_back(c);
            }
        }
    }
    const size_t d0 = hoff + hlen;
    if (h.find("'<c8'") != std::string::npos) {
        if (d0 + count * 8 > np.size())
            return false;
        out.resize(count);
        std::memcpy(out.data(), &np[d0], count * 8);
        return true;
    }
    if (h.find("'<c16'") != std::string::npos) {
        if (d0 + count * 16 > np.size())
            return false;
        out.resize(count);
        const double* d = reinterpret_cast<const double*>(&np[d0]);
        for (size_t i = 0; i < count; ++i)
            out[i] = cplx(static_cast<float>(d[2 * i]), static_cast<float>(d[2 * i + 1]));
        return true;
    }
    return false;
}

// Naive flat-JSON scalar lookup ("key": number) / bool.
bool json_number(const std::string& t, const std::string& key, double& v)
{
    const std::string pat = "\"" + key + "\"";
    size_t p = t.find(pat);
    if (p == std::string::npos)
        return false;
    p = t.find(':', p + pat.size());
    if (p == std::string::npos)
        return false;
    ++p;
    while (p < t.size() && (t[p] == ' ' || t[p] == '\t' || t[p] == '\n' || t[p] == '\r'))
        ++p;
    try {
        size_t used = 0;
        v = std::stod(t.substr(p), &used);
        return used > 0;
    } catch (...) {
        return false;
    }
}

// Naive flat-JSON string lookup ("key": "value").
bool json_string(const std::string& t, const std::string& key, std::string& v)
{
    const std::string pat = "\"" + key + "\"";
    size_t p = t.find(pat);
    if (p == std::string::npos)
        return false;
    p = t.find(':', p + pat.size());
    if (p == std::string::npos)
        return false;
    ++p;
    while (p < t.size() && (t[p] == ' ' || t[p] == '\t' || t[p] == '\n' || t[p] == '\r'))
        ++p;
    if (p >= t.size() || t[p] != '"')
        return false;
    ++p;
    const size_t e = t.find('"', p);
    if (e == std::string::npos)
        return false;
    v = t.substr(p, e - p);
    return true;
}

bool json_bool(const std::string& t, const std::string& key, bool& v)
{
    const std::string pat = "\"" + key + "\"";
    size_t p = t.find(pat);
    if (p == std::string::npos)
        return false;
    p = t.find(':', p + pat.size());
    if (p == std::string::npos)
        return false;
    const size_t tpos = t.find("true", p);
    const size_t fpos = t.find("false", p);
    if (tpos != std::string::npos && (fpos == std::string::npos || tpos < fpos)) {
        v = true;
        return true;
    }
    if (fpos != std::string::npos) {
        v = false;
        return true;
    }
    return false;
}

// Flat-JSON numeric array lookup ("key": [n, n, ...]) -> vector<double>.
// Used for the golden `expected.taus` cross-check (F1): the reference emits the
// selected fractional delays there, so the C++ selection can be asserted
// directly instead of only through an NMSE tolerance.
bool json_number_array(const std::string& t,
                       const std::string& key,
                       std::vector<double>& out)
{
    const std::string pat = "\"" + key + "\"";
    const size_t kp = t.find(pat);
    if (kp == std::string::npos)
        return false;
    const size_t lb = t.find('[', kp + pat.size());
    if (lb == std::string::npos)
        return false;
    const size_t rb = t.find(']', lb);
    if (rb == std::string::npos)
        return false;
    out.clear();
    std::string cur;
    for (size_t i = lb + 1; i <= rb; ++i) {
        const char c = (i == rb) ? ',' : t[i];
        if (c == ',') {
            size_t a = 0, b = cur.size();
            while (a < b &&
                   (cur[a] == ' ' || cur[a] == '\t' || cur[a] == '\n' || cur[a] == '\r'))
                ++a;
            while (b > a && (cur[b - 1] == ' ' || cur[b - 1] == '\t' ||
                             cur[b - 1] == '\n' || cur[b - 1] == '\r'))
                --b;
            if (b > a) {
                try {
                    out.push_back(std::stod(cur.substr(a, b - a)));
                } catch (...) {
                    return false;
                }
            }
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    return true;
}

std::string golden_root()
{
    const char* e = std::getenv("UWB_SIC_GOLDEN_DIR");
    if (e && *e)
        return e;
#ifdef UWB_TESTDATA_DIR
    return std::string(UWB_TESTDATA_DIR) + "/sic_preamble/golden";
#else
    return std::string();
#endif
}

bool dir_exists(const std::string& p)
{
    struct stat st;
    return ::stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

} // namespace

// ---------------------------------------------------------------------------
// S1 model-driver helpers (shared by the basics/model/perf shards).
// ---------------------------------------------------------------------------
// Compare EVERY frozen-model field (not just n_paths/valid) plus the terminal
// state and the field hash: the data-dependency tests must prove the whole
// model is independent of the tail, not just its size.
static void
check_models_identical(const sic::SicModel& a, const sic::SicModel& b, const char* what)
{
    BOOST_CHECK_MESSAGE(a.valid == b.valid, what << ": valid differs");
    BOOST_CHECK_MESSAGE(static_cast<int>(a.state) == static_cast<int>(b.state),
                        what << ": state differs");
    BOOST_CHECK_MESSAGE(a.model_id == b.model_id, what << ": model_id differs");
    BOOST_CHECK_MESSAGE(a.template_hash == b.template_hash,
                        what << ": template_hash differs");
    BOOST_CHECK_MESSAGE(a.config_hash == b.config_hash, what << ": config_hash differs");
    BOOST_CHECK_MESSAGE(a.origin == b.origin, what << ": origin differs");
    BOOST_CHECK_MESSAGE(a.cfo_hz == b.cfo_hz, what << ": cfo differs");
    BOOST_CHECK_MESSAGE(a.n_paths == b.n_paths, what << ": n_paths differs");
    BOOST_CHECK_MESSAGE(a.gain_condition == b.gain_condition,
                        what << ": gain_condition differs");
    BOOST_CHECK_MESSAGE(a.fit_residual == b.fit_residual,
                        what << ": fit_residual differs");
    BOOST_CHECK_MESSAGE(a.training_lo == b.training_lo && a.training_hi == b.training_hi,
                        what << ": training span differs");
    BOOST_CHECK_MESSAGE(a.observed_end_sample == b.observed_end_sample,
                        what << ": observed_end differs");
    BOOST_CHECK_MESSAGE(a.estimation_ready_sample == b.estimation_ready_sample,
                        what << ": estimation_ready differs");
    BOOST_CHECK_MESSAGE(a.predicted_support_lo == b.predicted_support_lo &&
                            a.predicted_support_hi == b.predicted_support_hi,
                        what << ": predicted support differs");
    BOOST_CHECK_MESSAGE(std::string(a.support_source) == std::string(b.support_source),
                        what << ": support_source differs");
    BOOST_CHECK_MESSAGE(a.sub_lo == b.sub_lo && a.sub_hi == b.sub_hi,
                        what << ": cancellation region differs");
    BOOST_CHECK_MESSAGE(a.training_reps_confirmed == b.training_reps_confirmed,
                        what << ": training_reps_confirmed differs");
    BOOST_CHECK_MESSAGE(a.nominal_sync_reps == b.nominal_sync_reps &&
                            a.samples_per_symbol == b.samples_per_symbol &&
                            a.sample_rate == b.sample_rate &&
                            a.phase_reference_sample == b.phase_reference_sample,
                        what << ": profile/phase reference differs");
    for (size_t k = 0; k < sic::kSicMaxPaths; ++k) {
        BOOST_CHECK_MESSAGE(a.tau_abs[k] == b.tau_abs[k], what << ": tau[" << k << "]");
        BOOST_CHECK_MESSAGE(a.gain[k] == b.gain[k], what << ": gain[" << k << "]");
    }
}

// Actual-output suppression (docs/sic/S1): 10*log10( energy(i_true) /
// energy(out - (y - i_true)) ) over [lo,hi).  This is the "real output"
// metric -- it compares the delivered window against the ideal clean signal,
// not against a reconstruction NMSE, so a bypass (out == y) scores 0 dB.
static double actual_suppression_db(const std::vector<cplx>& i_true,
                                    const std::vector<cplx>& y,
                                    const std::vector<cplx>& out,
                                    size_t lo,
                                    size_t hi)
{
    double ei = 0.0, er = 0.0;
    for (size_t k = lo; k < hi && k < y.size(); ++k) {
        ei += std::norm(i_true[k]);
        er += std::norm(out[k] - (y[k] - i_true[k]));
    }
    return (er > 0.0) ? 10.0 * std::log10(ei / er) : 300.0;
}

static double region_energy(const std::vector<cplx>& v, size_t lo, size_t hi)
{
    double e = 0.0;
    for (size_t k = lo; k < hi && k < v.size(); ++k)
        e += std::norm(v[k]);
    return e;
}

static SicModel estimate_model_checked(Fixture& f,
                                       const Signal& s,
                                       int64_t oe,
                                       SicResult* r_out = nullptr)
{
    SicModel m;
    SicResult r;
    const bool ok = sic::estimate_model(
        s.y.data(), s.y.size(), oe, f.templates(), f.cfg, f.scratch, m, nullptr, &r);
    BOOST_REQUIRE_MESSAGE(ok, "estimate_model returned false");
    if (r_out)
        *r_out = r;
    return m;
}

// Apply the frozen model to a tiling of [0, n) given as (offset, length) pairs.
static void apply_chunks(Fixture& f,
                         const SicModel& m,
                         const Signal& s,
                         const std::vector<std::pair<int64_t, size_t>>& chunks,
                         std::vector<cplx>& out,
                         std::vector<cplx>& ihat)
{
    out.assign(s.y.size(), cplx(0.f, 0.f));
    ihat.assign(s.y.size(), cplx(0.f, 0.f));
    for (const auto& c : chunks) {
        BOOST_REQUIRE(c.first >= 0 &&
                      static_cast<size_t>(c.first) + c.second <= s.y.size());
        const bool ok = sic::apply_model(m,
                                         s.y.data() + c.first,
                                         c.second,
                                         c.first,
                                         f.templates(),
                                         f.cfg,
                                         f.scratch,
                                         out.data() + c.first,
                                         ihat.data() + c.first);
        BOOST_REQUIRE(ok);
    }
}

#endif // OPENUWB_QA_PREAMBLE_SIC_TEST_HELPERS_H
