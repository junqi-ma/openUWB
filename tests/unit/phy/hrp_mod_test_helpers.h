#pragma once
// Shared testdata loaders for HRP modulator QA suites.

#include <boost/test/unit_test.hpp>

#include <gnuradio/uwb/uwb_gr_complex.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace openuwb::hrp_mod_test {

#ifdef UWB_TESTDATA_DIR
inline constexpr const char* kTestdata = UWB_TESTDATA_DIR;
#else
inline constexpr const char* kTestdata = "../../../testdata";
#endif

inline std::string golden_dir()
{
    return std::string(kTestdata) + "/realtime_demod_golden";
}

inline bool load_cf32(const std::string& path, std::vector<gr_complex>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const auto bytes = static_cast<size_t>(f.tellg());
    f.seekg(0);
    if (bytes == 0 || bytes % sizeof(gr_complex) != 0)
        return false;
    out.resize(bytes / sizeof(gr_complex));
    f.read(reinterpret_cast<char*>(out.data()),
           static_cast<std::streamsize>(bytes));
    return static_cast<bool>(f);
}

inline bool load_u8(const std::string& path, std::vector<uint8_t>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const auto bytes = static_cast<size_t>(f.tellg());
    f.seekg(0);
    out.resize(bytes);
    f.read(reinterpret_cast<char*>(out.data()),
           static_cast<std::streamsize>(bytes));
    return static_cast<bool>(f);
}

inline bool load_i8(const std::string& path, std::vector<int8_t>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const auto bytes = static_cast<size_t>(f.tellg());
    f.seekg(0);
    out.resize(bytes);
    f.read(reinterpret_cast<char*>(out.data()),
           static_cast<std::streamsize>(bytes));
    return static_cast<bool>(f);
}

inline bool load_f32(const std::string& path, std::vector<float>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const auto bytes = static_cast<size_t>(f.tellg());
    f.seekg(0);
    if (bytes == 0 || bytes % sizeof(float) != 0)
        return false;
    out.resize(bytes / sizeof(float));
    f.read(reinterpret_cast<char*>(out.data()),
           static_cast<std::streamsize>(bytes));
    return static_cast<bool>(f);
}

inline std::vector<gr_complex> load_template()
{
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(
        load_cf32(std::string(kTestdata) + "/reference_preamble.bin", tmpl));
    return tmpl;
}

inline double l2_rel(const std::vector<gr_complex>& a, const gr_complex* b,
                     size_t n)
{
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double dr = static_cast<double>(a[i].real()) - b[i].real();
        const double di = static_cast<double>(a[i].imag()) - b[i].imag();
        num += dr * dr + di * di;
        den += static_cast<double>(b[i].real()) * b[i].real() +
               static_cast<double>(b[i].imag()) * b[i].imag();
    }
    return std::sqrt(num / (den + 1e-30));
}

} // namespace openuwb::hrp_mod_test
