/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * TWR configuration, split by responsibility (Phase 2).
 * The aggregate entry point is <openuwb/twr/twr_config.h>; these pieces are
 * ordered, not independent: each one includes its predecessor.
 */
#ifndef INCLUDED_OPENUWB_TWR_CONFIG_VALIDATOR_H
#define INCLUDED_OPENUWB_TWR_CONFIG_VALIDATOR_H

#include <openuwb/twr/twr_config_fields.h>
namespace openuwb {
namespace twr {

// ===========================================================================
// 6. The validator
// ===========================================================================

// Every check appends a specific violation.  There is no early return: a
// config that is wrong in many places is reported in full so the operator
// fixes them in one pass.
class ConfigValidator
{
public:
    explicit ConfigValidator(const Capabilities& caps) : d_caps(caps) {}

    ValidationReport run(const TwrConfig& cfg) const;

private:
    // -- enum domains (N07) ----------------------------------------------
    //
    // One place that refuses an enum value the enum does not have, for EVERY
    // enum field in the config.  Fourteen of these were verified to be
    // accepted outright before this pass existed, and a further few were
    // refused only incidentally, with the violation naming an unrelated field.
    //
    // The check is deliberately unconditional: it does not matter whether the
    // value came from a cast, a deserialiser, an async-event mapping or a
    // future enum member whose gate was not updated.  `UnknownEnumValue` is the
    // same reason the JSON reader uses, so a caller sees one vocabulary
    // whichever path produced the bad value.
    void check_enum(const std::string& f, bool known, const char* what) const
    {
        if (!known)
            cfg_rej(f, ConfigReason::UnknownEnumValue,
                    std::string("value is not a member of ") + what +
                        "; an out-of-domain enum is refused rather than reasoned "
                        "about (it would serialise as \"invalid\" and be refused "
                        "by the reader, so accepting it would contradict this "
                        "layer's own contract)",
                    "REQ-API-01");
    }

    void check_enum_domains(const TwrConfig& c) const;

    void check_timed_enums(const std::string& f, const TimedField& t) const;

    // -- helpers ---------------------------------------------------------
    void rej(const std::string& f,
             ConfigReason r,
             ExchangeStatus s,
             const std::string& m,
             const char* req) const
    {
        d_report.add(f, r, s, m, req);
    }
    void cfg_rej(const std::string& f,
                 ConfigReason r,
                 const std::string& m,
                 const char* req = "REQ-API-01") const
    {
        rej(f, r, ExchangeStatus::ConfigRejected, m, req);
    }
    void unsup(const std::string& f,
               const std::string& m,
               ConfigReason r = ConfigReason::Unsupported,
               const char* req = "REQ-SCOPE-01") const
    {
        rej(f, r, ExchangeStatus::Unsupported, m, req);
    }
    // A float field that must be finite.
    void finite(const std::string& f, double v, const char* req = "REQ-API-01") const
    {
        if (!std::isfinite(v))
            cfg_rej(f, ConfigReason::NotFinite, "must be finite (NaN/Inf rejected)",
                req);
    }
    void check_finite_opt(const std::string& f, const Opt<double>& o) const
    {
        if (o.has() && !std::isfinite(o.value()))
            cfg_rej(f, ConfigReason::NotFinite,
                "must be finite when present (NaN/Inf rejected)");
    }
    void check_finite_timed(const std::string& f, const TimedField& t) const
    {
        finite(f + ".quantisation_hz", t.required_quantisation_hz);
        if (t.value.negative())
            cfg_rej(f, ConfigReason::NegativeValue,
                "a delay or window may not be negative");
    }

    // -- sections --------------------------------------------------------
    void check_meta(const TwrConfig& c) const;

    void check_session(const TwrConfig& c) const;

    void check_phy(const TwrConfig& c) const;

    void check_frame(const TwrConfig& c) const;

    void check_tx(const TwrConfig& c) const;

    void check_rx(const TwrConfig& c) const;

    void check_radio(const TwrConfig& c) const;

    // The generic timed-field contract.  Returns nothing; appends violations.
    void validate_timed_field(const std::string& f,
                              const TimedField& t,
                              bool required,
                              TimeReferenceEvent expect_ref,
                              bool require_marker) const
    {
        if (t.domain == TimeDomain::Unspecified) {
            cfg_rej(f + ".domain",
                    ConfigReason::EmptyValue,
                    "the clock domain must be stated (device_ticks or monotonic_host)",
                    "REQ-API-02");
            return;
        }
        if (t.value.negative())
            cfg_rej(f + ".ns", ConfigReason::NegativeValue, "may not be negative",
                "REQ-API-02");
        if (required && t.value.is_zero())
            cfg_rej(f + ".ns", ConfigReason::ZeroValue, "must be > 0", "REQ-API-02");
        if (expect_ref != TimeReferenceEvent::HostMonotonic && t.reference != expect_ref)
            cfg_rej(f + ".reference",
                    ConfigReason::FieldConflict,
                    std::string("must be ") + time_reference_event_to_string(expect_ref) +
                        std::string(", is ") + time_reference_event_to_string(t.reference),
                    "REQ-API-02");
        if (t.domain == TimeDomain::DeviceTicks) {
            if (require_marker && !t.marker.has())
                cfg_rej(f + ".marker",
                        ConfigReason::EmptyValue,
                        "a device-tick field must name the RF marker it is measured at "
                        "(REQ-TIME-02)",
                        "REQ-API-02");
            if (t.marker.has() && t.marker.value() == TimestampMarker::UhdRxFirstIqSample &&
                expect_ref == TimeReferenceEvent::FrameTail)
                cfg_rej(f + ".marker",
                        ConfigReason::FieldConflict,
                        "the first RX IQ sample is not the end of a frame",
                        "REQ-TIME-02");
        } else { // MonotonicHost
            if (t.marker.has())
                cfg_rej(f + ".marker",
                        ConfigReason::FieldConflict,
                        "a host monotonic field has no RF marker; a marker here would invite a "
                        "cross-domain subtraction",
                        "REQ-TIME-01");
            if (t.required_quantisation_hz != 1.0e9)
                cfg_rej(f + ".quantisation_hz",
                        ConfigReason::QuantisationUnsupported,
                        "a host monotonic field must be quantised at 1e9 (integer nanoseconds)",
                        "REQ-TIME-01");
        }
        if (!std::isfinite(t.required_quantisation_hz) || t.required_quantisation_hz <= 0.0) {
            cfg_rej(f + ".quantisation_hz",
                    ConfigReason::QuantisationUnsupported,
                    "the required quantisation rate must be finite and > 0",
                    "REQ-API-01");
            return;
        }
        if (!d_caps.quantisation_supported(t.required_quantisation_hz)) {
            cfg_rej(f + ".quantisation_hz",
                    ConfigReason::QuantisationUnsupported,
                    "quantisation " + twr_double_to_text(t.required_quantisation_hz) +
                        " Hz is not a rate this build can express (allowed: 1e9 ns, " +
                        twr_double_to_text(d_caps.work_sample_rate_hz) + " work rate, the native "
                        "device rates)",
                    "REQ-API-01");
            return;
        }
        if (t.max_quantisation_error_ns < 0)
            cfg_rej(f + ".max_quantisation_error_ns",
                    ConfigReason::OutOfRange,
                    "the tolerated quantisation error may not be negative");
        // The two distinct precision gates.
        if (!t.value.representable_at(t.required_quantisation_hz,
            t.max_quantisation_error_ns))
            cfg_rej(f + ".ns",
                    ConfigReason::TimingPrecisionInsufficient,
                    "a device tick at " + twr_double_to_text(t.required_quantisation_hz) +
                        " Hz cannot express this value with <= " +
                        twr_int_to_text(t.max_quantisation_error_ns) + " ns of error",
                    "REQ-API-01");
        bool ok = false;
        (void)quantise_duration(t.value, t.required_quantisation_hz, ok);
        if (!ok)
            cfg_rej(f + ".ns",
                    ConfigReason::TimingPrecisionInsufficient,
                    "the value does not fit the device tick range at " +
                        twr_double_to_text(t.required_quantisation_hz) + " Hz",
                    "REQ-API-01");
    }

    void check_timing(const TwrConfig& c) const;

    void check_timeouts(const TwrConfig& c) const;

    void check_calibration(const TwrConfig& c) const;

    void check_diagnostics(const TwrConfig& c) const;

    // Frame lengths: Poll / Response / Final must each fit the negotiated PSDU
    // (REQ-API-01 帧格式; REQ-PHY-01 standard-length PSDU).
    void check_frame_lengths(const TwrConfig& c) const;

    const Capabilities& d_caps;
    mutable ValidationReport d_report;
};

inline ValidationReport validate(const TwrConfig& cfg, const Capabilities& caps)
{
    ConfigValidator v(caps);
    return v.run(cfg);
}

inline ValidationReport validate(const TwrConfig& cfg) { return validate(cfg,
    capabilities()); }

} // namespace twr
} // namespace openuwb

#endif /* INCLUDED_OPENUWB_TWR_CONFIG_VALIDATOR_H */
