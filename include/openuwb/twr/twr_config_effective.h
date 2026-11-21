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
#ifndef INCLUDED_OPENUWB_TWR_CONFIG_EFFECTIVE_H
#define INCLUDED_OPENUWB_TWR_CONFIG_EFFECTIVE_H

#include <openuwb/twr/twr_config_validator.h>
namespace openuwb {
namespace twr {

// ===========================================================================
// 7. The frozen effective config
// ===========================================================================

// A frozen, self-describing snapshot (REQ-API-01: 完整保存 requested、effective、
// 硬件 readback).  It is a VALUE: an exchange in flight holds a copy, so a
// later change to the requested config cannot reach it (REQ-API-03).
struct EffectiveConfig {
    bool ok = false;
    ValidationReport validation;

    std::string schema_version;
    std::string profile_version;
    std::string calibration_version;
    std::string config_hash;

    TwrConfig requested;
    TwrConfig effective;
    std::vector<FieldChange> changes; // requested vs effective, per field

    // The hardware readback exactly as recorded at startup.
    RadioReadback readback;

    // Quantised timing actually handed to the device.  The REQUESTED delay and
    // the EFFECTIVE (tick-quantised) delay are both available so REQ-GR-04 can
    // report the margin.
    int64_t poll_to_response_ticks = 0;
    int64_t response_to_final_ticks = 0;
    int64_t post_tx_rx_enable_ticks = 0;
    int64_t poll_start_ticks = 0;
    Duration poll_to_response_effective{0};
    Duration response_to_final_effective{0};
    Duration post_tx_rx_enable_effective{0};
    Duration poll_start_effective{0};
    double tick_rate_hz = 0.0;

    // Frame budget, per frame type, as transmitted (including the FCS exactly
    // once).  These are the CODEC's numbers, not a re-summation of the
    // config's geometry claim: see effective_config().
    uint16_t poll_bytes = 0;
    uint16_t response_bytes = 0;
    uint16_t final_bytes = 0;
    uint16_t max_psdu_bytes = 0;
    uint32_t max_timestamp_count = 0;
    // The PHR of this profile, recorded so no consumer recomputes it: the
    // number of INFORMATION bytes and the number of transmitted coded bits.
    uint16_t phr_bytes = 0;
    uint16_t phr_coded_bits = 0;
};

// Build the frozen snapshot.  On failure `ok == false`, `validation` holds
// every reason, and `effective` is a value-initialised (NOT a default
// substituted) config: there is no way to accidentally use it.
inline EffectiveConfig
effective_config(const TwrConfig& cfg, const Capabilities& caps)
{
    EffectiveConfig out;
    out.requested = cfg;
    out.validation = validate(cfg, caps);
    out.ok = out.validation.ok();
    out.schema_version = caps.schema_version;
    out.profile_version = cfg.meta.profile_version;
    out.calibration_version = cfg.meta.calibration_version;
    if (!out.ok) {
        // A rejected config produces no effective values at all.  The hash of
        // the REQUESTED config is still recorded so a failure can be tied to
        // the exact input that caused it (REQ-OUT-01).
        out.effective = TwrConfig{};
        out.config_hash = config_hash(cfg);
        out.changes.clear();
        return out;
    }

    // Effective == requested plus the materialisations the validator proved
    // legal: nothing is invented, everything is recorded.
    out.effective = cfg;
    out.config_hash = config_hash(out.effective);
    out.readback = cfg.radio.readback;
    out.changes = diff_fields(out.requested, out.effective);
    out.max_psdu_bytes = caps.max_psdu_bytes();
    out.max_timestamp_count = frame_type_timestamp_count(FrameType::Final);

    const double rate = caps.native_rate_supported(cfg.radio.native_sample_rate_hz)
                            ? cfg.radio.native_sample_rate_hz
                            : caps.native_rates_hz.front();
    out.tick_rate_hz = rate;

    struct Q {
        const TimedField* in;
        int64_t* ticks;
        Duration* eff;
    } q[] = {
        { &cfg.timing.poll_start, &out.poll_start_ticks, &out.poll_start_effective },
        { &cfg.timing.poll_to_response, &out.poll_to_response_ticks,
          &out.poll_to_response_effective },
        { &cfg.timing.response_to_final, &out.response_to_final_ticks,
          &out.response_to_final_effective },
        { &cfg.timing.post_tx_rx_enable, &out.post_tx_rx_enable_ticks,
          &out.post_tx_rx_enable_effective },
    };
    for (const auto& item : q) {
        if (item.in->domain != TimeDomain::DeviceTicks)
            continue;
        bool ok = false;
        const int64_t ticks =
            quantise_duration(item.in->value, item.in->required_quantisation_hz, ok);
        *item.ticks = ok ? ticks : 0;
        Duration back;
        *item.eff = (ok && Duration::from_ticks(ticks, item.in->required_quantisation_hz, back))
                        ? back
                        : Duration{0};
        if (back != item.in->value) {
            // Record the quantisation as an explicit requested/effective
            // difference; the device gets the tick value, the operator sees both.
            FieldChange ch;
            ch.path = item.in->value.nanos() != 0 ? "timing.quantised_delay" : "timing.delay";
            ch.requested = twr_duration_to_text(item.in->value);
            ch.effective = twr_duration_to_text(back);
            ch.note = "quantised at " + twr_double_to_text(item.in->required_quantisation_hz) +
                      " Hz (ticks=" + twr_int_to_text(ticks) + ")";
            out.changes.push_back(ch);
        }
    }

    // The frame budget is the CODEC's, computed once through the authority.
    // M0 re-summed the config's own geometry claim here, which is how a 7-byte
    // header produced a 9-byte Poll and understated the on-air length; and it
    // added the FCS from a raw config field, so a wrong fcs_bytes silently
    // changed the answer.  `on_air_bytes()` is invariant under which layer
    // appends the FCS, which is exactly the property the budget should have.
    FrameProfileGeometry authority;
    FrameProfile profile;
    if (frame_geometry_for(cfg.frame.frame_profile, authority) &&
        frame_profile_for(cfg.frame.frame_profile, profile)) {
        out.poll_bytes = static_cast<uint16_t>(authority.on_air_bytes(FrameType::Poll, profile));
        out.response_bytes =
            static_cast<uint16_t>(authority.on_air_bytes(FrameType::Response, profile));
        out.final_bytes =
            static_cast<uint16_t>(authority.on_air_bytes(FrameType::Final, profile));
        // The PHR is present in every frame of this profile; recording its
        // information size here means a consumer never has to recompute it.
        out.phr_bytes = static_cast<uint16_t>(authority.phr_bytes());
        out.phr_coded_bits = authority.phr_coded_bits();
    }
    return out;
}

inline EffectiveConfig effective_config(const TwrConfig& cfg)
{
    return effective_config(cfg, capabilities());
}

// ---------------------------------------------------------------------------
// In-flight gate (phase 1: at most one exchange per endpoint)
// ---------------------------------------------------------------------------

// Enforces "同端初版一次仅一个在途 exchange" at RUN time, next to the config
// that promised it.  A second concurrent request is QueueFull, never a
// silently dropped or queued-forever exchange.
class ExchangeGate
{
public:
    explicit ExchangeGate(uint32_t max_in_flight = kTwrMaxInFlightExchanges)
        : d_max(max_in_flight == 0 ? 1u : max_in_flight)
    {
    }

    ValidationReport try_begin()
    {
        ValidationReport r;
        if (d_in_flight >= d_max) {
            ConfigViolation v;
            v.field = "session.max_in_flight_exchanges";
            v.reason = ConfigReason::ExchangeAlreadyInFlight;
            v.status = ExchangeStatus::QueueFull;
            v.message = "an exchange is already in flight; at most " + twr_int_to_text(d_max) +
                        " may be outstanding";
            v.requirement = "REQ-API-01";
            r.add(v);
            // A REFUSED request must not consume a slot, or a run of refusals
            // would deadlock the endpoint.
            return r;
        }
        ++d_in_flight;
        return r;
    }
    void end()
    {
        if (d_in_flight > 0)
            --d_in_flight;
    }
    uint32_t in_flight() const { return d_in_flight; }

private:
    uint32_t d_max;
    uint32_t d_in_flight = 0;
};

// ---------------------------------------------------------------------------
// REQ-API-03: per-message / per-exchange overrides on an immutable snapshot
// ---------------------------------------------------------------------------

// Only the values a single exchange may re-negotiate.  Anything that affects
// synchronisation (PHY, frequency, gain, addresses, calibration) is
// deliberately absent: it may only change at a session boundary, and the
// snapshot an exchange already holds is immutable (REQ-API-03).
struct PerMessageOverrides {
    // 0 = "not overridden"; the validator then requires the base config value.
    int64_t measurement_count = 0;
    int64_t measurement_interval_ns = 0;
    int64_t poll_to_response_ns = 0;
    int64_t response_to_final_ns = 0;
    int64_t post_tx_rx_enable_ns = 0;
    int64_t rx_timeout_ns = 0;
    int64_t exchange_timeout_ns = 0;

    bool any() const
    {
        return measurement_count != 0 || measurement_interval_ns != 0 ||
               poll_to_response_ns != 0 || response_to_final_ns != 0 ||
               post_tx_rx_enable_ns != 0 || rx_timeout_ns != 0 ||
                   exchange_timeout_ns != 0;
    }
};

// The frozen snapshot plus the exchange gate.  Copying it is the only way to
// hand an immutable config to a worker; there is no mutable accessor.
class TwrConfigSnapshot
{
public:
    TwrConfigSnapshot() : d_effective(), d_caps(&capabilities()) {}
    TwrConfigSnapshot(const EffectiveConfig& e, const Capabilities& caps)
        : d_effective(e), d_caps(&caps)
    {
    }

    const EffectiveConfig& get() const { return d_effective; }
    bool ok() const { return d_effective.ok; }

    uint32_t in_flight() const { return d_gate.in_flight(); }

    // Begin an exchange.  Returns Ok when the gate admitted it.
    ValidationReport begin_exchange() { return d_gate.try_begin(); }
    void end_exchange() { d_gate.end(); }

    // REQ-API-03: overrides take effect at an exchange BOUNDARY.  A call
    // while an exchange is in flight is rejected -- the running exchange keeps
    // the snapshot it started with.
    ValidationReport
    set_overrides(const PerMessageOverrides& o, TwrConfigSnapshot& out) const
    {
        ValidationReport r;
        if (d_gate.in_flight() > 0) {
            ConfigViolation v;
            v.field = "session";
            v.reason = ConfigReason::MidExchangeMutation;
            v.status = ExchangeStatus::ConfigRejected;
            v.message = "an exchange is in flight; overrides may only be applied at an exchange "
                        "boundary and the running exchange keeps its immutable snapshot";
            v.requirement = "REQ-API-03";
            r.add(v);
            return r;
        }
        if (!d_effective.ok) {
            ConfigViolation v;
            v.field = "";
            v.reason = ConfigReason::OutOfScope;
            v.status = ExchangeStatus::ConfigRejected;
            v.message = "the base config is not valid; overrides cannot be applied to it";
            v.requirement = "REQ-API-03";
            r.add(v);
            return r;
        }
        TwrConfig next = d_effective.effective;
        if (o.measurement_count != 0) {
            if (o.measurement_count < 0)
                next.session.measurement_count = 0; // let the validator reject
            else
                next.session.measurement_count = static_cast<uint32_t>(o.measurement_count);
        }
        if (o.measurement_interval_ns != 0)
            next.session.measurement_interval = Duration::from_nanos(o.measurement_interval_ns);
        if (o.poll_to_response_ns != 0)
            next.timing.poll_to_response.value = Duration::from_nanos(o.poll_to_response_ns);
        if (o.response_to_final_ns != 0)
            next.timing.response_to_final.value = Duration::from_nanos(o.response_to_final_ns);
        if (o.post_tx_rx_enable_ns != 0)
            next.timing.post_tx_rx_enable.value = Duration::from_nanos(o.post_tx_rx_enable_ns);
        if (o.rx_timeout_ns != 0)
            next.timeouts.rx_timeout.value = Duration::from_nanos(o.rx_timeout_ns);
        if (o.exchange_timeout_ns != 0)
            next.timeouts.exchange_timeout.value = Duration::from_nanos(o.exchange_timeout_ns);

        out = TwrConfigSnapshot(effective_config(next, *d_caps), *d_caps);
        if (!out.d_effective.ok) {
            r = out.d_effective.validation;
            return r;
        }
        out.d_effective.requested = d_effective.requested; // keep the original
        out.d_effective.changes = diff_fields(d_effective.requested,
            out.d_effective.effective);
        return r;
    }

private:
    EffectiveConfig d_effective;
    const Capabilities* d_caps;
    ExchangeGate d_gate;
};

// ---------------------------------------------------------------------------
// REQ-CAL-01: a calibration is applicable exactly once
// ---------------------------------------------------------------------------

// Apply the configured calibration to a measurement.  Idempotence is
// enforced, not assumed: a second application is rejected so a range can never
// be corrected twice.  `now_monotonic_ns` comes from the host monotonic clock
// and is only used for the expiry check.
inline ExchangeStatus
apply_calibration_once(TimestampCalibrationConfig& cal, int64_t now_monotonic_ns,
    std::string* why)
{
    if (cal.calibration_id.empty()) {
        if (why)
            *why = "no calibration id: an uncalibrated absolute range may not be claimed";
        return ExchangeStatus::CalibrationMissing;
    }
    if (cal.applied_count != 0) {
        if (why)
            *why = "calibration already applied " + twr_int_to_text(cal.applied_count) +
                   " time(s); a calibration is applicable exactly once";
        return ExchangeStatus::CalibrationMissing;
    }
    if (cal.record.valid_until_monotonic_ns != 0 && now_monotonic_ns != 0 &&
        now_monotonic_ns > cal.record.valid_until_monotonic_ns) {
        if (why)
            *why = "calibration " + cal.calibration_id + " expired at " +
                   twr_int_to_text(cal.record.valid_until_monotonic_ns);
        return ExchangeStatus::CalibrationExpired;
    }
    cal.applied_count = 1;
    if (why)
        *why = "calibration " + cal.calibration_id + " applied once";
    return ExchangeStatus::Ok;
}

// ===========================================================================
// 8. JSON import / export (dependency-free)
// ===========================================================================
//
// The repository has NO JSON library dependency (verified: nothing under
// gr-uwb/include or gr-uwb/lib provides one) and this milestone must not add
// one, so a minimal reader/writer for exactly this schema lives here.
//
// SCHEMA NOTES (REQ-OUT-01)
//   * Every value is emitted, including `null` for an absent optional field,
//     so `from_json(to_json(c)) == c` is exact.
//   * 64-bit integers are emitted UNQUOTED when |v| <= 2^53-1 and as a
//     QUOTED DECIMAL STRING otherwise.  The reader accepts both forms for any
//     integer field, so a timestamp beyond the JSON safe range survives a
//     round trip.  A decimal fraction or exponent on an integer field is
//     REJECTED (`type_mismatch`) rather than silently truncated.
//   * Doubles are emitted with %.17g, which round-trips every finite IEEE 754
//     double.  NaN / +Inf / -Inf have no JSON representation and are
//     REJECTED on export, so a corrupt config cannot be serialised into a
//     file that a later run would parse as a valid one.
//   * Durations are emitted as a plain integer number of NANOSECONDS under a
//     `_ns` key.  Native-tick fields carry a `_native_ticks` key and the
//     calibration additionally carries an explicit `link_delay_unit`, so ns,
//     seconds and native ticks can never be confused.
//   * Number parsing uses strtod and therefore expects the C numeric locale,
//     which is what GNU Radio and the Python host run with.
namespace json {

enum class Type { Null, Bool, Int, Double, String, Array, Object };

class Value
{
public:
    Type type = Type::Null;
    bool boolean = false;
    int64_t integer = 0;
    double number = 0.0;
    std::string text;
    std::vector<Value> items;                            // Array
    std::vector<std::pair<std::string, Value>> members;  // Object, ordered

    static Value make_null() { return Value(); }
    static Value make_bool(bool v)
    {
        Value x;
        x.type = Type::Bool;
        x.boolean = v;
        return x;
    }
    static Value make_int(int64_t v)
    {
        Value x;
        x.type = Type::Int;
        x.integer = v;
        return x;
    }
    static Value make_double(double v)
    {
        Value x;
        x.type = Type::Double;
        x.number = v;
        return x;
    }
    static Value make_string(std::string v)
    {
        Value x;
        x.type = Type::String;
        x.text = std::move(v);
        return x;
    }
    static Value make_array()
    {
        Value x;
        x.type = Type::Array;
        return x;
    }
    static Value make_object()
    {
        Value x;
        x.type = Type::Object;
        return x;
    }

    const Value* find(const std::string& key) const
    {
        for (const auto& m : members)
            if (m.first == key)
                return &m.second;
        return nullptr;
    }
    Value* find(const std::string& key)
    {
        for (auto& m : members)
            if (m.first == key)
                return &m.second;
        return nullptr;
    }
    void set(const std::string& key, Value v)
    {
        if (type == Type::Null)
            type = Type::Object;
        for (auto& m : members) {
            if (m.first == key) {
                m.second = std::move(v);
                return;
            }
        }
        members.emplace_back(key, std::move(v));
    }
    void push(Value v)
    {
        if (type == Type::Null)
            type = Type::Array;
        items.push_back(std::move(v));
    }

    bool is_number() const { return type == Type::Int || type == Type::Double; }
    double as_double() const
    {
        return type == Type::Int ? static_cast<double>(integer) : number;
    }
    const char* type_name() const
    {
        switch (type) {
        case Type::Null:
            return "null";
        case Type::Bool:
            return "bool";
        case Type::Int:
            return "integer";
        case Type::Double:
            return "number";
        case Type::String:
            return "string";
        case Type::Array:
            return "array";
        case Type::Object:
            return "object";
        }
        return "invalid";
    }
};

inline constexpr size_t kJsonMaxDepth = 32;
inline constexpr size_t kJsonMaxInputBytes = 4u * 1024u * 1024u;

class Parser
{
public:
    explicit Parser(const std::string& text) : s(text) {}

    bool parse(Value& out, std::string& err)
    {
        if (s.size() > kJsonMaxInputBytes) {
            err = "json input above " + twr_int_to_text(static_cast<int64_t>(kJsonMaxInputBytes)) +
                  " bytes";
            return false;
        }
        skip_ws();
        if (!parse_value(out, err))
            return false;
        skip_ws();
        if (p != s.size()) {
            err = "trailing characters at offset " + twr_int_to_text(static_cast<int64_t>(p));
            return false;
        }
        return true;
    }

private:
    bool fail(const std::string& m, std::string& err)
    {
        err = m + " at offset " + twr_int_to_text(static_cast<int64_t>(p));
        return false;
    }
    void skip_ws()
    {
        while (p < s.size() &&
               (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' || s[p] == '\r'))
            ++p;
    }
    bool parse_value(Value& out, std::string& err)
    {
        if (d_depth >= kJsonMaxDepth)
            return fail("json nesting too deep", err);
        if (p >= s.size())
            return fail("unexpected end of input", err);
        switch (s[p]) {
        case '{':
            return parse_object(out, err);
        case '[':
            return parse_array(out, err);
        case '"': {
            std::string t;
            if (!parse_string(t, err))
                return false;
            out = Value::make_string(std::move(t));
            return true;
        }
        case 't':
            if (s.compare(p, 4, "true") != 0)
                return fail("bad literal", err);
            p += 4;
            out = Value::make_bool(true);
            return true;
        case 'f':
            if (s.compare(p, 5, "false") != 0)
                return fail("bad literal", err);
            p += 5;
            out = Value::make_bool(false);
            return true;
        case 'n':
            if (s.compare(p, 4, "null") != 0)
                return fail("bad literal", err);
            p += 4;
            out = Value::make_null();
            return true;
        default:
            return parse_number(out, err);
        }
    }
    bool parse_object(Value& out, std::string& err)
    {
        ++p; // '{'
        ++d_depth;
        out = Value::make_object();
        skip_ws();
        if (p < s.size() && s[p] == '}') {
            ++p;
            --d_depth;
            return true;
        }
        while (true) {
            skip_ws();
            if (p >= s.size() || s[p] != '"')
                return fail("expected object key", err);
            std::string key;
            if (!parse_string(key, err))
                return false;
            skip_ws();
            if (p >= s.size() || s[p] != ':')
                return fail("expected ':'", err);
            ++p;
            skip_ws();
            Value v;
            if (!parse_value(v, err))
                return false;
            out.set(key, std::move(v));
            skip_ws();
            if (p < s.size() && s[p] == ',') {
                ++p;
                continue;
            }
            if (p < s.size() && s[p] == '}') {
                ++p;
                --d_depth;
                return true;
            }
            return fail("expected ',' or '}'", err);
        }
    }
    bool parse_array(Value& out, std::string& err)
    {
        ++p; // '['
        ++d_depth;
        out = Value::make_array();
        skip_ws();
        if (p < s.size() && s[p] == ']') {
            ++p;
            --d_depth;
            return true;
        }
        while (true) {
            skip_ws();
            Value v;
            if (!parse_value(v, err))
                return false;
            out.push(std::move(v));
            skip_ws();
            if (p < s.size() && s[p] == ',') {
                ++p;
                continue;
            }
            if (p < s.size() && s[p] == ']') {
                ++p;
                --d_depth;
                return true;
            }
            return fail("expected ',' or ']'", err);
        }
    }
    bool parse_string(std::string& out, std::string& err)
    {
        ++p; // opening quote
        out.clear();
        while (p < s.size()) {
            const char c = s[p++];
            if (c == '"')
                return true;
            if (c == '\\') {
                if (p >= s.size())
                    return fail("unterminated escape", err);
                const char e = s[p++];
                switch (e) {
                case '"':
                    out.push_back('"');
                    break;
                case '\\':
                    out.push_back('\\');
                    break;
                case '/':
                    out.push_back('/');
                    break;
                case 'b':
                    out.push_back('\b');
                    break;
                case 'f':
                    out.push_back('\f');
                    break;
                case 'n':
                    out.push_back('\n');
                    break;
                case 'r':
                    out.push_back('\r');
                    break;
                case 't':
                    out.push_back('\t');
                    break;
                case 'u': {
                    if (p + 4 > s.size())
                        return fail("short \\u escape", err);
                    unsigned cp = 0;
                    for (int i = 0; i < 4; ++i) {
                        const char h = s[p + static_cast<size_t>(i)];
                        cp <<= 4;
                        if (h >= '0' && h <= '9')
                            cp |= static_cast<unsigned>(h - '0');
                        else if (h >= 'a' && h <= 'f')
                            cp |= static_cast<unsigned>(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F')
                            cp |= static_cast<unsigned>(h - 'A' + 10);
                        else
                            return fail("bad \\u escape", err);
                    }
                    p += 4;
                    // UTF-8 encode (surrogate pairs are passed through as-is;
                    // this config only ever carries ASCII identifiers).
                    if (cp < 0x80) {
                        out.push_back(static_cast<char>(cp));
                    } else if (cp < 0x800) {
                        out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
                        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
                    } else {
                        out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
                        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
                    }
                    break;
                }
                default:
                    return fail("unknown escape", err);
                }
                continue;
            }
            out.push_back(c);
        }
        return fail("unterminated string", err);
    }
    bool parse_number(Value& out, std::string& err)
    {
        const size_t start = p;
        if (p < s.size() && s[p] == '-')
            ++p;
        if (p >= s.size() || s[p] < '0' || s[p] > '9')
            return fail("expected a number", err);
        if (s[p] == '0') {
            ++p;
        } else {
            while (p < s.size() && s[p] >= '0' && s[p] <= '9')
                ++p;
        }
        bool integral = true;
        if (p < s.size() && s[p] == '.') {
            integral = false;
            ++p;
            if (p >= s.size() || s[p] < '0' || s[p] > '9')
                return fail("expected digits after '.'", err);
            while (p < s.size() && s[p] >= '0' && s[p] <= '9')
                ++p;
        }
        if (p < s.size() && (s[p] == 'e' || s[p] == 'E')) {
            integral = false;
            ++p;
            if (p < s.size() && (s[p] == '+' || s[p] == '-'))
                ++p;
            if (p >= s.size() || s[p] < '0' || s[p] > '9')
                return fail("expected digits in exponent", err);
            while (p < s.size() && s[p] >= '0' && s[p] <= '9')
                ++p;
        }
        const std::string tok = s.substr(start, p - start);
        if (integral) {
            errno = 0;
            char* end = nullptr;
            const long long v = std::strtoll(tok.c_str(), &end, 10);
            if (errno == 0 && end != nullptr && *end == '\0') {
                out = Value::make_int(static_cast<int64_t>(v));
                return true;
            }
            // An integer literal outside int64 is kept as a Double so the
            // reader for an integer field can reject it explicitly.
        }
        errno = 0;
        char* end = nullptr;
        const double d = std::strtod(tok.c_str(), &end);
        if (end == nullptr || *end != '\0')
            return fail("malformed number", err);
        out = Value::make_double(d);
        return true;
    }

    const std::string& s;
    size_t p = 0;
    size_t d_depth = 0;
};

inline bool parse(const std::string& text, Value& out, std::string& err)
{
    Parser parser(text);
    return parser.parse(out, err);
}

inline void escape_into(const std::string& in, std::string& out)
{
    out.push_back('"');
    for (unsigned char c : in) {
        switch (c) {
        case '"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (c < 0x20) {
                char buf[8];
                std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned>(c));
                out += buf;
            } else {
                out.push_back(static_cast<char>(c));
            }
        }
    }
    out.push_back('"');
}

// 64-bit integer as a JSON number when safe, else as a decimal string.
inline void put_int(int64_t v, std::string& out)
{
    if (v <= kTwrJsonMaxSafeInteger && v >= -kTwrJsonMaxSafeInteger) {
        out += twr_int_to_text(v);
        return;
    }
    escape_into(twr_int_to_text(v), out);
}

// -- dump ------------------------------------------------------------------

class Dumper
{
public:
    Dumper(const Value& v, bool pretty) : d_v(v), d_pretty(pretty) {}

    bool run(std::string& out, std::string& err)
    {
        m_buf.clear();
        if (!dump(d_v, 0, err))
            return false;
        out = m_buf;
        if (d_pretty)
            out += "\n";
        return true;
    }

private:
    void indent(int depth)
    {
        if (!d_pretty)
            return;
        m_buf += "\n";
        m_buf.append(static_cast<size_t>(depth) * 2, ' ');
    }
    bool dump(const Value& v, int depth, std::string& err)
    {
        switch (v.type) {
        case Type::Null:
            m_buf += "null";
            return true;
        case Type::Bool:
            m_buf += v.boolean ? "true" : "false";
            return true;
        case Type::Int:
            put_int(v.integer, m_buf);
            return true;
        case Type::Double:
            if (!std::isfinite(v.number)) {
                err = "non-finite double cannot be represented in JSON";
                return false;
            }
            m_buf += twr_double_to_text(v.number);
            return true;
        case Type::String:
            escape_into(v.text, m_buf);
            return true;
        case Type::Array: {
            if (v.items.empty()) {
                m_buf += "[]";
                return true;
            }
            m_buf += "[";
            for (size_t i = 0; i < v.items.size(); ++i) {
                if (i)
                    m_buf += ",";
                indent(depth + 1);
                if (!dump(v.items[i], depth + 1, err))
                    return false;
            }
            indent(depth);
            m_buf += "]";
            return true;
        }
        case Type::Object: {
            if (v.members.empty()) {
                m_buf += "{}";
                return true;
            }
            m_buf += "{";
            for (size_t i = 0; i < v.members.size(); ++i) {
                if (i)
                    m_buf += ",";
                indent(depth + 1);
                escape_into(v.members[i].first, m_buf);
                m_buf += ":";
                if (d_pretty)
                    m_buf += " ";
                if (!dump(v.members[i].second, depth + 1, err))
                    return false;
            }
            indent(depth);
            m_buf += "}";
            return true;
        }
        }
        err = "unreachable json type";
        return false;
    }

    const Value& d_v;
    bool d_pretty;
    std::string m_buf;
};

inline bool dump(const Value& v, std::string& out, std::string& err, bool pretty = true)
{
    Dumper d(v, pretty);
    return d.run(out, err);
}

} // namespace json

// ===========================================================================
// 9. TwrConfig <-> JSON
// ===========================================================================

// Small ordered object writer.  Emits every field, so the round trip is
// exact and a config file documents itself.
class ConfigJsonWriter
{
public:
    ConfigJsonWriter() { m_stack.push_back(&m_root); }

    bool ok() const { return m_err.empty(); }
    const std::string& error() const { return m_err; }
    json::Value value() const { return m_root; }
    std::string to_string(bool pretty = true) const
    {
        std::string out;
        std::string err;
        if (!json::dump(m_root, out, err, pretty)) {
            return std::string();
        }
        return out;
    }

    // The root is already an object; this exists so the document shape reads
    // the same as the reader.
    void current_object() {}

    void begin_group(const std::string& key)
    {
        m_stack.back()->set(key, json::Value::make_object());
        m_stack.push_back(m_stack.back()->find(key));
    }
    void begin_array(const std::string& key)
    {
        m_stack.back()->set(key, json::Value::make_array());
        m_stack.push_back(m_stack.back()->find(key));
    }
    void end() { m_stack.pop_back(); }

    void push_object()
    {
        json::Value* parent = m_stack.back();
        parent->push(json::Value::make_object());
        m_stack.push_back(&parent->items.back());
    }
    void pop_object() { m_stack.pop_back(); }

    void put(const std::string& key, const std::string& v)
    {
        current()->set(key, json::Value::make_string(v));
    }
    void put(const std::string& key, int64_t v)
    {
        current()->set(key, json::Value::make_int(v));
    }
    void put(const std::string& key, uint32_t v)
    {
        current()->set(key, json::Value::make_int(static_cast<int64_t>(v)));
    }
    void put(const std::string& key, uint16_t v)
    {
        current()->set(key, json::Value::make_int(static_cast<int64_t>(v)));
    }
    void put(const std::string& key, uint8_t v)
    {
        current()->set(key, json::Value::make_int(static_cast<int64_t>(v)));
    }
    void put(const std::string& key, int v)
    {
        current()->set(key, json::Value::make_int(static_cast<int64_t>(v)));
    }
    void put(const std::string& key, uint64_t v)
    {
        // uint64 above int64 max is not representable; the schema has no such
        // field, but keep the check explicit rather than wrapping silently.
        if (v > static_cast<uint64_t>(INT64_MAX)) {
            fail("field " + key + " exceeds the int64 range of the JSON schema");
            current()->set(key, json::Value::make_null());
            return;
        }
        current()->set(key, json::Value::make_int(static_cast<int64_t>(v)));
    }
    void put(const std::string& key, bool v)
    {
        current()->set(key, json::Value::make_bool(v));
    }
    void put(const std::string& key, double v)
    {
        if (!std::isfinite(v)) {
            fail("field " + key + " is not finite and cannot be written as JSON");
            current()->set(key, json::Value::make_null());
            return;
        }
        current()->set(key, json::Value::make_double(v));
    }
    void put(const std::string& key, const char* v)
    {
        current()->set(key, json::Value::make_string(v ? v : ""));
    }
    void put_ns(const std::string& key, const Duration& d)
    {
        current()->set(key, json::Value::make_int(d.nanos()));
    }
    void put_opt(const std::string& key, const Opt<int64_t>& o)
    {
        current()->set(key,
                       o.has() ? json::Value::make_int(o.value()) : json::Value::make_null());
    }
    void put_opt(const std::string& key, const Opt<uint32_t>& o)
    {
        current()->set(key,
                       o.has() ? json::Value::make_int(static_cast<int64_t>(o.value()))
                               : json::Value::make_null());
    }
    void put_opt(const std::string& key, const Opt<double>& o)
    {
        if (!o.has()) {
            current()->set(key, json::Value::make_null());
            return;
        }
        if (!std::isfinite(o.value())) {
            fail("field " + key + " is not finite and cannot be written as JSON");
            current()->set(key, json::Value::make_null());
            return;
        }
        current()->set(key, json::Value::make_double(o.value()));
    }
    void put_timed(const std::string& key, const TimedField& t)
    {
        json::Value o = json::Value::make_object();
        o.set("ns", json::Value::make_int(t.value.nanos()));
        o.set("domain", json::Value::make_string(time_domain_to_string(t.domain)));
        o.set("reference",
            json::Value::make_string(time_reference_event_to_string(t.reference)));
        if (t.marker.has())
            o.set("marker",
                json::Value::make_string(timestamp_marker_to_string(t.marker.value())));
        else
            o.set("marker", json::Value::make_null());
        o.set("quantisation_hz", json::Value::make_double(t.required_quantisation_hz));
        o.set("max_quantisation_error_ns",
              json::Value::make_int(t.max_quantisation_error_ns));
        o.set("note", json::Value::make_string(t.note));
        current()->set(key, std::move(o));
    }

private:
    json::Value* current() { return m_stack.back(); }
    void fail(const std::string& m)
    {
        if (m_err.empty())
            m_err = m;
    }

    json::Value m_root;
    std::vector<json::Value*> m_stack;
    std::string m_err;
};

void write_config_json(const TwrConfig& c, ConfigJsonWriter& w);

// Strict reader: every failure is a machine-readable violation, never a
// silently skipped key and never a default.
class ConfigJsonReader
{
public:
    void bad(const std::string& field,
             ConfigReason r,
             const std::string& msg,
             ExchangeStatus st = ExchangeStatus::ConfigRejected)
    {
        m_report.add(field, r, st, msg, "REQ-API-01");
    }

    ValidationReport report() const { return m_report; }
    bool ok() const { return m_report.ok(); }

    // -- object navigation ------------------------------------------------
    const json::Value* group(const std::string& path, const json::Value* parent)
    {
        const json::Value* g = parent ? parent->find(path) : nullptr;
        if (!g) {
            bad(path, ConfigReason::MissingKey, "missing required group");
            return nullptr;
        }
        if (g->type != json::Type::Object) {
            bad(path, ConfigReason::TypeMismatch, "expected an object, got " +
                                                         std::string(g->type_name()));
            return nullptr;
        }
        return g;
    }

    // -- typed getters ----------------------------------------------------
    std::string get_str(const json::Value* o, const std::string& path, const char* key)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return std::string();
        if (v->type == json::Type::Null)
            return std::string();
        if (v->type != json::Type::String) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected a string, got " + std::string(v->type_name()));
            return std::string();
        }
        return v->text;
    }
    bool get_bool(const json::Value* o, const std::string& path, const char* key)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return false;
        if (v->type != json::Type::Bool) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected a bool, got " + std::string(v->type_name()));
            return false;
        }
        return v->boolean;
    }
    int64_t get_i64(const json::Value* o, const std::string& path, const char* key)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return 0;
        if (v->type == json::Type::String) {
            // The decimal-string form used for values beyond 2^53.
            int64_t out = 0;
            if (parse_int64_strict(v->text, out))
                return out;
            bad(path + "." + key, ConfigReason::IntegerPrecisionLoss,
                "string is not a decimal int64: '" + v->text + "'");
            return 0;
        }
        if (v->type != json::Type::Int) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected an integer, got " + std::string(v->type_name()) +
                    " (a decimal fraction or exponent is never truncated silently)");
            return 0;
        }
        return v->integer;
    }

    // -- range-checked unsigned narrowing (M0.1 / N01) --------------------
    //
    // The reader used to assign through `static_cast<uint16_t>(get_i64(...))`,
    // which WRAPS an out-of-range JSON number instead of refusing it: channel
    // 261 became channel 5, session 2^32+1 became session 1, local_address
    // 65537 became 1.  The validator then never saw the value the document
    // actually carried, and two different documents could produce the same
    // runtime configuration -- exactly the silent fallback REQ-CFG forbids and
    // the reason `session_id_to_wire()` exists.
    //
    // A C++ cast is not a specification.  These getters are: the value must be
    // in [0, max] or it is refused with `out_of_range` and stored as 0.  They
    // keep `get_i64`'s shape (same arguments) so the call sites read the same
    // way, but the range check happens BEFORE any narrowing.
    bool check_unsigned(const std::string& path, const char* key, int64_t v,
                        uint64_t max_value, uint64_t& out)
    {
        // Formatted unsigned: `twr_int_to_text` takes an int64_t, so passing a
        // uint64_t max would print a negative bound.
        const std::string bound = std::to_string(
            static_cast<unsigned long long>(max_value));
        if (v < 0) {
            bad(path + "." + key, ConfigReason::OutOfRange,
                "must be in [0, " + bound + "], got " + twr_int_to_text(v) +
                    " (a negative value is refused, never wrapped to fit the field)");
            out = 0u;
            return false;
        }
        const uint64_t u = static_cast<uint64_t>(v);
        if (u > max_value) {
            bad(path + "." + key, ConfigReason::OutOfRange,
                "must be in [0, " + bound + "], got " + twr_int_to_text(v) +
                    " (out of range for this field; refused, not masked)");
            out = 0u;
            return false;
        }
        out = u;
        return true;
    }

    uint64_t get_u64(const json::Value* o, const std::string& path, const char* key,
                     uint64_t max_value)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return 0u;
        int64_t raw = 0;
        if (v->type == json::Type::String) {
            if (!parse_int64_strict(v->text, raw)) {
                bad(path + "." + key, ConfigReason::IntegerPrecisionLoss,
                    "string is not a decimal int64: '" + v->text + "'");
                return 0u;
            }
        } else if (v->type != json::Type::Int) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected an integer, got " + std::string(v->type_name()) +
                    " (a decimal fraction or exponent is never truncated silently)");
            return 0u;
        } else {
            raw = v->integer;
        }
        uint64_t out = 0u;
        check_unsigned(path, key, raw, max_value, out);
        return out;
    }

    uint8_t get_u8(const json::Value* o, const std::string& path, const char* key)
    {
        return static_cast<uint8_t>(get_u64(o, path, key, 0xFFu));
    }
    uint16_t get_u16(const json::Value* o, const std::string& path, const char* key)
    {
        return static_cast<uint16_t>(get_u64(o, path, key, 0xFFFFu));
    }
    uint32_t get_u32(const json::Value* o, const std::string& path, const char* key)
    {
        return static_cast<uint32_t>(get_u64(o, path, key, 0xFFFFFFFFu));
    }
    double get_double(const json::Value* o, const std::string& path, const char* key)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return 0.0;
        if (!v->is_number()) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected a number, got " + std::string(v->type_name()));
            return 0.0;
        }
        return v->as_double();
    }
    Duration get_dur(const json::Value* o, const std::string& path, const char* key)
    {
        return Duration::from_nanos(get_i64(o, path, key));
    }
    // -- optional getters ------------------------------------------------
    Opt<int64_t> get_opt_i64(const json::Value* o, const std::string& path,
        const char* key)
    {
        const json::Value* v = key_of(o, path, key, /*required=*/false);
        if (!v || v->type == json::Type::Null)
            return Opt<int64_t>::none();
        const std::string p = path + "." + key;
        Opt<int64_t> out;
        if (v->type == json::Type::String) {
            int64_t parsed = 0;
            if (!parse_int64_strict(v->text, parsed)) {
                bad(p, ConfigReason::IntegerPrecisionLoss,
                    "string is not a decimal int64: '" + v->text + "'");
                return Opt<int64_t>::none();
            }
            out = parsed;
        } else if (v->type == json::Type::Int) {
            out = v->integer;
        } else {
            bad(p, ConfigReason::TypeMismatch,
                "expected an integer or null, got " + std::string(v->type_name()));
        }
        return out;
    }
    Opt<uint32_t> get_opt_u32(const json::Value* o, const std::string& path,
        const char* key)
    {
        const Opt<int64_t> v = get_opt_i64(o, path, key);
        if (!v.has())
            return Opt<uint32_t>::none();
        if (v.value() < 0 || v.value() > 0xffffffffLL) {
            bad(path + "." + key, ConfigReason::OutOfRange, "outside the uint32 range");
            return Opt<uint32_t>::none();
        }
        return Opt<uint32_t>(static_cast<uint32_t>(v.value()));
    }
    Opt<double> get_opt_double(const json::Value* o, const std::string& path,
        const char* key)
    {
        const json::Value* v = key_of(o, path, key, /*required=*/false);
        if (!v || v->type == json::Type::Null)
            return Opt<double>::none();
        if (!v->is_number()) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected a number or null, got " + std::string(v->type_name()));
            return Opt<double>::none();
        }
        return Opt<double>(v->as_double());
    }
    // -- enum getters ----------------------------------------------------
    template <typename E, typename F>
    E get_enum(const json::Value* o, const std::string& path, const char* key,
        F from_string)
    {
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return static_cast<E>(0);
        if (v->type != json::Type::String) {
            bad(path + "." + key, ConfigReason::TypeMismatch,
                "expected a string enum, got " + std::string(v->type_name()));
            return static_cast<E>(0);
        }
        E out = static_cast<E>(0);
        if (!from_string(v->text, out)) {
            bad(path + "." + key, ConfigReason::UnknownEnumValue,
                "'" + v->text + "' is not a known value of this enum");
            return static_cast<E>(0);
        }
        return out;
    }
    TimedField get_timed(const json::Value* o, const std::string& path, const char* key)
    {
        TimedField t;
        const std::string p = path + "." + key;
        const json::Value* v = key_of(o, path, key);
        if (!v)
            return t;
        if (v->type != json::Type::Object) {
            bad(p, ConfigReason::TypeMismatch,
                "expected an object, got " + std::string(v->type_name()));
            return t;
        }
        t.value = Duration::from_nanos(get_i64(v, p, "ns"));
        t.domain = get_enum<TimeDomain>(v, p, "domain", time_domain_from_string);
        t.reference = get_enum<TimeReferenceEvent>(v, p, "reference",
                                                   time_reference_event_from_string);
        const json::Value* mk = key_of(v, p, "marker", /*required=*/false);
        if (mk && mk->type == json::Type::String) {
            bool found = false;
            for (int i = 0; i <= static_cast<int>(TimestampMarker::AntennaPlane); ++i) {
                const TimestampMarker m = static_cast<TimestampMarker>(i);
                if (mk->text == timestamp_marker_to_string(m)) {
                    t.marker = m;
                    found = true;
                    break;
                }
            }
            if (!found)
                bad(p + ".marker", ConfigReason::UnknownEnumValue,
                    "'" + mk->text + "' is not a known timestamp marker");
        } else if (mk && mk->type != json::Type::Null) {
            bad(p + ".marker", ConfigReason::TypeMismatch,
                "expected a marker name or null, got " + std::string(mk->type_name()));
        }
        t.required_quantisation_hz = get_double(v, p, "quantisation_hz");
        t.max_quantisation_error_ns = get_i64(v, p, "max_quantisation_error_ns");
        t.note = get_str(v, p, "note");
        return t;
    }

    static bool parse_int64_strict(const std::string& s, int64_t& out)
    {
        if (s.empty())
            return false;
        size_t i = 0;
        if (s[0] == '-')
            i = 1;
        if (i >= s.size())
            return false;
        for (size_t k = i; k < s.size(); ++k)
            if (s[k] < '0' || s[k] > '9')
                return false;
        errno = 0;
        char* end = nullptr;
        const long long v = std::strtoll(s.c_str(), &end, 10);
        if (errno != 0 || end == nullptr || *end != '\0')
            return false;
        out = static_cast<int64_t>(v);
        return true;
    }

    // Reject keys the schema does not define: an unknown key is a typo or an
    // unsupported feature, and either way it must not be ignored.
    void reject_unknown_keys(const json::Value* o, const std::string& path,
        const char* const* known)
    {
        if (!o || o->type != json::Type::Object)
            return;
        for (const auto& m : o->members) {
            bool found = false;
            for (size_t i = 0; known[i] != nullptr; ++i) {
                if (m.first == known[i]) {
                    found = true;
                    break;
                }
            }
            if (!found)
                bad(path + "." + m.first, ConfigReason::UnknownKey,
                    "key is not part of schema twr-config/2");
        }
    }

private:
    const json::Value* key_of(const json::Value* o,
                              const std::string& path,
                              const char* key,
                              bool required = true)
    {
        if (!o)
            return nullptr;
        const json::Value* v = o->find(key);
        if (!v) {
            if (required)
                bad(path + "." + key, ConfigReason::MissingKey, "missing required key");
            return nullptr;
        }
        return v;
    }
    ValidationReport m_report;
};

// ---------------------------------------------------------------------------

bool to_json_string(const TwrConfig& c, std::string& out, std::string& err);

// Convenience: the frozen snapshot plus its report, as one JSON document.
bool to_json_string(const EffectiveConfig& e, std::string& out, std::string& err);

ValidationReport from_json_string(const std::string& text, TwrConfig& out);

} // namespace twr
} // namespace openuwb

#endif /* INCLUDED_OPENUWB_TWR_CONFIG_EFFECTIVE_H */
