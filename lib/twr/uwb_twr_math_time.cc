// TWR math/time implementation moved out of the public headers.
//
// Plan: docs/performance/多Subagent并行开发方案_按功能裁剪与编译测试再加速_2026-11-21.md
// (B5) DEC-B half: twr_math / twr_protocol_time / twr_tof_input /
// twr_capability_evidence.  Short accessors stay inline in the headers.
#include <openuwb/twr/twr_capability_evidence.h>
#include <openuwb/twr/twr_math.h>
#include <openuwb/twr/twr_protocol_time.h>
#include <openuwb/twr/twr_tof_input.h>

namespace openuwb {
namespace twr {


TofResult compute_ss_tof(const AdmittedRangingInterval& ra,
                         const AdmittedRangingInterval& db,
                         const ClockRatio& k_ab)
{
    TofResult r;
    r.protocol = Protocol::Ss;
    r.computed_at = ComputedAt::InitiatorA; // REQ-PROTO-05: baseline at A
    r.ratio = k_ab;
    r.domain = k_ab.domain_a();

    const TofStatus ks = k_ab.check();
    if (ks != TofStatus::Ok)
        return detail::tof_fail(r, ks, "clock ratio is unusable: " + k_ab.to_string());

    if (!detail::tof_domain_check(r, "ra", ra.domain(), k_ab.domain_a()))
        return r;
    if (!detail::tof_domain_check(r, "db", db.domain(), k_ab.domain_b()))
        return r;
    // REQ-TIME-04: the ratio must be in force across the A-domain interval.
    if (!detail::tof_window_check(r, "ra", ra, k_ab))
        return r;

    detail::Rat128 RA;
    detail::Rat128 DB;
    if (!detail::tof_load_interval(r, "ra", ra, RA))
        return r;
    if (!detail::tof_load_interval(r, "db", db, DB))
        return r;

    const detail::Rat128 k = detail::rat128_make(k_ab.k_num(), k_ab.k_den());
    detail::Rat128 kdb;
    detail::Rat128 half;
    if (!detail::rat128_mul(k, DB, kdb))
        return detail::tof_fail(r, TofStatus::Overflow, "k*DB left 128 bits");
    // The formula itself lives in the ONE shared kernel (M1-B §3.3); this
    // strict entry point keeps its own admission checks above and maps the
    // kernel status back onto its own vocabulary below.
    if (detail::ss_tof_kernel(RA, DB, k, half) != detail::KernelStatus::Ok)
        return detail::tof_fail(r, TofStatus::Overflow, "SS exact kernel left 128 bits");
    if (!detail::rat128_fits_int64(half))
        return detail::tof_fail(
            r, TofStatus::Overflow, "ToF numerator/denominator do not fit int64 ticks");

    r.ra_raw = detail::rat128_to_ticks(RA);
    r.db_raw = detail::rat128_to_ticks(DB);
    r.ra = r.ra_raw;                     // already A ticks
    r.db = detail::rat128_to_ticks(kdb); // k*DB, converted to A ticks
    r.tof = detail::rat128_to_ticks(half);
    r.display_seconds = r.tof.to_seconds(r.domain.tick_rate_hz);

    if (r.tof.num < 0) {
        r.ok = false;
        r.status = TofStatus::NegativeTof;
        r.exchange_status = tof_status_to_exchange_status(TofStatus::NegativeTof);
        r.detail = "SS ToF is negative (" + r.tof.to_string() +
                   " A ticks): retained with its sign, never clamped to 0";
        return r;
    }
    r.ok = true;
    r.status = TofStatus::Ok;
    r.exchange_status = ExchangeStatus::Ok;
    r.detail = "SS ToF = (RA - k*DB)/2";
    return r;
}


TofResult compute_ds_tof(const AdmittedRangingInterval& ra,
                         const AdmittedRangingInterval& rb,
                         const AdmittedRangingInterval& da,
                         const AdmittedRangingInterval& db,
                         const ClockRatio& k_ab)
{
    TofResult r;
    r.protocol = Protocol::Ds;
    r.computed_at = ComputedAt::ResponderB; // REQ-PROTO-05: baseline at B
    r.ratio = k_ab;
    r.domain = k_ab.domain_a(); // the common unit is A ticks

    const TofStatus ks = k_ab.check();
    if (ks != TofStatus::Ok)
        return detail::tof_fail(r, ks, "clock ratio is unusable: " + k_ab.to_string());

    if (!detail::tof_domain_check(r, "ra", ra.domain(), k_ab.domain_a()))
        return r;
    if (!detail::tof_domain_check(r, "da", da.domain(), k_ab.domain_a()))
        return r;
    if (!detail::tof_domain_check(r, "rb", rb.domain(), k_ab.domain_b()))
        return r;
    if (!detail::tof_domain_check(r, "db", db.domain(), k_ab.domain_b()))
        return r;
    if (!detail::tof_window_check(r, "ra", ra, k_ab))
        return r;
    if (!detail::tof_window_check(r, "da", da, k_ab))
        return r;

    detail::Rat128 RA;
    detail::Rat128 RB;
    detail::Rat128 DA;
    detail::Rat128 DB;
    if (!detail::tof_load_interval(r, "ra", ra, RA))
        return r;
    if (!detail::tof_load_interval(r, "rb", rb, RB))
        return r;
    if (!detail::tof_load_interval(r, "da", da, DA))
        return r;
    if (!detail::tof_load_interval(r, "db", db, DB))
        return r;

    const detail::Rat128 k = detail::rat128_make(k_ab.k_num(), k_ab.k_den());
    detail::Rat128 rba;
    detail::Rat128 dba;
    if (!detail::rat128_mul(k, RB, rba))
        return detail::tof_fail(r, TofStatus::Overflow, "k*RB left 128 bits");
    if (!detail::rat128_mul(k, DB, dba))
        return detail::tof_fail(r, TofStatus::Overflow, "k*DB left 128 bits");

    // The formula itself lives in the ONE shared kernel (M1-B §3.3); the
    // strict entry point keeps its own admission checks above.
    detail::Rat128 tof;
    const detail::KernelStatus kst = detail::ds_tof_kernel(RA, RB, DA, DB, k, tof);
    if (kst == detail::KernelStatus::ZeroDenominator) {
        // Four non-negative intervals sum to something non-positive: only the
        // all-zero case is reachable, and the quotient is undefined.
        r.ok = false;
        r.status = TofStatus::ZeroDenominator;
        r.exchange_status = tof_status_to_exchange_status(TofStatus::ZeroDenominator);
        r.detail = "DS denominator is not positive: RA+k*RB+DA+k*DB is non-positive";
        return r;
    }
    if (kst != detail::KernelStatus::Ok)
        return detail::tof_fail(r, TofStatus::Overflow, "DS exact kernel left 128 bits");
    if (!detail::rat128_fits_int64(tof))
        return detail::tof_fail(
            r, TofStatus::Overflow, "ToF numerator/denominator do not fit int64 ticks");

    r.ra_raw = detail::rat128_to_ticks(RA);
    r.rb_raw = detail::rat128_to_ticks(RB);
    r.da_raw = detail::rat128_to_ticks(DA);
    r.db_raw = detail::rat128_to_ticks(DB);
    r.ra = r.ra_raw;                     // A ticks already
    r.da = r.da_raw;                     // A ticks already
    r.rb = detail::rat128_to_ticks(rba); // converted to A ticks
    r.db = detail::rat128_to_ticks(dba); // converted to A ticks
    r.tof = detail::rat128_to_ticks(tof);
    r.display_seconds = r.tof.to_seconds(r.domain.tick_rate_hz);

    if (r.tof.num < 0) {
        r.ok = false;
        r.status = TofStatus::NegativeTof;
        r.exchange_status = tof_status_to_exchange_status(TofStatus::NegativeTof);
        r.detail = "DS ToF is negative (" + r.tof.to_string() +
                   " A ticks): retained with its sign, never clamped to 0";
        return r;
    }
    r.ok = true;
    r.status = TofStatus::Ok;
    r.exchange_status = ExchangeStatus::Ok;
    r.detail = "DS ToF = (RA*k*RB - DA*k*DB)/(RA+k*RB+DA+k*DB), all in A ticks";
    return r;
}


namespace detail {

TofResult tof_fail(TofResult r, TofStatus s, const std::string& detail)
{
    r.ok = false;
    r.status = s;
    r.exchange_status = tof_status_to_exchange_status(s);
    r.detail = detail;
    return r;
}


bool tof_domain_check(TofResult& r,
                      const char* what,
                      const ClockDomain& got,
                      const ClockDomain& want)
{
    if (!got.is_valid()) {
        r = tof_fail(r,
                     TofStatus::InvalidClockDomain,
                     std::string(what) + ": interval carries an invalid clock domain");
        return false;
    }
    if (!clock_domain_same_identity(got, want)) {
        r = tof_fail(r,
                     TofStatus::ClockRatioDomainMismatch,
                     std::string(what) + " is in " + got.to_string() +
                         " but the clock ratio declares " + want.to_string());
        return false;
    }
    if (!clock_domain_same_epoch(got, want)) {
        r = tof_fail(r,
                     TofStatus::ClockRatioEpochMismatch,
                     std::string(what) + " is epoch " + std::to_string(got.epoch_id) +
                         " but the clock ratio declares epoch " +
                         std::to_string(want.epoch_id));
        return false;
    }
    return true;
}


bool tof_window_check(TofResult& r,
                      const char* what,
                      const AdmittedRangingInterval& iv,
                      const ClockRatio& k)
{
    if (!k.has_validity_window())
        return true;
    const int64_t a = iv.earlier().ticks();
    const int64_t b = iv.later().ticks();
    const int64_t lo = (a < b) ? a : b;
    const int64_t hi = (a < b) ? b : a;
    if (!k.covers_ticks(lo) || !k.covers_ticks(hi)) {
        r = tof_fail(r,
                     TofStatus::ClockRatioNotValidAtTime,
                     std::string(what) + " spans A ticks [" + std::to_string(lo) + "," +
                         std::to_string(hi) + "] which is outside the clock ratio's " +
                         "validity window [" + std::to_string(k.valid_from_ticks()) +
                         "," + std::to_string(k.valid_until_ticks()) + "]");
        return false;
    }
    return true;
}


bool tof_load_interval(TofResult& r,
                       const char* what,
                       const AdmittedRangingInterval& iv,
                       Rat128& out)
{
    if (!iv.interval().is_valid()) {
        r = tof_fail(r,
                     TofStatus::InvalidInput,
                     std::string(what) + ": admitted interval is not structurally valid");
        return false;
    }
    if (!rat128_from_interval(iv.interval(), out)) {
        r = tof_fail(r,
                     TofStatus::Overflow,
                     std::string(what) + ": exact_ratio() did not fit int64");
        return false;
    }
    return true;
}


} // namespace detail


ProtocolTofEstimate compute_protocol_ss_tof(const ProtocolInterval& ra,
                                            const ProtocolInterval& db,
                                            const ClockRatio& k_ab)
{
    ProtocolTofEstimate e;
    e.protocol = Protocol::Ss;
    e.computed_at = ComputedAt::InitiatorA;
    e.ratio = k_ab;
    e.domain = k_ab.domain_a();

    const TofStatus ks = k_ab.check();
    if (ks != TofStatus::Ok)
        return detail::protocol_fail(
            e, ks, "clock ratio is unusable: " + k_ab.to_string());

    if (!detail::protocol_domain_ok(e, "ra", ra, k_ab.domain_a()))
        return e;
    if (!detail::protocol_domain_ok(e, "db", db, k_ab.domain_b()))
        return e;
    if (!detail::protocol_window_ok(e, "ra", ra, k_ab))
        return e;

    detail::Rat128 RA;
    detail::Rat128 DB;
    if (!detail::protocol_load_interval(e, "ra", ra, RA))
        return e;
    if (!detail::protocol_load_interval(e, "db", db, DB))
        return e;

    const detail::Rat128 k = detail::rat128_make(k_ab.k_num(), k_ab.k_den());
    detail::Rat128 half;
    const detail::KernelStatus st = detail::ss_tof_kernel(RA, DB, k, half);
    if (st != detail::KernelStatus::Ok)
        return detail::protocol_fail(
            e, detail::protocol_kernel_to_tof_status(st), "SS exact kernel failed");
    if (!detail::rat128_fits_int64(half))
        return detail::protocol_fail(
            e, TofStatus::Overflow, "ToF numerator/denominator do not fit int64 ticks");

    e.tof = detail::rat128_to_ticks(half);
    e.math_status = TofStatus::Ok;
    e.local_evidence_complete = true;
    e.peer_evidence_is_wire_claim = ra.is_peer_wire_claim() || db.is_peer_wire_claim();
    e.estimate_available = true;
    e.completion = (e.tof.num < 0) ? ProtocolCompletionStatus::Failed
                                   : ProtocolCompletionStatus::Complete;
    if (e.tof.num < 0) {
        e.estimate_available = false;
        e.math_status = TofStatus::NegativeTof;
        e.set_failure_reason(tof_status_to_exchange_status(TofStatus::NegativeTof));
        e.detail = "protocol SS ToF is negative (" + e.tof.to_string() +
                   " A ticks): retained with its sign";
    } else {
        e.detail = "protocol SS ToF = (RA - k*DB)/2 over wire/local inputs";
    }
    return e;
}


ProtocolTofEstimate compute_protocol_ds_tof(const ProtocolInterval& ra,
                                            const ProtocolInterval& rb,
                                            const ProtocolInterval& da,
                                            const ProtocolInterval& db,
                                            const ClockRatio& k_ab)
{
    ProtocolTofEstimate e;
    e.protocol = Protocol::Ds;
    e.computed_at = ComputedAt::ResponderB;
    e.ratio = k_ab;
    e.domain = k_ab.domain_a(); // the common unit is A ticks

    const TofStatus ks = k_ab.check();
    if (ks != TofStatus::Ok)
        return detail::protocol_fail(
            e, ks, "clock ratio is unusable: " + k_ab.to_string());

    if (!detail::protocol_domain_ok(e, "ra", ra, k_ab.domain_a()))
        return e;
    if (!detail::protocol_domain_ok(e, "da", da, k_ab.domain_a()))
        return e;
    if (!detail::protocol_domain_ok(e, "rb", rb, k_ab.domain_b()))
        return e;
    if (!detail::protocol_domain_ok(e, "db", db, k_ab.domain_b()))
        return e;
    if (!detail::protocol_window_ok(e, "ra", ra, k_ab))
        return e;
    if (!detail::protocol_window_ok(e, "da", da, k_ab))
        return e;

    detail::Rat128 RA;
    detail::Rat128 RB;
    detail::Rat128 DA;
    detail::Rat128 DB;
    if (!detail::protocol_load_interval(e, "ra", ra, RA))
        return e;
    if (!detail::protocol_load_interval(e, "rb", rb, RB))
        return e;
    if (!detail::protocol_load_interval(e, "da", da, DA))
        return e;
    if (!detail::protocol_load_interval(e, "db", db, DB))
        return e;

    const detail::Rat128 k = detail::rat128_make(k_ab.k_num(), k_ab.k_den());
    detail::Rat128 tof;
    const detail::KernelStatus st = detail::ds_tof_kernel(RA, RB, DA, DB, k, tof);
    if (st != detail::KernelStatus::Ok)
        return detail::protocol_fail(
            e, detail::protocol_kernel_to_tof_status(st), "DS exact kernel failed");
    if (!detail::rat128_fits_int64(tof))
        return detail::protocol_fail(
            e, TofStatus::Overflow, "ToF numerator/denominator do not fit int64 ticks");

    e.tof = detail::rat128_to_ticks(tof);
    e.math_status = TofStatus::Ok;
    e.local_evidence_complete = true;
    e.peer_evidence_is_wire_claim = ra.is_peer_wire_claim() || rb.is_peer_wire_claim() ||
                                    da.is_peer_wire_claim() || db.is_peer_wire_claim();
    e.estimate_available = true;
    e.completion = (e.tof.num < 0) ? ProtocolCompletionStatus::Failed
                                   : ProtocolCompletionStatus::Complete;
    if (e.tof.num < 0) {
        e.estimate_available = false;
        e.math_status = TofStatus::NegativeTof;
        e.set_failure_reason(tof_status_to_exchange_status(TofStatus::NegativeTof));
        e.detail = "protocol DS ToF is negative (" + e.tof.to_string() +
                   " A ticks): retained with its sign";
    } else {
        e.detail =
            "protocol DS ToF = (RA*k*RB - DA*k*DB)/(RA+k*RB+DA+k*DB), all in A ticks";
    }
    return e;
}


bool wire_timestamp_binding_is_well_formed(const WireTimestampBinding& b,
                                           std::string& why)
{
    if (!b.peer_domain.is_valid()) {
        why = "peer_domain is not a valid clock domain";
        return false;
    }
    if (!timestamp_marker_is_known(b.peer_marker)) {
        why = "peer_marker is not a known TimestampMarker";
        return false;
    }
    if (b.session_generation == 0) {
        why = "session_generation is 0 (unset)";
        return false;
    }
    if (b.binding_generation == 0) {
        why = "binding_generation is 0 (unset)";
        return false;
    }
    if (b.unit_convention.empty()) {
        why = "unit_convention is empty";
        return false;
    }
    if (b.calibration_convention_id.empty()) {
        why = "calibration_convention_id is empty";
        return false;
    }
    if (b.max_interval_ticks == 0) {
        why = "max_interval_ticks is 0";
        return false;
    }
    if (!sequence_modulus_is_supported(b.sequence_modulus)) {
        why = "sequence_modulus is not one of 4/16/64/256";
        return false;
    }
    const uint64_t period = clock_domain_wrap_period(b.peer_domain);
    if (period != 0u && b.max_interval_ticks > period / 2u) {
        why = "max_interval_ticks exceeds the domain's unambiguous half period";
        return false;
    }
    return true;
}


uint64_t frame_message_identity(const Frame& f)
{
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](uint64_t v) {
        for (int i = 0; i < 8; ++i) {
            h ^= (v >> (8 * i)) & 0xffu;
            h *= 1099511628211ull;
        }
    };
    mix(static_cast<uint64_t>(f.version));
    mix(static_cast<uint64_t>(f.function_code));
    mix(f.session_id);
    mix(f.seq);
    mix(f.pan_id);
    mix(f.src_addr);
    mix(f.dst_addr);
    mix(f.flags);
    for (size_t i = 0; i < f.timestamp_count(); ++i)
        mix(f.timestamps[i]);
    return h;
}


namespace detail {

ProtocolTofEstimate
protocol_fail(ProtocolTofEstimate e, TofStatus s, const std::string& why)
{
    e.completion = ProtocolCompletionStatus::Failed;
    e.estimate_available = false;
    e.math_status = s;
    e.set_failure_reason(tof_status_to_exchange_status(s));
    e.detail = why;
    return e;
}


bool protocol_domain_ok(ProtocolTofEstimate& e,
                        const char* what,
                        const ProtocolInterval& iv,
                        const ClockDomain& want)
{
    if (!iv.is_valid()) {
        e = protocol_fail(e,
                          TofStatus::InvalidInput,
                          std::string(what) + ": protocol interval is not valid");
        return false;
    }
    if (!protocol_interval_origin_is_known(iv.origin())) {
        e = protocol_fail(e,
                          TofStatus::InvalidInput,
                          std::string(what) +
                              ": protocol interval origin is out of domain");
        return false;
    }
    if (!iv.domain().is_valid()) {
        e = protocol_fail(e,
                          TofStatus::InvalidClockDomain,
                          std::string(what) +
                              ": interval carries an invalid clock domain");
        return false;
    }
    if (!clock_domain_same_identity(iv.domain(), want)) {
        e = protocol_fail(e,
                          TofStatus::ClockRatioDomainMismatch,
                          std::string(what) + " is in " + iv.domain().to_string() +
                              " but the clock ratio declares " + want.to_string());
        return false;
    }
    if (!clock_domain_same_epoch(iv.domain(), want)) {
        e = protocol_fail(
            e,
            TofStatus::ClockRatioEpochMismatch,
            std::string(what) + " is epoch " + std::to_string(iv.domain().epoch_id) +
                " but the clock ratio declares epoch " + std::to_string(want.epoch_id));
        return false;
    }
    return true;
}


bool protocol_load_interval(ProtocolTofEstimate& e,
                            const char* what,
                            const ProtocolInterval& iv,
                            Rat128& out)
{
    if (!iv.is_valid()) {
        e = protocol_fail(e,
                          TofStatus::InvalidInput,
                          std::string(what) + ": protocol interval is not valid");
        return false;
    }
    int64_t num = 0;
    int64_t den = 0;
    if (!iv.exact_ratio(num, den)) {
        e = protocol_fail(e,
                          TofStatus::Overflow,
                          std::string(what) + ": exact_ratio() did not fit int64");
        return false;
    }
    out = rat128_make(static_cast<i128>(num), static_cast<i128>(den));
    return true;
}


bool protocol_window_ok(ProtocolTofEstimate& e,
                        const char* what,
                        const ProtocolInterval& iv,
                        const ClockRatio& k)
{
    if (!k.has_validity_window())
        return true;
    if (!iv.is_valid()) {
        e = protocol_fail(e,
                          TofStatus::InvalidInput,
                          std::string(what) + ": protocol interval is not valid");
        return false;
    }
    const int64_t a = iv.earlier_tick();
    const int64_t b = iv.later_tick();
    const int64_t lo = (a < b) ? a : b;
    const int64_t hi = (a < b) ? b : a;
    if (!k.covers_ticks(lo) || !k.covers_ticks(hi)) {
        e = protocol_fail(e,
                          TofStatus::ClockRatioNotValidAtTime,
                          std::string(what) + " spans ticks [" + std::to_string(lo) +
                              "," + std::to_string(hi) +
                              "], outside the clock ratio validity window [" +
                              std::to_string(k.valid_from_ticks()) + "," +
                              std::to_string(k.valid_until_ticks()) + "]");
        return false;
    }
    return true;
}


} // namespace detail


RangeAdmission admit_range_capable_time(const Timestamp& ts,
                                        const RangeAdmissionContext& ctx)
{
    RangeAdmission r;

    // --- 1. REQ-TIME-01: the timestamp must be a well-formed instant. -----
    if (!timestamp_is_self_consistent(ts)) {
        detail::reject(r,
                       RangeAdmissionReason::InvalidTimestamp,
                       std::string("timestamp_rejection_reason=") +
                           timestamp_rejection_reason(ts),
                       RangeTimeSide::Subject,
                       ts.marker);
        return r;
    }

    // --- 2. REQ-TIME-02: the marker must denote an on-air instant. --------
    // Checked BEFORE the correction bits, and independently of them: a raw
    // capture coordinate is not an RMARKER however many stages claim to have
    // run, which is the R6 point.
    if (!detail::tof_marker_is_range_capable(ts.marker)) {
        detail::reject(
            r,
            RangeAdmissionReason::NotRangeCapableMarker,
            std::string("marker=") + timestamp_marker_to_string(ts.marker) + " class=" +
                timestamp_marker_class_to_string(timestamp_marker_class(ts.marker)) +
                " is not an on-air instant",
            RangeTimeSide::Subject,
            ts.marker);
        return r;
    }

    // --- 3. REQ-TIME-03: a real measurement, or an EVIDENCED schedule. ----
    //
    // This step used to refuse every source but `HardwareMeasured`, which also
    // refused the deterministic timed TX that REQ-TIME-03 explicitly PERMITS
    // ("可通过确定性时序加标定取得合格 TX timestamp，但必须标注
    // scheduled_calibrated 来源...不能冒充 hardware_measured").  The forbidden
    // thing is presenting a schedule AS a measurement, not using a calibrated
    // schedule -- and the planned X410 timed-TX route depends on the permitted
    // case, so refusing it would have forced the pipeline to mislabel its own
    // timestamps to make progress.
    //
    // The rule is therefore split by marker class:
    //   * a RECEIVE instant must be hardware-measured.  A scheduled receive
    //     time is a PREDICTION, which is the R6 point and stays refused;
    //   * a TRANSMIT RMARKER may be `ScheduledCalibrated`, but only with the
    //     send evidence REQ-TIME-03 lists.  It is never relabelled: the
    //     admitted value keeps source == ScheduledCalibrated, so a caller can
    //     never mistake a planned instant for a chip timestamp.
    //
    // A planned instant may be STORED before the exchange converges, but it
    // must not become a valid distance until the outcome is known; that is why
    // `Unknown` is a refusal here rather than a "probably fine".
    if (ts.source == TimestampSource::ScheduledCalibrated &&
        timestamp_marker_class(ts.marker) == TimestampMarkerClass::WaveformTx) {
        const TxSendEvidence* ev = ctx.tx_evidence;
        if (ev == nullptr) {
            detail::reject(
                r,
                RangeAdmissionReason::ScheduledTxEvidenceMissing,
                "a scheduled TX needs TxSendEvidence in the admission context: "
                "send() returning is not evidence that the burst happened as planned",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        const uint32_t missing = ev->missing_plan_record();
        if (missing != 0u) {
            uint32_t first = 0u;
            for (uint32_t bit = 1u; bit != 0u; bit <<= 1) {
                if ((missing & bit) != 0u) {
                    first = bit;
                    break;
                }
            }
            detail::reject(
                r,
                RangeAdmissionReason::ScheduledTxEvidenceMissing,
                std::string("the scheduled TX plan is not fully recorded: missing=") +
                    tx_plan_record_name(first) + " (mask=" + std::to_string(missing) +
                    ", " + ev->to_string() + "); send_accepted alone is never sufficient",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        // N07: the domain test comes FIRST.  `TxOutcome` is a PUBLIC enum, so a
        // variable -- set by a caller, a deserialiser or an async-event mapping
        // -- can hold a value the enum does not have.  The switch below covers
        // every enumerator, but it has no `default`, so an unknown value would
        // match no case, fall THROUGH it and reach the admitted path: fail-open
        // in the gate that decides whether a distance may be published.  The
        // domain test is the switch's partner, not a substitute for it.
        if (!tx_outcome_is_known(ev->outcome)) {
            detail::reject(
                r,
                RangeAdmissionReason::ScheduledTxEvidenceMissing,
                std::string("the TX outcome is not a member of TxOutcome (raw=") +
                    std::to_string(static_cast<int>(ev->outcome)) +
                    "): an out-of-domain outcome cannot be reasoned about, and "
                    "only Completed may publish a ranging instant",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        switch (ev->outcome) {
        case TxOutcome::Completed:
            break; // the burst completed as planned
        case TxOutcome::Late:
            detail::reject(r,
                           RangeAdmissionReason::TxOutcomeLate,
                           "the timed TX missed its deadline (ERR_TX_LATE_DATA): the "
                           "frame was not sent at the planned instant",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        case TxOutcome::Underflow:
            detail::reject(r,
                           RangeAdmissionReason::TxOutcomeUnderflow,
                           "the TX stream underran: the frame is not on the air",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        case TxOutcome::Cancelled:
            detail::reject(r,
                           RangeAdmissionReason::TxOutcomeCancelled,
                           "the TX burst was cancelled before it completed",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        case TxOutcome::Unknown:
            detail::reject(r,
                           RangeAdmissionReason::ScheduledTxEvidenceMissing,
                           std::string("the TX outcome is Unknown (") + ev->to_string() +
                               "): send_accepted does not establish that the burst "
                               "completed, and a planned instant is not a measurement "
                               "until it has",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
    } else if (ts.source != TimestampSource::HardwareMeasured) {
        const bool scheduled_rx = ts.source == TimestampSource::ScheduledCalibrated;
        detail::reject(
            r,
            RangeAdmissionReason::NotMeasured,
            std::string("source=") + timestamp_source_to_string(ts.source) +
                " marker=" + timestamp_marker_to_string(ts.marker) +
                (scheduled_rx
                     ? " is a scheduled RECEIVE time: a prediction, not a measured "
                       "receive instant"
                     : " is not a hardware measurement"),
            RangeTimeSide::Subject,
            ts.marker);
        return r;
    }

    // --- 4. REQ-TIME-02/03: every correction the marker requires. --------
    // The refusal names the FIRST missing stage, so the pipeline is told what
    // to run rather than just that something is wrong.
    if (!timestamp_corrections_satisfied(ts.marker, ts.applied_corrections)) {
        const uint32_t need = timestamp_required_corrections(ts.marker);
        uint32_t missing = 0u;
        for (uint32_t bit = 1u; bit != 0u; bit <<= 1) {
            if ((need & bit) != 0u && (ts.applied_corrections & bit) == 0u) {
                missing = bit;
                break;
            }
        }
        detail::reject(
            r,
            RangeAdmissionReason::MissingCorrection,
            std::string("missing=") + timestamp_correction_name(missing) +
                " marker=" + timestamp_marker_to_string(ts.marker) +
                " applied=" + timestamp_corrections_to_string(ts.applied_corrections),
            RangeTimeSide::Subject,
            ts.marker,
            missing);
        return r;
    }

    // --- 5. REQ-CAL-01: present, exactly once, current. ------------------
    {
        if (ts.calibration_id.empty()) {
            detail::reject(r,
                           RangeAdmissionReason::NoCalibration,
                           "timestamp carries no calibration_id",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        if (ctx.calibration == nullptr) {
            detail::reject(r,
                           RangeAdmissionReason::NoCalibration,
                           "no calibration set supplied in the admission context",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        const CalibrationStamp& cal = *ctx.calibration;
        if (cal.id.empty()) {
            detail::reject(r,
                           RangeAdmissionReason::NoCalibration,
                           "the calibration set in force has an empty id",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        // A zero-length or inverted window cannot show the calibration to be
        // current for any instant, so it is reported as STALENESS rather than
        // as an absent calibration: the set exists, it just does not cover
        // this exchange.
        if (!cal.is_well_formed() ||
            !detail::calibration_window_fits_domain(ts.domain, cal)) {
            detail::reject(
                r,
                RangeAdmissionReason::CalibrationStale,
                "calibration validity window is degenerate or does not fit the "
                "clock domain",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        if (!ctx.reference_ticks_recorded) {
            detail::reject(
                r,
                RangeAdmissionReason::NoCalibration,
                "no reference instant recorded: calibration currency cannot be decided",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        // Applied exactly once, by this set, with no conflict.
        if (cal.applications.size() != 1u) {
            detail::reject(
                r,
                cal.applications.empty() ? RangeAdmissionReason::NoCalibration
                                         : RangeAdmissionReason::CalibrationAppliedTwice,
                "calibration stage records " + std::to_string(cal.applications.size()) +
                    " applications; exactly one Applied is required",
                RangeTimeSide::Subject,
                ts.marker);
            return r;
        }
        const CalibrationApplication& app = cal.applications.front();
        if (app.result != CalibrationResult::Applied) {
            detail::reject(r,
                           RangeAdmissionReason::CalibrationAppliedTwice,
                           std::string("calibration stage returned ") +
                               calibration_result_to_string(app.result) +
                               "; the calibration was not applied exactly once",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        if (app.calibration_id != cal.id) {
            detail::reject(r,
                           RangeAdmissionReason::CalibrationAppliedTwice,
                           "calibration stage record names a different set: " +
                               app.calibration_id + " vs " + cal.id,
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        // The timestamp must be the one THIS set was applied to.
        if (ts.calibration_id != cal.id) {
            detail::reject(r,
                           RangeAdmissionReason::CalibrationMismatch,
                           "timestamp calibration_id=" + ts.calibration_id +
                               " but the exchange runs on " + cal.id,
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        // Currency: same epoch (a reboot invalidates the old constant) and
        // the exchange's own instant inside the window.
        if (cal.calibrated_epoch != ts.domain.epoch_id) {
            char buf[192];
            std::snprintf(buf,
                          sizeof(buf),
                          "calibration %s was established in epoch %llu, timestamp is in "
                          "epoch %llu",
                          cal.id.c_str(),
                          static_cast<unsigned long long>(cal.calibrated_epoch),
                          static_cast<unsigned long long>(ts.domain.epoch_id));
            detail::reject(r,
                           RangeAdmissionReason::CalibrationStale,
                           std::string(buf),
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        if (ctx.reference_ticks < cal.valid_from_ticks ||
            ctx.reference_ticks >= cal.valid_until_ticks) {
            char buf[256];
            std::snprintf(
                buf,
                sizeof(buf),
                "calibration %s valid on [%lld,%lld) but the exchange instant is "
                "%lld",
                cal.id.c_str(),
                static_cast<long long>(cal.valid_from_ticks),
                static_cast<long long>(cal.valid_until_ticks),
                static_cast<long long>(ctx.reference_ticks));
            detail::reject(r,
                           RangeAdmissionReason::CalibrationStale,
                           std::string(buf),
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
    }

    // --- 6. REQ-TIME-05: the first-path quality DECISION. ----------------
    // A default-constructed `FirstPathQuality` is `NotRecorded` and is
    // refused, so "we never ran the gate" cannot be read as "the gate passed".
    if (detail::tof_marker_requires_first_path(ts.marker)) {
        if (!ctx.rx_first_path.is_recorded()) {
            detail::reject(r,
                           RangeAdmissionReason::FirstPathQualityNotRecorded,
                           "no first-path quality decision recorded for marker " +
                               std::string(timestamp_marker_to_string(ts.marker)),
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        if (!ctx.rx_first_path.measurements_are_sane()) {
            detail::reject(r,
                           RangeAdmissionReason::FirstPathQualityInvalid,
                           "first-path quality decision carries unusable measurements",
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
        if (!ctx.rx_first_path.has_passed()) {
            detail::reject(r,
                           RangeAdmissionReason::FirstPathUnreliable,
                           "first-path quality gate failed: " +
                               ctx.rx_first_path.to_string(),
                           RangeTimeSide::Subject,
                           ts.marker);
            return r;
        }
    }

    // --- admitted --------------------------------------------------------
    r.admitted = true;
    r.reason = RangeAdmissionReason::Admitted;
    r.status = ExchangeStatus::Ok;
    r.side = RangeTimeSide::Subject;
    r.marker = ts.marker;
    r.missing_correction = 0u;
    r.detail.clear();
    r.value = RangeCapableTime(ts,
                               detail::tof_marker_requires_first_path(ts.marker)
                                   ? ctx.rx_first_path
                                   : FirstPathQuality::not_recorded());
    return r;
}


RangingIntervalAdmission admit_ranging_interval(const Timestamp& later,
                                                const Timestamp& earlier,
                                                const RangeAdmissionContext& ctx)
{
    RangingIntervalAdmission r;

    // --- 1. structural validity of both numbers --------------------------
    if (!timestamp_is_self_consistent(later) || !timestamp_is_self_consistent(earlier)) {
        const bool later_ok = timestamp_is_self_consistent(later);
        const Timestamp& bad = later_ok ? earlier : later;
        r.admitted = false;
        r.reason = RangeAdmissionReason::InvalidTimestamp;
        r.status = range_admission_reason_to_exchange_status(r.reason);
        r.side = later_ok ? RangeTimeSide::Earlier : RangeTimeSide::Later;
        r.marker = bad.marker;
        r.detail =
            std::string("timestamp_rejection_reason=") + timestamp_rejection_reason(bad);
        return r;
    }

    // --- 2. REQ-TIME-01: one timeline. ----------------------------------
    // Reuses the header's comparability rule (`ClockDomain` identity plus
    // epoch) -- there is deliberately no second notion of domain identity in
    // this file, and no cross-domain conversion is offered.
    if (!clock_domain_is_comparable(later.domain, earlier.domain)) {
        const bool same_epoch = clock_domain_same_epoch(later.domain, earlier.domain);
        r.admitted = false;
        r.side = RangeTimeSide::Earlier;
        r.marker = earlier.marker;
        r.reason = same_epoch ? RangeAdmissionReason::CrossDomain
                              : RangeAdmissionReason::WrongEpoch;
        r.status = range_admission_reason_to_exchange_status(r.reason);
        r.interval_status = same_epoch ? TimeIntervalStatus::DomainNameMismatch
                                       : TimeIntervalStatus::EpochMismatch;
        r.detail = same_epoch ? std::string("not one clock domain: ") +
                                    later.domain.to_string() + " vs " +
                                    earlier.domain.to_string()
                              : std::string("different epoch: ") +
                                    std::to_string(later.domain.epoch_id) + " vs " +
                                    std::to_string(earlier.domain.epoch_id);
        return r;
    }

    // --- 3. each side on its own ----------------------------------------
    const RangeAdmission a_later = admit_range_capable_time(later, ctx);
    if (!a_later.admitted) {
        r.admitted = false;
        r.reason = a_later.reason;
        r.status = a_later.status;
        r.side = RangeTimeSide::Later;
        r.marker = a_later.marker;
        r.missing_correction = a_later.missing_correction;
        r.detail = std::string("later: ") + a_later.detail;
        return r;
    }
    const RangeAdmission a_earlier = admit_range_capable_time(earlier, ctx);
    if (!a_earlier.admitted) {
        r.admitted = false;
        r.reason = a_earlier.reason;
        r.status = a_earlier.status;
        r.side = RangeTimeSide::Earlier;
        r.marker = a_earlier.marker;
        r.missing_correction = a_earlier.missing_correction;
        r.detail = std::string("earlier: ") + a_earlier.detail;
        return r;
    }

    // REQ-TIME-01: one timeline.  Reuses the header's comparability rule
    // (`ClockDomain` identity + epoch) -- there is deliberately no second
    // notion of domain identity in this file.
    if (!clock_domain_is_comparable(later.domain, earlier.domain)) {
        const bool same_epoch = clock_domain_same_epoch(later.domain, earlier.domain);
        r.admitted = false;
        r.side = RangeTimeSide::Earlier;
        r.marker = earlier.marker;
        r.reason = same_epoch ? RangeAdmissionReason::CrossDomain
                              : RangeAdmissionReason::WrongEpoch;
        r.status = range_admission_reason_to_exchange_status(r.reason);
        r.interval_status = same_epoch ? TimeIntervalStatus::DomainNameMismatch
                                       : TimeIntervalStatus::EpochMismatch;
        r.detail = same_epoch ? std::string("not one clock domain: ") +
                                    later.domain.to_string() + " vs " +
                                    earlier.domain.to_string()
                              : std::string("different epoch: ") +
                                    std::to_string(later.domain.epoch_id) + " vs " +
                                    std::to_string(earlier.domain.epoch_id);
        return r;
    }

    // --- 4. REQ-TIME-02: the marker PAIR must denote a protocol interval. -
    // Both sides can individually be on-air instants and still not be
    // subtractable (e.g. the antenna plane minus the protocol RMARKER is a
    // bookkeeping offset, not a propagation time).
    if (!timestamp_marker_interval_allowed(later.marker, earlier.marker)) {
        r.admitted = false;
        r.reason = RangeAdmissionReason::MarkerPairNotAllowed;
        r.status = range_admission_reason_to_exchange_status(r.reason);
        r.side = RangeTimeSide::Later;
        r.marker = later.marker;
        r.interval_status = TimeIntervalStatus::MarkerMismatch;
        r.detail = std::string("marker pair (") +
                   timestamp_marker_to_string(later.marker) + ", " +
                   timestamp_marker_to_string(earlier.marker) +
                   ") does not denote a protocol interval";
        return r;
    }

    // --- 5. REQ-TIME-04: the interval must actually be formable (wrap
    // ambiguity, modular order, exact sub-tick combination).
    // `timestamp_relative_interval()` applies the same gate as
    // `timestamp_interval()` and returns the exact tick-space value with no
    // nanosecond projection.
    RelativeTickInterval ri;
    TimeIntervalStatus istatus = TimeIntervalStatus::InvalidTimestamp;
    if (!timestamp_relative_interval(later, earlier, ri, istatus)) {
        r.admitted = false;
        r.reason = RangeAdmissionReason::IntervalNotFormable;
        r.status = range_admission_reason_to_exchange_status(r.reason);
        r.side = RangeTimeSide::Unknown;
        r.marker = later.marker;
        r.interval_status = istatus;
        r.detail = std::string("interval not formable: ") +
                   time_interval_status_to_string(istatus);
        return r;
    }

    r.admitted = true;
    r.reason = RangeAdmissionReason::Admitted;
    r.status = ExchangeStatus::Ok;
    r.side = RangeTimeSide::Unknown;
    r.marker = later.marker;
    r.missing_correction = 0u;
    r.interval_status = TimeIntervalStatus::Ok;
    r.detail.clear();
    r.value = AdmittedRangingInterval(ri, *a_later.value, *a_earlier.value);
    return r;
}


const char* range_admission_reason_to_string(RangeAdmissionReason r)
{
    switch (r) {
    case RangeAdmissionReason::Admitted:
        return "admitted";
    case RangeAdmissionReason::InvalidTimestamp:
        return "invalid_timestamp";
    case RangeAdmissionReason::NotRangeCapableMarker:
        return "not_range_capable_marker";
    case RangeAdmissionReason::NotMeasured:
        return "not_measured";
    case RangeAdmissionReason::MissingCorrection:
        return "missing_correction";
    case RangeAdmissionReason::NoCalibration:
        return "no_calibration";
    case RangeAdmissionReason::CalibrationMismatch:
        return "calibration_mismatch";
    case RangeAdmissionReason::CalibrationAppliedTwice:
        return "calibration_applied_twice";
    case RangeAdmissionReason::CalibrationStale:
        return "calibration_stale";
    case RangeAdmissionReason::FirstPathQualityNotRecorded:
        return "first_path_quality_not_recorded";
    case RangeAdmissionReason::FirstPathQualityInvalid:
        return "first_path_quality_invalid";
    case RangeAdmissionReason::FirstPathUnreliable:
        return "first_path_unreliable";
    case RangeAdmissionReason::MarkerPairNotAllowed:
        return "marker_pair_not_allowed";
    case RangeAdmissionReason::WrongEpoch:
        return "wrong_epoch";
    case RangeAdmissionReason::CrossDomain:
        return "cross_domain";
    case RangeAdmissionReason::IntervalNotFormable:
        return "interval_not_formable";
    case RangeAdmissionReason::ScheduledTxEvidenceMissing:
        return "scheduled_tx_evidence_missing";
    case RangeAdmissionReason::TxOutcomeLate:
        return "tx_outcome_late";
    case RangeAdmissionReason::TxOutcomeUnderflow:
        return "tx_outcome_underflow";
    case RangeAdmissionReason::TxOutcomeCancelled:
        return "tx_outcome_cancelled";
    }
    return "invalid";
}


ExchangeStatus range_admission_reason_to_exchange_status(RangeAdmissionReason r)
{
    switch (r) {
    case RangeAdmissionReason::Admitted:
        return ExchangeStatus::Ok;
    case RangeAdmissionReason::NoCalibration:
    case RangeAdmissionReason::CalibrationMismatch:
    case RangeAdmissionReason::CalibrationAppliedTwice:
        return ExchangeStatus::CalibrationMissing;
    case RangeAdmissionReason::CalibrationStale:
        return ExchangeStatus::CalibrationExpired;
    case RangeAdmissionReason::FirstPathQualityNotRecorded:
    case RangeAdmissionReason::FirstPathQualityInvalid:
    case RangeAdmissionReason::FirstPathUnreliable:
        return ExchangeStatus::FirstPathUnreliable;
    case RangeAdmissionReason::InvalidTimestamp:
    case RangeAdmissionReason::NotRangeCapableMarker:
    case RangeAdmissionReason::NotMeasured:
    case RangeAdmissionReason::MissingCorrection:
    case RangeAdmissionReason::MarkerPairNotAllowed:
    case RangeAdmissionReason::WrongEpoch:
    case RangeAdmissionReason::CrossDomain:
    case RangeAdmissionReason::IntervalNotFormable:
    case RangeAdmissionReason::ScheduledTxEvidenceMissing:
        return ExchangeStatus::InvalidTimeDomain;
    // A timed TX that did not complete is NOT an invalid time domain: the
    // domain is fine, the burst is the problem.  Each of these maps to its own
    // terminal status so a soak report can separate them without parsing a
    // detail string (REQ-ERR-01).
    case RangeAdmissionReason::TxOutcomeLate:
        return ExchangeStatus::TxLate;
    case RangeAdmissionReason::TxOutcomeUnderflow:
        return ExchangeStatus::TxUnderflow;
    case RangeAdmissionReason::TxOutcomeCancelled:
        return ExchangeStatus::Cancelled;
    }
    return ExchangeStatus::InternalError;
}


namespace detail {

void reject(RangeAdmission& r,
            RangeAdmissionReason reason,
            std::string detail,
            RangeTimeSide side,
            TimestampMarker marker,
            uint32_t missing)
{
    r.admitted = false;
    r.reason = reason;
    r.status = range_admission_reason_to_exchange_status(reason);
    r.value.reset();
    r.side = side;
    r.marker = marker;
    r.missing_correction = missing;
    r.detail = std::move(detail);
}


} // namespace detail

} // namespace twr
} // namespace openuwb
