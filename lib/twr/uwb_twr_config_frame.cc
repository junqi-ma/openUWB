// TWR configuration/frame implementation moved out of the public headers.
//
// Plan: docs/performance/多Subagent并行开发方案_按功能裁剪与编译测试再加速_2026-11-21.md
// (B5) DEC-A half: twr_timestamp / twr_config_effective / twr_config_validator /
// twr_config_util / twr_frame.  Short accessors stay inline in the headers.

#include <openuwb/twr/twr_config_effective.h>
#include <openuwb/twr/twr_config_util.h>
#include <openuwb/twr/twr_config_validator.h>
#include <openuwb/twr/twr_frame.h>
#include <openuwb/twr/twr_timestamp.h>

namespace openuwb {
namespace twr {

ValidationReport ConfigValidator::run(const TwrConfig& cfg) const
{
    d_report = ValidationReport{};
    // FIRST: an out-of-domain enum makes every later check meaningless --
    // a value the enum does not have cannot be reasoned about.  Catching it
    // here also means the refusal NAMES the enum field instead of surfacing
    // as a confusing complaint about some field that depends on it (an
    // unknown `session.protocol`, for instance, used to be reported as a
    // problem with `timing.response_to_final.ns`).
    check_enum_domains(cfg);
    check_meta(cfg);
    check_session(cfg);
    check_phy(cfg);
    check_frame(cfg);
    check_tx(cfg);
    check_rx(cfg);
    check_radio(cfg);
    check_timing(cfg);
    check_timeouts(cfg);
    check_calibration(cfg);
    check_diagnostics(cfg);
    check_frame_lengths(cfg);
    return d_report;
}


void ConfigValidator::check_enum_domains(const TwrConfig& c) const
{
    check_enum("session.protocol", protocol_is_known(c.session.protocol), "Protocol");
    check_enum("session.role", role_is_known(c.session.role), "Role");
    check_enum("phy.prf_class", prf_class_is_known(c.phy.prf_class), "PrfClass");
    check_enum("phy.data_rate", data_rate_is_known(c.phy.data_rate), "DataRate");
    check_enum("phy.phr_rate", phr_rate_is_known(c.phy.phr_rate), "PhrRate");
    check_enum("frame.frame_profile",
               frame_profile_id_is_supported(c.frame.frame_profile),
               "FrameProfileId");
    check_enum("frame.sfd_mode", sfd_mode_is_known(c.frame.sfd_mode), "SfdMode");
    check_enum("frame.phr_mode", phr_mode_is_known(c.frame.phr_mode), "PhrMode");
    check_enum(
        "frame.fcs_append", fcs_appender_is_known(c.frame.fcs_append), "FcsAppender");
    check_enum("frame.sts_mode", sts_mode_is_known(c.frame.sts_mode), "StsMode");
    check_enum(
        "tx.power_policy", tx_power_policy_is_known(c.tx.power_policy), "TxPowerPolicy");
    check_enum(
        "tx.pulse_shaping", pulse_shaping_is_known(c.tx.pulse_shaping), "PulseShaping");
    check_enum("rx.agc", agc_mode_is_known(c.rx.agc), "AgcMode");
    check_enum("calibration.link_delay_unit",
               time_unit_is_known(c.calibration.link_delay_unit),
               "TimeUnit");
    check_enum("calibration.first_path_algorithm",
               first_path_algorithm_is_known(c.calibration.first_path_algorithm),
               "FirstPathAlgorithm");
    check_enum("calibration.cfo_compensation",
               compensation_flag_is_known(c.calibration.cfo_compensation),
               "CompensationFlag");
    check_enum("calibration.sfo_compensation",
               compensation_flag_is_known(c.calibration.sfo_compensation),
               "CompensationFlag");

    // Every TimedField carries a (domain, reference) pair; both are enums,
    // and both are reachable by direct construction.
    check_timed_enums("frame.sfd_timeout", c.frame.sfd_timeout);
    check_timed_enums("timing.poll_start", c.timing.poll_start);
    check_timed_enums("timing.poll_to_response", c.timing.poll_to_response);
    check_timed_enums("timing.response_to_final", c.timing.response_to_final);
    check_timed_enums("timing.final_to_report", c.timing.final_to_report);
    check_timed_enums("timing.post_tx_rx_enable", c.timing.post_tx_rx_enable);
    check_timed_enums("timing.min_tx_lead_time", c.timing.min_tx_lead_time);
    check_timed_enums("timeouts.poll_rx_window", c.timeouts.poll_rx_window);
    check_timed_enums("timeouts.response_rx_window", c.timeouts.response_rx_window);
    check_timed_enums("timeouts.final_rx_window", c.timeouts.final_rx_window);
    check_timed_enums("timeouts.report_rx_window", c.timeouts.report_rx_window);
    check_timed_enums("timeouts.rx_timeout", c.timeouts.rx_timeout);
    check_timed_enums("timeouts.exchange_timeout", c.timeouts.exchange_timeout);
    check_timed_enums("timeouts.retry_interval", c.timeouts.retry_interval);
    check_timed_enums("diagnostics.stats_cadence", c.diagnostics.stats_cadence);

    // The peer bindings repeat one of the same enums, once per endpoint.
    for (size_t i = 0; i < c.radio.peers.size(); ++i) {
        const std::string p =
            "radio.peers[" + twr_int_to_text(static_cast<int64_t>(i)) + "]";
        check_enum(p + ".role", role_is_known(c.radio.peers[i].role), "Role");
    }
}


void ConfigValidator::check_timed_enums(const std::string& f, const TimedField& t) const
{
    check_enum(f + ".domain", time_domain_is_known(t.domain), "TimeDomain");
    check_enum(f + ".reference",
               time_reference_event_is_known(t.reference),
               "TimeReferenceEvent");
    // `marker` is OPTIONAL -- absent is a legal state (there is no RF marker
    // on a host-monotonic clock) -- but when it is present it must be a
    // value the enum has.  Missing this third enum was worth twelve more
    // accepted configurations (an exhaustive sweep over every enum-typed
    // field is what found it; covering domain+reference and assuming the
    // pair was complete was the mistake).
    if (t.marker.has())
        check_enum(f + ".marker",
                   timestamp_marker_is_known(t.marker.value()),
                   "TimestampMarker");
}


void ConfigValidator::check_meta(const TwrConfig& c) const
{
    if (c.meta.schema_version.empty())
        cfg_rej("meta.schema_version",
                ConfigReason::EmptyValue,
                "schema_version is required");
    if (c.meta.schema_version != d_caps.schema_version) {
        // M0.1: twr-config/2 is a BREAKING change -- frame.frame_profile
        // became required and phy.phr_rate changed value domain.  A v1
        // document is refused with a migration message, never silently
        // misread as a v2 one.
        if (c.meta.schema_version == "twr-config/1")
            unsup("meta.schema_version",
                  "twr-config/1 is superseded by twr-config/2: "
                  "frame.frame_profile is now REQUIRED (frame_v1) and "
                  "phy.phr_rate now takes 850k (the PHR rate the modulator "
                  "actually uses) rather than the payload rate; the frame "
                  "geometry moved from the old 7-byte header to the frozen "
                  "14-byte header, so on-air sizes are 16/26/31 B and "
                  "session.session_id must fit 16 bits. This document is "
                  "refused rather than reinterpreted");
        else
            unsup("meta.schema_version",
                  "this build implements " + d_caps.schema_version +
                      ", config asks for '" + c.meta.schema_version + "'");
    }
    if (c.meta.profile_version.empty())
        cfg_rej(
            "meta.profile_version",
            ConfigReason::EmptyValue,
            "profile_version is required (REQ-OUT-01 calibration/profile/config hash)");
}


void ConfigValidator::check_session(const TwrConfig& c) const
{
    const SessionConfig& s = c.session;
    if (s.local_address == s.peer_address)
        cfg_rej("session.peer_address",
                ConfigReason::FieldConflict,
                "local and peer address are identical (would make every frame a match)");
    if (s.local_address == 0x0000 || s.local_address == 0xffff)
        cfg_rej("session.local_address",
                ConfigReason::FieldConflict,
                "0x0000 and 0xffff are reserved/invalid local addresses");
    if (s.peer_address == 0x0000 || s.peer_address == 0xffff)
        cfg_rej("session.peer_address",
                ConfigReason::FieldConflict,
                "0x0000 and 0xffff are reserved/invalid peer addresses");
    if (s.pan_id == 0x0000 || s.pan_id == 0xffff)
        cfg_rej("session.pan_id",
                ConfigReason::EmptyValue,
                "PAN 0x0000/0xffff is not usable");
    if (s.session_id == 0)
        cfg_rej("session.session_id",
                ConfigReason::ZeroValue,
                "0 is the frame profile's 'no session' marker");
    else {
        // The wire field at kOffSessionId is 16 bits.  session_id_to_wire()
        // is identity-or-refuse: in range -> exact, out of range -> an
        // explicit error.  It is NOT folded, masked or hashed, because a
        // lossy map makes two different local sessions produce a
        // byte-identical session field and frame_match() would then accept
        // the wrong session's reply as this one's.
        uint16_t wire = 0;
        std::string sid_err_text;
        SessionIdError sid_err = SessionIdError::None;
        if (!session_id_to_wire(s.session_id, wire, sid_err_text, &sid_err)) {
            rej("session.session_id",
                ConfigReason::OutOfRange,
                session_id_error_to_exchange_status(sid_err),
                sid_err_text,
                "REQ-PROTO-01");
        }
    }
    if (s.sequence_modulus == 0)
        cfg_rej(
            "session.sequence_modulus", ConfigReason::ZeroValue, "wrap modulus is zero");
    else if (s.sequence_modulus != 4 && s.sequence_modulus != 16 &&
             s.sequence_modulus != 64 && s.sequence_modulus != 256)
        cfg_rej("session.sequence_modulus",
                ConfigReason::OutOfRange,
                "must be 4, 16, 64 or 256 (2/4/6/8 sequence bits)");
    if (static_cast<uint16_t>(s.sequence) >= s.sequence_modulus)
        cfg_rej("session.sequence",
                ConfigReason::IndexOutOfRange,
                "sequence does not fit the wrap modulus");
    if (s.measurement_count == 0 || s.measurement_count > kTwrMaxMeasurementCount)
        cfg_rej("session.measurement_count",
                ConfigReason::OutOfRange,
                "must be 1.." + twr_int_to_text(kTwrMaxMeasurementCount));
    if (s.measurement_interval.negative())
        cfg_rej("session.measurement_interval",
                ConfigReason::NegativeValue,
                "negative interval");
    if (s.max_attempts_per_exchange == 0 ||
        s.max_attempts_per_exchange > kTwrMaxAttemptsPerExchange)
        cfg_rej("session.max_attempts_per_exchange",
                ConfigReason::OutOfRange,
                "must be 1.." + twr_int_to_text(kTwrMaxAttemptsPerExchange));
    if (s.retry_backoff.negative())
        cfg_rej("session.retry_backoff", ConfigReason::NegativeValue, "negative backoff");
    if (s.max_attempts_per_exchange > 1) {
        if (s.retry_backoff.is_zero())
            cfg_rej("session.retry_backoff",
                    ConfigReason::FieldConflict,
                    "a retry policy with more than one attempt needs a non-zero backoff");
    } else if (!s.retry_backoff.is_zero()) {
        cfg_rej("session.retry_backoff",
                ConfigReason::FieldConflict,
                "retry_backoff is set but max_attempts_per_exchange is 1");
    }
    if (s.max_in_flight_exchanges != kTwrMaxInFlightExchanges)
        unsup("session.max_in_flight_exchanges",
              "phase 1 allows at most " + twr_int_to_text(kTwrMaxInFlightExchanges) +
                  " in-flight exchange per endpoint; requested " +
                  twr_int_to_text(s.max_in_flight_exchanges),
              ConfigReason::InFlightNotSupported,
              "REQ-API-01");
}


void ConfigValidator::check_phy(const TwrConfig& c) const
{
    const PhyConfig& p = c.phy;
    if (p.channel > 16u)
        cfg_rej("phy.channel", ConfigReason::IndexOutOfRange, "channel must be 0..16");
    else if (!d_caps.channel_supported(p.channel))
        unsup("phy.channel",
              "channel " + twr_int_to_text(p.channel) +
                  " is not in this build's capability whitelist; the phase-1 common "
                  "profile "
                  "is channel 5 (REQ-PHY-01/03) and other channels are per-channel "
                  "acceptance items");

    // Channel vs centre frequency: a conflict is rejected, never ignored.
    finite("phy.center_frequency_hz", p.center_frequency_hz, "REQ-PHY-01");
    if (std::isfinite(p.center_frequency_hz) && p.center_frequency_hz > 0.0) {
        const double plan = uwb_channel_center_frequency_hz(p.channel);
        if (plan <= 0.0) {
            cfg_rej("phy.center_frequency_hz",
                    ConfigReason::IndexOutOfRange,
                    "no channel plan entry for this channel");
        } else if (std::fabs(p.center_frequency_hz - plan) >
                   kTwrChannelFrequencyToleranceHz) {
            cfg_rej("phy.center_frequency_hz",
                    ConfigReason::ChannelFrequencyMismatch,
                    "channel " + twr_int_to_text(p.channel) + " is " +
                        twr_double_to_text(plan) + " Hz, config states " +
                        twr_double_to_text(p.center_frequency_hz) + " Hz",
                    "REQ-PHY-01");
        }
    }

    // code / PRF / channel / profile are cross-validated together.
    for (int dir = 0; dir < 2; ++dir) {
        const std::string f = dir == 0 ? "phy.tx_preamble_code" : "phy.rx_preamble_code";
        const uint8_t code = dir == 0 ? p.tx_preamble_code : p.rx_preamble_code;
        if (code == 0)
            cfg_rej(f, ConfigReason::ZeroValue, "preamble code is required");
        else if (!d_caps.code_index_supported(code))
            unsup(f,
                  "preamble code " + twr_int_to_text(code) +
                      " is not in the capability whitelist (" +
                      twr_int_to_text(kTwrCodeIndexMin) + ".." +
                      twr_int_to_text(kTwrCodeIndexMax) + ")");
    }
    if (p.tx_preamble_code != 0 && p.rx_preamble_code != 0 &&
        p.tx_preamble_code != p.rx_preamble_code)
        cfg_rej("phy.rx_preamble_code",
                ConfigReason::FieldConflict,
                "TX and RX preamble codes differ (" +
                    twr_int_to_text(p.tx_preamble_code) + " vs " +
                    twr_int_to_text(p.rx_preamble_code) + "); one profile uses one code",
                "REQ-PHY-01");
    if (!d_caps.prf_class_supported(p.prf_class))
        unsup("phy.prf_class",
              std::string("PRF class ") + prf_class_to_string(p.prf_class) +
                  " is not in the capability whitelist");
    // The PHR rate is validated on its OWN axis.  It is NOT compared to
    // the payload data rate: the PHR is a different modulation at a
    // different rate, and the data-rate field INSIDE the PHR describes
    // the payload, not the PHR's own transmission.  M0's `phr_rate ==
    // data_rate` rule therefore accepted the 6.81 Mb/s PHR this modulator
    // cannot produce and rejected the 0.85 Mb/s PHR it does produce.
    if (!phr_rate_is_implemented(p.phr_rate)) {
        unsup("phy.phr_rate",
              "PHR rate " + std::string(phr_rate_to_string(p.phr_rate)) +
                  " is not implemented: " + phr_rate_unsupported_reason(p.phr_rate) +
                  " It is NOT compared with phy.data_rate (" +
                  data_rate_to_string(p.data_rate) +
                  "): the PHR and the payload are separate modulations, and the "
                  "data-rate field carried inside the PHR describes the payload",
              ConfigReason::Unsupported,
              "REQ-PHY-02");
    }

    // MEASURED-PENDING: sync repetitions / SFD mode / data rate.
    if (!d_caps.data_rate_supported(p.data_rate))
        unsup("phy.data_rate",
              std::string("data rate ") + data_rate_to_string(p.data_rate) +
                  " is not in the capability whitelist: " + d_caps.pending_reason,
              ConfigReason::Unsupported,
              "REQ-PHY-01");
    if (p.preamble_symbols == 0)
        cfg_rej("phy.preamble_symbols",
                ConfigReason::ZeroValue,
                "preamble length is required");
    else {
        // Gate 1: the length must be a member of PreambleLength, i.e. one
        // this repository's demodulator / CFO path can actually run.  The
        // reason names THAT limit and never a chip's capability: the
        // Qorvo parts' own API accepts 64 through 2048, so "the hardware
        // cannot do 128" would be a false statement about hardware this
        // project has never spoken to.
        PreambleLength plen = PreambleLength::Sym64;
        if (!preamble_length_from_symbols(p.preamble_symbols, plen)) {
            unsup("phy.preamble_symbols",
                  preamble_length_unsupported_reason(p.preamble_symbols),
                  ConfigReason::Unsupported,
                  "REQ-PHY-01");
        }
        // Gate 2: the length must also be MEASURED for this (native rate,
        // code, SFD, PSDU, ranging) combination.  Separate from gate 1 on
        // purpose: a length can be implemented and still unmeasured.
        if (!d_caps.sync_repetitions_supported(p.preamble_symbols)) {
            // M0.1 R7: state the SCOPE of the refusal.  For the lengths
            // Qorvo parts do support, this build's decoder / CFO path is
            // the limit -- that is not a hardware claim, and the note
            // carries the vendor citation so C++ and the CSV agree.
            std::string note =
                evidence::preamble_length_supported_by_vendor(p.preamble_symbols)
                    ? std::string(" [") +
                          evidence::preamble_length_reject_note(p.preamble_symbols) + "]"
                    : std::string();
            unsup("phy.preamble_symbols",
                  "preamble length " + twr_int_to_text(p.preamble_symbols) +
                      " is not in the capability whitelist: " + d_caps.pending_reason +
                      note,
                  ConfigReason::Unsupported,
                  "REQ-PHY-01");
        }
    }
}


void ConfigValidator::check_frame(const TwrConfig& c) const
{
    const FrameFormatConfig& f = c.frame;

    // ---- which layout? -------------------------------------------------
    // The profile is an ENUM, so "an arbitrary non-empty profile string"
    // is not expressible.  An id this build has no implementation for is
    // refused with its name, never defaulted to frame v1.
    if (!frame_profile_id_is_supported(f.frame_profile)) {
        unsup(
            "frame.frame_profile",
            "frame profile " + std::string(frame_profile_id_to_string(f.frame_profile)) +
                " is not implemented: only the ids named by "
                "frame_profile_id_to_string() have a codec, a geometry authority and an "
                "encoder in this build. A free-form profile name is never an executable "
                "geometry",
            ConfigReason::Unsupported,
            "REQ-PROTO-06");
    }
    // The authority this config is checked against, and the codec profile
    // whose timestamp width it reports.  Both come from the enum; neither
    // is read out of the config.
    FrameProfileGeometry authority;
    FrameProfile profile;
    const bool have_profile = frame_profile_id_is_supported(f.frame_profile) &&
                              frame_geometry_for(f.frame_profile, authority) &&
                              frame_profile_for(f.frame_profile, profile);

    if (!d_caps.sfd_mode_supported(f.sfd_mode))
        unsup("frame.sfd_mode",
              std::string("SFD mode ") + sfd_mode_to_string(f.sfd_mode) +
                  " is not in the capability whitelist: " + d_caps.pending_reason,
              ConfigReason::Unsupported,
              "REQ-PHY-01");
    // The raw SFD length is a RESTATEMENT of the mode, never an
    // independent number: "4z2" and 4 symbols cannot both be stated.
    if (f.sfd_symbols == 0)
        cfg_rej("frame.sfd_symbols", ConfigReason::ZeroValue, "SFD length is required");
    else if (f.sfd_symbols != sfd_mode_symbols(f.sfd_mode))
        cfg_rej("frame.sfd_symbols",
                ConfigReason::FieldConflict,
                std::string("SFD mode ") + sfd_mode_to_string(f.sfd_mode) + " is " +
                    twr_int_to_text(sfd_mode_symbols(f.sfd_mode)) +
                    " symbols, config states " + twr_int_to_text(f.sfd_symbols) +
                    "; the length is derived from the mode, not chosen beside it",
                "REQ-PHY-01");

    // PHR presence/form.  Each unimplemented member carries its own reason
    // (phr_mode_unsupported_reason); M0 caught only `none` and let
    // `extended` through as a valid profile.
    if (f.phr_mode != PhrMode::Standard) {
        unsup("frame.phr_mode",
              std::string("PHR mode ") + phr_mode_to_string(f.phr_mode) +
                  " is not implemented: " + phr_mode_unsupported_reason(f.phr_mode),
              f.phr_mode == PhrMode::None ? ConfigReason::OutOfScope
                                          : ConfigReason::Unsupported,
              "REQ-PHY-01");
    }
    if (d_caps.ranging_bit_required && !f.ranging_bit)
        cfg_rej(
            "frame.ranging_bit",
            ConfigReason::FieldConflict,
            "a TWR MAC frame must set the ranging bit; ranging_bit=false is not a valid "
            "TWR profile",
            "REQ-PROTO-05");
    if (!d_caps.ranging_bit_required && f.ranging_bit)
        cfg_rej("frame.ranging_bit",
                ConfigReason::FieldConflict,
                "this build's whitelist carries no ranging frame",
                "REQ-PROTO-05");

    check_finite_timed("frame.sfd_timeout", f.sfd_timeout);
    validate_timed_field("frame.sfd_timeout",
                         f.sfd_timeout,
                         /*required=*/false,
                         TimeReferenceEvent::RxEnable,
                         /*require_marker=*/true);

    // Exactly one layer appends the FCS, and WHICH one is not free.
    // A config that names the MAC as the appender is DESCRIBABLE -- the
    // geometry authority has a variant for it, frame_geometry_check()
    // reports the mac_fcs_bytes that layout requires, and the on-air
    // length does not change -- but it is not ENCODABLE here:
    // FrameProfile::fcs_appended_by_modulation_layer is true, so
    // frame_profile_validate() rejects it, and
    // FrameProfileGeometry::executable() is false for it.  M0 accepted
    // such a configuration and then reported a "MAC PSDU" no encoder can
    // produce, which is the R3 class of defect one layer down.
    if (f.fcs_append == FcsAppender::None)
        cfg_rej("frame.fcs_append",
                ConfigReason::FieldConflict,
                "no layer appends the FCS; exactly one of mac/phy must");
    else if (f.fcs_append == FcsAppender::MacLayer && have_profile) {
        unsup(
            "frame.fcs_append",
            "fcs_append=mac asks for a MAC that builds and appends the FCS itself. The "
            "geometry is DESCRIBABLE and the on-air length is unchanged, but this codec "
            "only encodes the PHY-appends form: FrameProfileGeometry::executable() is "
            "false "
            "for the MAC-appends variant and frame_profile_validate() rejects "
            "fcs_appended_by_modulation_layer=false",
            ConfigReason::Unsupported,
            "REQ-API-01");
    }
    if (f.fcs_bytes == 0)
        cfg_rej("frame.fcs_bytes", ConfigReason::ZeroValue, "FCS length is required");
    else if (have_profile &&
             f.fcs_bytes != static_cast<uint16_t>(authority.fcs_bytes_on_air()))
        unsup("frame.fcs_bytes",
              "frame profile " +
                  std::string(frame_profile_id_to_string(f.frame_profile)) + " has a " +
                  twr_int_to_text(authority.fcs_bytes_on_air()) +
                  "-byte FCS (IEEE 802.15.4), the config states " +
                  twr_int_to_text(f.fcs_bytes) + " bytes",
              ConfigReason::Unsupported,
              "REQ-API-01");
    const bool declared_includes = (f.fcs_append == FcsAppender::MacLayer);
    if (f.mac_psdu_includes_fcs != declared_includes)
        cfg_rej("frame.mac_psdu_includes_fcs",
                ConfigReason::FieldConflict,
                std::string("fcs_append=") + fcs_appender_to_string(f.fcs_append) +
                    " but mac_psdu_includes_fcs=" +
                    twr_bool_to_text(f.mac_psdu_includes_fcs) +
                    "; exactly one layer appends the FCS",
                "REQ-API-01");

    // ---- the geometry CLAIM, field by field, against the codec ---------
    if (!f.geometry.any_set()) {
        cfg_rej(
            "frame.geometry",
            ConfigReason::EmptyValue,
            "the frame profile geometry must be stated explicitly; it is never assumed. "
            "Fill it with frame_geometry_of_profile(frame.frame_profile) so the claim "
            "comes from the codec rather than from a literal");
    } else if (have_profile) {
        // A claim that disagrees with the codec is a REJECTION PER FIELD,
        // not a silently-corrected value: an operator who believes the
        // header is 7 bytes must be told it is 14, not handed a working
        // config that says something else.
        //
        // The authority the claim is checked against follows the FCS owner
        // the config NAMES, not the one it should name: a MAC-appends
        // layout is described by frame_geometry_mac_appends_fcs(), which
        // expects the 2 FCS bytes inside the MAC PSDU.  That way a
        // MAC-appends claim that forgets to reserve them is caught as the
        // field disagreement it is, and the separate `fcs_append`
        // rejection above says the layout is not encodable here at all.
        const FrameProfileGeometry& named_authority =
            (f.fcs_append == FcsAppender::MacLayer)
                ? frame_geometry_mac_appends_fcs(f.frame_profile)
                : authority;
        const GeometryCheckResult geo = frame_geometry_check(
            frame_geometry_claim(f.geometry), named_authority, profile);
        for (size_t i = 0; i < geo.fields.count; ++i) {
            const GeometryMismatch& m = geo.fields.items[i];
            cfg_rej(std::string("frame.geometry.") + geometry_field_to_string(m.field),
                    ConfigReason::FieldConflict,
                    std::string("the frame codec is the geometry authority: profile ") +
                        frame_profile_id_to_string(f.frame_profile) +
                        " with FCS appended by " +
                        fcs_owner_to_string(named_authority.fcs_owner) + " has " +
                        geometry_field_to_string(m.field) + " = " +
                        twr_int_to_text(static_cast<int64_t>(m.expected)) +
                        ", the config claims " +
                        twr_int_to_text(static_cast<int64_t>(m.actual)) +
                        ". A claim does not decide the layout: fill the geometry from "
                        "frame_geometry_of_profile() so there is one source",
                    "REQ-PROTO-06");
        }
    }

    // STS: out of scope for phase 1, explicitly.
    if (f.sts_mode != StsMode::Off)
        unsup("frame.sts_mode",
              std::string("STS mode ") + sts_mode_to_string(f.sts_mode) + ": " +
                  d_caps.unsupported_sts_reason,
              ConfigReason::OutOfScope,
              "REQ-SCOPE-04");
    if (f.sts_length_symbols != 0)
        cfg_rej("frame.sts_length_symbols",
                ConfigReason::FieldConflict,
                "sts_mode is off but an STS length is requested");
}


void ConfigValidator::check_tx(const TwrConfig& c) const
{
    const TransmitConfig& t = c.tx;
    if (static_cast<uint32_t>(t.port) >= kTwrMaxPhysicalChannels)
        cfg_rej("tx.port",
                ConfigReason::IndexOutOfRange,
                "physical TX channel index too large");
    check_finite_opt("tx.gain_db", t.gain_db);
    check_finite_opt("tx.iq_amplitude", t.iq_amplitude);
    check_finite_opt("tx.calibrated_tx_power_dbm", t.calibrated_tx_power_dbm);

    if (t.gain_db.has() && (t.gain_db.value() < 0.0 || t.gain_db.value() > 120.0))
        cfg_rej(
            "tx.gain_db", ConfigReason::OutOfRange, "TX gain must be within 0..120 dB");
    if (t.iq_amplitude.has() &&
        (t.iq_amplitude.value() <= 0.0 || t.iq_amplitude.value() > 1.0))
        cfg_rej(
            "tx.iq_amplitude",
            ConfigReason::OutOfRange,
            "digital IQ amplitude must be within (0, 1] for a full-scale SC16 stream");
    if (t.calibrated_tx_power_dbm.has() && (t.calibrated_tx_power_dbm.value() < -30.0 ||
                                            t.calibrated_tx_power_dbm.value() > 30.0))
        cfg_rej("tx.calibrated_tx_power_dbm",
                ConfigReason::OutOfRange,
                "calibrated antenna-plane power must be within -30..30 dBm");

    // Exactly one power authority, and the authority must be the field that
    // is actually present.  The three quantities are never conflated.
    int present = 0;
    if (t.gain_db.has())
        present++;
    if (t.iq_amplitude.has())
        present++;
    if (t.calibrated_tx_power_dbm.has())
        present++;
    const bool authority_present =
        (t.power_policy == TxPowerPolicy::LeaveUntouched) ||
        (t.power_policy == TxPowerPolicy::ManualGainDb && t.gain_db.has()) ||
        (t.power_policy == TxPowerPolicy::IqAmplitude && t.iq_amplitude.has()) ||
        (t.power_policy == TxPowerPolicy::CalibratedDbm &&
         t.calibrated_tx_power_dbm.has());
    if (!authority_present)
        cfg_rej("tx.power_policy",
                ConfigReason::FieldConflict,
                std::string("power_policy=") + tx_power_policy_to_string(t.power_policy) +
                    " but the field it names is not set (gain dB, IQ amplitude and "
                    "calibrated "
                    "dBm are three separate fields)");
    if (t.power_policy == TxPowerPolicy::LeaveUntouched && present != 0)
        cfg_rej("tx.power_policy",
                ConfigReason::FieldConflict,
                "power_policy is leave_untouched but " + twr_int_to_text(present) +
                    " power field(s) are set");
    if (t.power_policy == TxPowerPolicy::CalibratedDbm &&
        t.calibrated_tx_power_dbm.has() && !t.gain_db.has())
        cfg_rej("tx.gain_db",
                ConfigReason::FieldConflict,
                "a calibrated dBm is only meaningful together with the gain it was "
                "measured at");

    // A vendor power word is a fourth, separately named field and may only
    // exist with the backend capability that defines its bits.
    if (t.vendor_power_word.has()) {
        if (t.vendor_power_word_backend.empty())
            cfg_rej("tx.vendor_power_word",
                    ConfigReason::CalibrationMismatch,
                    "a vendor power word is meaningless without the backend capability "
                    "id that "
                    "defines its bits");
    } else if (!t.vendor_power_word_backend.empty()) {
        cfg_rej("tx.vendor_power_word_backend",
                ConfigReason::FieldConflict,
                "a backend capability id is set but no vendor power word is");
    }
}


void ConfigValidator::check_rx(const TwrConfig& c) const
{
    const ReceiveConfig& r = c.rx;
    if (static_cast<uint32_t>(r.port) >= kTwrMaxPhysicalChannels)
        cfg_rej("rx.port",
                ConfigReason::IndexOutOfRange,
                "physical RX channel index too large");
    check_finite_opt("rx.gain_db", r.gain_db);
    check_finite_opt("rx.bandwidth_hz", r.bandwidth_hz);
    finite("rx.detection_threshold", r.detection_threshold);
    finite("rx.correlation_threshold", r.correlation_threshold);
    finite("rx.first_path_threshold", r.first_path_threshold);
    check_finite_opt("rx.vendor_pac_applied_step", r.vendor_pac_applied_step);

    if (r.gain_db.has() && (r.gain_db.value() < 0.0 || r.gain_db.value() > 120.0))
        cfg_rej(
            "rx.gain_db", ConfigReason::OutOfRange, "RX gain must be within 0..120 dB");
    if (r.agc == AgcMode::Manual && !r.gain_db.has())
        cfg_rej("rx.gain_db",
                ConfigReason::FieldConflict,
                "agc=manual needs an explicit RX gain in dB");
    if (r.bandwidth_hz.has() && r.bandwidth_hz.value() <= 0.0)
        cfg_rej(
            "rx.bandwidth_hz",
            ConfigReason::OutOfRange,
            "RF bandwidth must be > 0; it is a different quantity from the sample rate "
            "(REQ-PHY-02)");
    if (r.bandwidth_hz.has() && r.bandwidth_hz.value() > 1.0e9)
        cfg_rej(
            "rx.bandwidth_hz", ConfigReason::OutOfRange, "RF bandwidth is implausible");
    if (r.detection_threshold < 0.0 || r.detection_threshold > 1.0)
        cfg_rej(
            "rx.detection_threshold", ConfigReason::OutOfRange, "must be within 0..1");
    if (r.correlation_threshold < 0.0 || r.correlation_threshold > 1.0)
        cfg_rej(
            "rx.correlation_threshold", ConfigReason::OutOfRange, "must be within 0..1");
    if (r.first_path_threshold < 0.0 || r.first_path_threshold > 1.0)
        cfg_rej(
            "rx.first_path_threshold", ConfigReason::OutOfRange, "must be within 0..1");
    if (r.correlation_threshold < r.detection_threshold)
        cfg_rej("rx.correlation_threshold",
                ConfigReason::FieldConflict,
                "the fine correlation gate must not be looser than the coarse detection "
                "gate");
    if (r.first_path_window == 0)
        cfg_rej("rx.first_path_window",
                ConfigReason::ZeroValue,
                "the first-path search window must be at least one tap");
    else if (r.first_path_index >= r.first_path_window)
        cfg_rej("rx.first_path_index",
                ConfigReason::IndexOutOfRange,
                "first usable CIR index lies outside the search window");

    // A vendor PAC value must come from a backend capability query AND
    // carry the software step that is actually applied.  Mapping a PAC
    // value onto an arbitrary software detection step is forbidden.
    if (r.vendor_pac_value.has()) {
        if (r.vendor_pac_backend.empty())
            cfg_rej("rx.vendor_pac_value",
                    ConfigReason::CalibrationMismatch,
                    "a vendor PAC value is meaningless without the backend capability id "
                    "that "
                    "defines it");
        if (!r.vendor_pac_applied_step.has())
            cfg_rej("rx.vendor_pac_applied_step",
                    ConfigReason::CalibrationMismatch,
                    "a vendor PAC value may not be mapped onto an arbitrary software "
                    "detection step; the backend must state the step it applies");
        else if (r.vendor_pac_applied_step.value() <= 0.0)
            cfg_rej("rx.vendor_pac_applied_step",
                    ConfigReason::OutOfRange,
                    "the applied software step must be > 0");
    } else {
        if (!r.vendor_pac_backend.empty())
            cfg_rej("rx.vendor_pac_backend",
                    ConfigReason::FieldConflict,
                    "a backend capability id is set but no vendor PAC value is");
        if (r.vendor_pac_applied_step.has())
            cfg_rej("rx.vendor_pac_applied_step",
                    ConfigReason::FieldConflict,
                    "a software step is declared without a vendor PAC value");
    }
}


void ConfigValidator::check_radio(const TwrConfig& c) const
{
    const RadioConfig& r = c.radio;
    if (static_cast<uint32_t>(r.tx_channel) >= kTwrMaxPhysicalChannels)
        cfg_rej("radio.tx_channel",
                ConfigReason::IndexOutOfRange,
                "TX channel index too large");
    if (static_cast<uint32_t>(r.rx_channel) >= kTwrMaxPhysicalChannels)
        cfg_rej("radio.rx_channel",
                ConfigReason::IndexOutOfRange,
                "RX channel index too large");
    finite("radio.native_sample_rate_hz", r.native_sample_rate_hz, "REQ-PHY-02");
    if (std::isfinite(r.native_sample_rate_hz) && r.native_sample_rate_hz > 0.0 &&
        !d_caps.native_rate_supported(r.native_sample_rate_hz))
        unsup("radio.native_sample_rate_hz",
              "native sample rate " + twr_double_to_text(r.native_sample_rate_hz) +
                  " is not supported; this build accepts 737280000 or 491520000",
              ConfigReason::Unsupported,
              "REQ-PHY-02");

    if (r.peers.size() > kTwrMaxPeersPerConfig)
        cfg_rej("radio.peers",
                ConfigReason::OverCapacity,
                "at most " +
                    twr_int_to_text(static_cast<int64_t>(kTwrMaxPeersPerConfig)) +
                    " peer endpoints per config");

    // Resource conflicts on the same physical channel.
    for (size_t i = 0; i < r.peers.size(); ++i) {
        const EndpointBinding& p = r.peers[i];
        const std::string pf =
            "radio.peers[" + twr_int_to_text(static_cast<int64_t>(i)) + "]";
        if (p.id.empty())
            cfg_rej(pf + ".id", ConfigReason::EmptyValue, "endpoint id is required");
        if (static_cast<uint32_t>(p.tx_channel) >= kTwrMaxPhysicalChannels)
            cfg_rej(pf + ".tx_channel",
                    ConfigReason::IndexOutOfRange,
                    "TX channel too large");
        if (static_cast<uint32_t>(p.rx_channel) >= kTwrMaxPhysicalChannels)
            cfg_rej(pf + ".rx_channel",
                    ConfigReason::IndexOutOfRange,
                    "RX channel too large");
        finite(pf + ".native_sample_rate_hz", p.native_sample_rate_hz);
        if (!p.occupies_resources)
            continue;
        if (p.tx_channel == r.tx_channel)
            cfg_rej(pf + ".tx_channel",
                    ConfigReason::DuplicateResource,
                    "physical TX channel " + twr_int_to_text(p.tx_channel) +
                        " is already claimed by this endpoint; one device resource, one "
                        "owner",
                    "REQ-BASE-01");
        if (p.rx_channel == r.rx_channel)
            cfg_rej(pf + ".rx_channel",
                    ConfigReason::DuplicateResource,
                    "physical RX channel " + twr_int_to_text(p.rx_channel) +
                        " is already claimed by this endpoint; one device resource, one "
                        "owner",
                    "REQ-BASE-01");
        if (p.role == c.session.role)
            cfg_rej(pf + ".role",
                    ConfigReason::FieldConflict,
                    std::string("both endpoints claim role ") + role_to_string(p.role) +
                        "; a two-endpoint session needs one initiator and one responder",
                    "REQ-PROTO-01");
        if (p.native_sample_rate_hz > 0.0 && r.native_sample_rate_hz > 0.0 &&
            !twr_rate_matches_strict(
                p.native_sample_rate_hz, r.native_sample_rate_hz, 1e-9))
            cfg_rej(pf + ".native_sample_rate_hz",
                    ConfigReason::FieldConflict,
                    "endpoints sharing one device must share the native sample rate (" +
                        twr_double_to_text(r.native_sample_rate_hz) + " vs " +
                        twr_double_to_text(p.native_sample_rate_hz) + ")",
                    "REQ-PHY-02");
        for (size_t j = 0; j < i; ++j) {
            const EndpointBinding& q = r.peers[j];
            if (!q.occupies_resources)
                continue;
            if (q.id == p.id)
                cfg_rej(
                    pf + ".id", ConfigReason::DuplicateResource, "duplicate endpoint id");
            if (q.tx_channel == p.tx_channel)
                cfg_rej(pf + ".tx_channel",
                        ConfigReason::DuplicateResource,
                        "physical TX channel " + twr_int_to_text(p.tx_channel) +
                            " claimed by two endpoints",
                        "REQ-BASE-01");
            if (q.rx_channel == p.rx_channel)
                cfg_rej(pf + ".rx_channel",
                        ConfigReason::DuplicateResource,
                        "physical RX channel " + twr_int_to_text(p.rx_channel) +
                            " claimed by two endpoints",
                        "REQ-BASE-01");
            if (q.role == p.role)
                cfg_rej(pf + ".role",
                        ConfigReason::FieldConflict,
                        std::string("two peers claim the same role ") +
                            role_to_string(p.role),
                        "REQ-PROTO-01");
        }
    }

    // Readback: a mismatch is a STARTUP FAILURE, never a soft warning.
    if (r.require_readback && !r.readback.present)
        rej("radio.readback.present",
            ConfigReason::CalibrationMissing,
            ExchangeStatus::CalibrationMissing,
            "require_readback is set but no hardware readback was recorded",
            "REQ-PHY-02");
    if (r.readback.present) {
        finite("radio.readback.sample_rate_hz", r.readback.sample_rate_hz, "REQ-PHY-02");
        finite("radio.readback.center_freq_hz", r.readback.center_freq_hz, "REQ-PHY-01");
        if (!d_caps.native_rate_supported(r.readback.sample_rate_hz))
            cfg_rej("radio.readback.sample_rate_hz",
                    ConfigReason::OutOfRange,
                    "readback native rate " +
                        twr_double_to_text(r.readback.sample_rate_hz) +
                        " Hz is not a supported rate");
        else if (r.native_sample_rate_hz > 0.0 &&
                 !twr_rate_matches_strict(
                     r.native_sample_rate_hz, r.readback.sample_rate_hz, 1e-9))
            cfg_rej("radio.readback.sample_rate_hz",
                    ConfigReason::FieldConflict,
                    "readback native rate " +
                        twr_double_to_text(r.readback.sample_rate_hz) +
                        " Hz differs from the requested " +
                        twr_double_to_text(r.native_sample_rate_hz) + " Hz",
                    "REQ-PHY-02");
        if (r.readback.tx_channel != r.tx_channel)
            cfg_rej("radio.readback.tx_channel",
                    ConfigReason::FieldConflict,
                    "readback TX channel differs from the requested one",
                    "REQ-PHY-02");
        if (r.readback.rx_channel != r.rx_channel)
            cfg_rej("radio.readback.rx_channel",
                    ConfigReason::FieldConflict,
                    "readback RX channel differs from the requested one",
                    "REQ-PHY-02");
        if (c.phy.center_frequency_hz > 0.0 && r.readback.center_freq_hz > 0.0 &&
            std::fabs(c.phy.center_frequency_hz - r.readback.center_freq_hz) >
                kTwrChannelFrequencyToleranceHz)
            cfg_rej("radio.readback.center_freq_hz",
                    ConfigReason::ChannelFrequencyMismatch,
                    "readback centre frequency " +
                        twr_double_to_text(r.readback.center_freq_hz) +
                        " Hz differs from the requested " +
                        twr_double_to_text(c.phy.center_frequency_hz) + " Hz",
                    "REQ-PHY-01");
        if (!r.clock_source.empty() && !r.readback.clock_source.empty() &&
            r.clock_source != r.readback.clock_source)
            cfg_rej("radio.readback.clock_source",
                    ConfigReason::FieldConflict,
                    "readback clock source differs from the requested one",
                    "REQ-PHY-02");
        if (!r.time_source.empty() && !r.readback.time_source.empty() &&
            r.time_source != r.readback.time_source)
            cfg_rej("radio.readback.time_source",
                    ConfigReason::FieldConflict,
                    "readback time source differs from the requested one",
                    "REQ-PHY-02");
    }
}


void ConfigValidator::check_timing(const TwrConfig& c) const
{
    const PerMessageTiming& t = c.timing;
    const bool ss = c.session.protocol == Protocol::Ss;
    const bool initiator = c.session.role == Role::Initiator;

    check_finite_timed("timing.poll_start", t.poll_start);
    check_finite_timed("timing.poll_to_response", t.poll_to_response);
    check_finite_timed("timing.response_to_final", t.response_to_final);
    check_finite_timed("timing.final_to_report", t.final_to_report);
    check_finite_timed("timing.post_tx_rx_enable", t.post_tx_rx_enable);
    check_finite_timed("timing.min_tx_lead_time", t.min_tx_lead_time);

    // poll_start: absolute, on the device clock, at the Poll TX RMARKER.
    // Only the initiator transmits a Poll, so a responder's poll_start must
    // be zero rather than silently reusing the initiator's value.
    validate_timed_field("timing.poll_start",
                         t.poll_start,
                         initiator,
                         TimeReferenceEvent::PollTransmitRmarker,
                         true);
    if (!initiator && !t.poll_start.value.is_zero())
        cfg_rej("timing.poll_start.ns",
                ConfigReason::FieldConflict,
                "the responder transmits no Poll, so the Poll start time must be zero",
                "REQ-API-02");

    // The three distinct quantities.
    if (initiator) {
        // The initiator does not send a Response, so it owns no reply delay.
        if (!t.poll_to_response.value.is_zero())
            cfg_rej(
                "timing.poll_to_response",
                ConfigReason::FieldConflict,
                "the initiator never sends a Response, so this reply delay must be zero",
                "REQ-API-02");
        if (t.poll_to_response.reference != TimeReferenceEvent::PollReceiveRmarker)
            cfg_rej("timing.poll_to_response.reference",
                    ConfigReason::FieldConflict,
                    "must be poll_rx_rmarker",
                    "REQ-API-02");
        if (ss) {
            if (!t.response_to_final.value.is_zero())
                cfg_rej("timing.response_to_final",
                        ConfigReason::FieldConflict,
                        "SS-TWR has no Final message, so this delay must be zero",
                        "REQ-API-02");
        } else {
            validate_timed_field("timing.response_to_final",
                                 t.response_to_final,
                                 true,
                                 TimeReferenceEvent::ResponseReceiveRmarker,
                                 true);
        }
    } else {
        if (!t.response_to_final.value.is_zero())
            cfg_rej("timing.response_to_final",
                    ConfigReason::FieldConflict,
                    "the responder never sends a Final, so this delay must be zero",
                    "REQ-API-02");
        if (t.response_to_final.reference != TimeReferenceEvent::FinalTransmitRmarker)
            cfg_rej("timing.response_to_final.reference",
                    ConfigReason::FieldConflict,
                    "must be final_tx_rmarker",
                    "REQ-API-02");
        validate_timed_field("timing.poll_to_response",
                             t.poll_to_response,
                             true,
                             TimeReferenceEvent::PollReceiveRmarker,
                             true);
    }

    // post-TX RX enable: a device-tick delay from the end of OUR OWN frame.
    validate_timed_field("timing.post_tx_rx_enable",
                         t.post_tx_rx_enable,
                         true,
                         TimeReferenceEvent::FrameTail,
                         true);
    validate_timed_field("timing.min_tx_lead_time",
                         t.min_tx_lead_time,
                         true,
                         TimeReferenceEvent::ResponseTransmitRmarker,
                         true);

    // Report is reserved in phase 1.
    if (!t.final_to_report.value.is_zero())
        unsup("timing.final_to_report",
              "FrameType::Report is reserved in phase 1; a Report message would "
              "invalidate the "
              "three-message interop claim (REQ-PROTO-05)",
              ConfigReason::OutOfScope,
              "REQ-PROTO-05");

    // REQ-GR-04 feasibility: the reply delay must cover the device lead
    // time.  There is no default lead time -- the operator must supply the
    // measured value, and a delay shorter than it is rejected instead of
    // being sent late.
    const TimedField& reply = initiator ? t.response_to_final : t.poll_to_response;
    if (reply.value.negative() || t.min_tx_lead_time.value.negative() ||
        reply.value.is_zero())
        return; // the structural violation was already reported above
    if (reply.value > t.min_tx_lead_time.value)
        return; // the budget is feasible
    rej(initiator ? "timing.response_to_final.ns" : "timing.poll_to_response.ns",
        ConfigReason::TimingBudgetInfeasible,
        ExchangeStatus::DeadlineMissed,
        "reply delay " + twr_duration_to_text(reply.value) +
            " does not cover the measured minimum TX lead time " +
            twr_duration_to_text(t.min_tx_lead_time.value) +
            " (a scheduled TX whose deadline cannot be met must fail, not be sent late)",
        "REQ-GR-04");
}


void ConfigValidator::check_timeouts(const TwrConfig& c) const
{
    const TimeoutConfig& t = c.timeouts;
    const bool ss = c.session.protocol == Protocol::Ss;
    const bool initiator = c.session.role == Role::Initiator;

    check_finite_timed("timeouts.poll_rx_window", t.poll_rx_window);
    check_finite_timed("timeouts.response_rx_window", t.response_rx_window);
    check_finite_timed("timeouts.final_rx_window", t.final_rx_window);
    check_finite_timed("timeouts.report_rx_window", t.report_rx_window);
    check_finite_timed("timeouts.rx_timeout", t.rx_timeout);
    check_finite_timed("timeouts.exchange_timeout", t.exchange_timeout);
    check_finite_timed("timeouts.retry_interval", t.retry_interval);

    // The RX timeout is measured FROM RX ENABLE, on the device clock.
    validate_timed_field(
        "timeouts.rx_timeout", t.rx_timeout, true, TimeReferenceEvent::RxEnable, true);

    // The whole-exchange deadline and the retry gap live on the host
    // monotonic clock and must carry NO RF marker (REQ-TIME-01).
    validate_timed_field("timeouts.exchange_timeout",
                         t.exchange_timeout,
                         true,
                         TimeReferenceEvent::HostMonotonic,
                         false);
    validate_timed_field("timeouts.retry_interval",
                         t.retry_interval,
                         false,
                         TimeReferenceEvent::HostMonotonic,
                         false);

    // Per-frame RX windows, armed at RX enable.  Only the frames this role
    // actually waits for may be non-zero.
    validate_timed_field("timeouts.poll_rx_window",
                         t.poll_rx_window,
                         !initiator,
                         TimeReferenceEvent::RxEnable,
                         true);
    validate_timed_field("timeouts.response_rx_window",
                         t.response_rx_window,
                         initiator,
                         TimeReferenceEvent::RxEnable,
                         true);
    if (initiator && !t.poll_rx_window.value.is_zero())
        cfg_rej("timeouts.poll_rx_window.ns",
                ConfigReason::FieldConflict,
                "the initiator never receives a Poll, so its Poll RX window must be zero",
                "REQ-API-02");
    if (!initiator && !t.response_rx_window.value.is_zero())
        cfg_rej(
            "timeouts.response_rx_window.ns",
            ConfigReason::FieldConflict,
            "the responder never receives a Response, so its Response RX window must be "
            "zero",
            "REQ-API-02");
    if (ss) {
        if (!t.final_rx_window.value.is_zero())
            cfg_rej("timeouts.final_rx_window.ns",
                    ConfigReason::FieldConflict,
                    "SS-TWR has no Final message, so the Final RX window must be zero",
                    "REQ-API-02");
    } else {
        validate_timed_field("timeouts.final_rx_window",
                             t.final_rx_window,
                             !initiator,
                             TimeReferenceEvent::RxEnable,
                             true);
    }
    if (!t.report_rx_window.value.is_zero())
        unsup("timeouts.report_rx_window",
              "FrameType::Report is reserved in phase 1 (REQ-PROTO-05)",
              ConfigReason::OutOfScope,
              "REQ-PROTO-05");

    // Every active window must fit inside the RX timeout measured from the
    // same RX enable event, otherwise the timeout can fire before the
    // window closes.
    const struct {
        const char* path;
        const TimedField& f;
    } windows[] = {
        { "timeouts.poll_rx_window", t.poll_rx_window },
        { "timeouts.response_rx_window", t.response_rx_window },
        { "timeouts.final_rx_window", t.final_rx_window },
    };
    for (const auto& w : windows) {
        if (w.f.value.is_zero())
            continue;
        if (w.f.value > t.rx_timeout.value)
            cfg_rej(std::string(w.path) + ".ns",
                    ConfigReason::TimingOrderViolation,
                    "the RX window " + twr_duration_to_text(w.f.value) +
                        " is longer than the RX timeout " +
                        twr_duration_to_text(t.rx_timeout.value) +
                        " measured from the same RX enable",
                    "REQ-API-02");
    }

    // The whole-exchange deadline must cover the chain this endpoint runs.
    int64_t chain_ns = t.rx_timeout.value.nanos();
    if (!initiator)
        chain_ns += c.timing.poll_to_response.value.nanos();
    if (!ss && initiator)
        chain_ns += c.timing.response_to_final.value.nanos();
    if (chain_ns > 0 && t.exchange_timeout.value.nanos() < chain_ns)
        cfg_rej("timeouts.exchange_timeout.ns",
                ConfigReason::TimingOrderViolation,
                "the whole-exchange deadline " +
                    twr_duration_to_text(t.exchange_timeout.value) +
                    " is shorter than the chain this endpoint must run (" +
                    twr_duration_to_text(Duration::from_nanos(chain_ns)) + ")",
                "REQ-API-02");

    // Consecutive measurements may not overlap: at most one exchange is in
    // flight, so the interval must be at least one whole exchange.
    if (c.session.measurement_count > 1 &&
        c.session.measurement_interval.nanos() < t.exchange_timeout.value.nanos())
        cfg_rej(
            "session.measurement_interval.ns",
            ConfigReason::TimingOrderViolation,
            "at most one exchange is in flight, so the measurement interval must be at "
            "least one whole-exchange timeout",
            "REQ-API-01");
}


void ConfigValidator::check_calibration(const TwrConfig& c) const
{
    const TimestampCalibrationConfig& k = c.calibration;
    if (k.tx_link_delay.negative() || k.rx_link_delay.negative() ||
        k.antenna_delay.negative() || k.cable_delay.negative())
        cfg_rej("calibration.tx_link_delay",
                ConfigReason::NegativeValue,
                "a link/antenna/cable delay may not be negative",
                "REQ-CAL-01");
    finite("calibration.native_sample_rate_hz", k.native_sample_rate_hz, "REQ-PHY-02");
    check_finite_opt("calibration.record.gain_db", k.record.gain_db);
    if (k.native_sample_rate_hz > 0.0 && !is_allowed_native_rate(k.native_sample_rate_hz))
        cfg_rej("calibration.native_sample_rate_hz",
                ConfigReason::OutOfRange,
                "calibration is versioned by sample rate; this build only calibrates "
                "at 737280000 or 491520000",
                "REQ-CAL-01");
    if (k.tx_link_delay_native_ticks < 0 || k.rx_link_delay_native_ticks < 0)
        cfg_rej("calibration.tx_link_delay_native_ticks",
                ConfigReason::NegativeValue,
                "a native-tick delay may not be negative",
                "REQ-CAL-01");
    if (k.link_delay_unit == TimeUnit::Seconds)
        cfg_rej("calibration.link_delay_unit",
                ConfigReason::UnitMismatch,
                "the delay fields are integer nanoseconds; use link_delay_unit=ns or "
                "native_ticks and state the seconds value elsewhere",
                "REQ-TIME-01");
    if (k.link_delay_unit == TimeUnit::NativeTicks && k.native_sample_rate_hz <= 0.0)
        cfg_rej("calibration.native_sample_rate_hz",
                ConfigReason::UnitMismatch,
                "native_ticks needs the native sample rate it refers to",
                "REQ-TIME-01");
    if (k.link_delay_unit == TimeUnit::NativeTicks &&
        (k.tx_link_delay_native_ticks == 0 && k.rx_link_delay_native_ticks == 0))
        cfg_rej("calibration.link_delay_unit",
                ConfigReason::FieldConflict,
                "link_delay_unit=native_ticks but no native-tick delay is stated",
                "REQ-TIME-01");

    if (k.applied_count != 0)
        rej("calibration.applied_count",
            ConfigReason::CalibrationAlreadyApplied,
            ExchangeStatus::CalibrationMissing,
            "a calibration is applicable exactly once; this one has already been "
            "applied " +
                twr_int_to_text(k.applied_count) + " time(s)",
            "REQ-CAL-01");

    if (!k.calibration_required)
        return;
    if (k.calibration_id.empty()) {
        rej("calibration.calibration_id",
            ConfigReason::CalibrationMissing,
            ExchangeStatus::CalibrationMissing,
            "a calibration id is required; an uncalibrated absolute range may not be "
            "claimed",
            "REQ-CAL-01");
        return;
    }
    if (k.record.calibration_id != k.calibration_id)
        cfg_rej("calibration.record.calibration_id",
                ConfigReason::FieldConflict,
                "the applicability record names a different calibration id",
                "REQ-CAL-01");
    if (k.record.channel != c.phy.channel)
        cfg_rej("calibration.record.channel",
                ConfigReason::CalibrationMismatch,
                "calibration channel " + twr_int_to_text(k.record.channel) +
                    " does not match phy.channel " + twr_int_to_text(c.phy.channel),
                "REQ-CAL-01");
    // M0.1: the applicability rule below is guarded by `> 0.0`, and a NaN
    // fails that comparison, so a NaN here used to sail through the whole
    // validator.  The record is compared against a rate, so it must be
    // finite whether or not the match ends up being attempted.
    // (record.gain_db is already covered by check_finite_opt above.)
    finite("calibration.record.native_sample_rate_hz",
           k.record.native_sample_rate_hz,
           "REQ-CAL-01");
    if (k.record.native_sample_rate_hz > 0.0 && c.radio.native_sample_rate_hz > 0.0 &&
        !twr_rate_matches_strict(
            k.record.native_sample_rate_hz, c.radio.native_sample_rate_hz, 1e-9))
        cfg_rej("calibration.record.native_sample_rate_hz",
                ConfigReason::CalibrationMismatch,
                "calibration sample rate does not match the radio native rate",
                "REQ-CAL-01");
    if (!k.record.profile_version.empty() && !c.meta.profile_version.empty() &&
        k.record.profile_version != c.meta.profile_version)
        cfg_rej("calibration.record.profile_version",
                ConfigReason::CalibrationMismatch,
                "calibration profile_version '" + k.record.profile_version +
                    "' does not match the session profile_version '" +
                    c.meta.profile_version + "'",
                "REQ-CAL-01");
    if (!c.meta.calibration_version.empty() &&
        k.calibration_id != c.meta.calibration_version)
        cfg_rej("meta.calibration_version",
                ConfigReason::CalibrationMismatch,
                "meta.calibration_version does not name the calibration that is loaded",
                "REQ-OUT-01");
    if (c.tx.gain_db.has() && k.record.gain_db.has() &&
        k.record.gain_db.value() != c.tx.gain_db.value())
        cfg_rej("calibration.record.gain_db",
                ConfigReason::CalibrationMismatch,
                "the calibration was taken at a different TX gain",
                "REQ-CAL-01");
}


void ConfigValidator::check_diagnostics(const TwrConfig& c) const
{
    const DiagnosticsConfig& d = c.diagnostics;
    const struct {
        const char* on;
        const char* bytes;
        const char* stride;
        bool enabled;
        uint64_t bytes_v;
        uint32_t stride_v;
    } buffers[] = {
        { "diagnostics.cir_capture_enabled",
          "diagnostics.cir_capture_max_bytes",
          "diagnostics.cir_capture_stride",
          d.cir_capture_enabled,
          d.cir_capture_max_bytes,
          d.cir_capture_stride },
        { "diagnostics.short_iq_enabled",
          "diagnostics.short_iq_max_bytes",
          "diagnostics.short_iq_stride",
          d.short_iq_enabled,
          d.short_iq_max_bytes,
          d.short_iq_stride },
        { "diagnostics.raw_frame_dump",
          "diagnostics.raw_frame_max_bytes",
          nullptr,
          d.raw_frame_dump,
          d.raw_frame_max_bytes,
          0u },
    };
    for (const auto& b : buffers) {
        if (b.bytes_v > kTwrMaxDiagnosticBytes)
            cfg_rej(b.bytes,
                    ConfigReason::OverCapacity,
                    "diagnostic buffer above the " +
                        twr_int_to_text(static_cast<int64_t>(kTwrMaxDiagnosticBytes)) +
                        " byte bound (bounded memory, REQ-GR-03)");
        if (b.enabled && b.bytes_v == 0)
            cfg_rej(b.bytes,
                    ConfigReason::ZeroValue,
                    "the diagnostic is enabled but no byte budget was stated");
        if (!b.enabled && b.bytes_v != 0)
            cfg_rej(b.bytes,
                    ConfigReason::FieldConflict,
                    "a byte budget is stated but the diagnostic is disabled");
        if (b.stride != nullptr && b.enabled && b.stride_v == 0)
            cfg_rej(b.stride,
                    ConfigReason::ZeroValue,
                    "a diagnostic stride of 0 would capture nothing");
        if (b.stride != nullptr && !b.enabled && b.stride_v != 0)
            cfg_rej(b.stride,
                    ConfigReason::FieldConflict,
                    "a diagnostic stride is stated but the diagnostic is disabled");
    }
    if (d.result_queue_capacity == 0)
        cfg_rej("diagnostics.result_queue_capacity",
                ConfigReason::ZeroValue,
                "the result queue must hold at least one entry");
    else if (d.result_queue_capacity > kTwrMaxQueueEntries)
        cfg_rej("diagnostics.result_queue_capacity",
                ConfigReason::OverCapacity,
                "above the bounded queue limit");
    if (d.event_queue_capacity == 0)
        cfg_rej("diagnostics.event_queue_capacity",
                ConfigReason::ZeroValue,
                "the event queue must hold at least one entry");
    else if (d.event_queue_capacity > kTwrMaxQueueEntries)
        cfg_rej("diagnostics.event_queue_capacity",
                ConfigReason::OverCapacity,
                "above the bounded queue limit");
    if (!d.result_output_path.empty() && d.result_output_path.size() > 4096)
        cfg_rej("diagnostics.result_output_path",
                ConfigReason::OverCapacity,
                "absurdly long output path");
    if (d.io_on_realtime_thread)
        cfg_rej("diagnostics.io_on_realtime_thread",
                ConfigReason::FieldConflict,
                "diagnostic I/O must not run on the realtime processing thread",
                "REQ-GR-02");
    check_finite_timed("diagnostics.stats_cadence", d.stats_cadence);
    validate_timed_field("diagnostics.stats_cadence",
                         d.stats_cadence,
                         false,
                         TimeReferenceEvent::HostMonotonic,
                         false);
}


void ConfigValidator::check_frame_lengths(const TwrConfig& c) const
{
    const FrameFormatConfig& f = c.frame;
    if (f.mac_psdu_bytes == 0)
        return; // already reported as a missing value elsewhere
    const uint16_t cap = d_caps.max_psdu_bytes();
    if (cap != 0 && f.mac_psdu_bytes > cap)
        rej("frame.mac_psdu_bytes",
            ConfigReason::FrameLengthOverflow,
            ExchangeStatus::Unsupported,
            "negotiated PSDU " + twr_int_to_text(f.mac_psdu_bytes) +
                " bytes exceeds the capability maximum " + twr_int_to_text(cap) +
                " bytes",
            "REQ-PHY-01");

    // Every length below comes from the codec, through
    // frame_psdu_bytes() -> FrameProfileGeometry::mac_payload_bytes() ->
    // frame_length_for().  The config's geometry is a claim, already
    // checked field by field in check_frame(), and is used here ONLY to
    // describe the arithmetic back to the operator.
    FrameProfileGeometry authority;
    FrameProfile profile;
    const bool have_profile = frame_profile_id_is_supported(f.frame_profile) &&
                              frame_geometry_for(f.frame_profile, authority) &&
                              frame_profile_for(f.frame_profile, profile);

    const FrameType types[] = { FrameType::Poll, FrameType::Response, FrameType::Final };
    for (FrameType t : types) {
        // `mac_psdu_bytes` is the negotiated MAC PSDU and, when the PHY
        // appends the FCS, EXCLUDES it (see the fcs_append /
        // mac_psdu_includes_fcs cross-check above).  So the comparison
        // below must be against the MAC bytes only.  The 127-byte IEEE
        // limit is separately checked against the on-air length.
        const uint32_t need = frame_psdu_bytes(t, f.fcs_append);
        const uint32_t on_air =
            have_profile
                ? static_cast<uint32_t>(authority.on_air_bytes(t, profile))
                : need + ((f.fcs_append == FcsAppender::PhyLayer) ? f.fcs_bytes : 0u);
        const std::string path =
            std::string("frame.mac_psdu_bytes[") + frame_type_to_string(t) + "]";
        if (need == 0)
            continue; // the profile itself was already reported as unknown
        if (need > static_cast<uint32_t>(f.mac_psdu_bytes))
            rej(path,
                ConfigReason::FrameLengthOverflow,
                ExchangeStatus::Unsupported,
                std::string("the ") + frame_type_to_string(t) + " frame needs " +
                    twr_int_to_text(static_cast<int64_t>(need)) +
                    " MAC PSDU bytes (as frame profile " +
                    frame_profile_id_to_string(f.frame_profile) + " builds it: header " +
                    twr_int_to_text(authority.mac_header_bytes()) + " + " +
                    twr_int_to_text(frame_type_timestamp_count(t)) + "x" +
                    twr_int_to_text(authority.timestamp_bytes(profile)) +
                    " timestamps + footer " +
                    twr_int_to_text(authority.mac_footer_bytes()) +
                    (f.fcs_append == FcsAppender::MacLayer
                         ? std::string(" + FCS ") + twr_int_to_text(f.fcs_bytes)
                         : std::string("; the PHY appends the FCS separately")) +
                    ") but only " + twr_int_to_text(f.mac_psdu_bytes) + " are negotiated",
                "REQ-PHY-01");
        // The IEEE 127-byte limit applies to the ON-AIR PSDU, i.e. with
        // the FCS the PHY appends (UwbHrpPacketSource appends and then
        // checks 127), so this is the check that uses `on_air`.
        if (cap != 0 && on_air > static_cast<uint32_t>(cap))
            rej(path,
                ConfigReason::FrameLengthOverflow,
                ExchangeStatus::Unsupported,
                std::string("the ") + frame_type_to_string(t) + " frame needs " +
                    twr_int_to_text(static_cast<int64_t>(on_air)) +
                    " on-air PSDU bytes (MAC " +
                    twr_int_to_text(static_cast<int64_t>(need)) + " + FCS " +
                    twr_int_to_text(authority.fcs_bytes_on_air()) +
                    ") and does not fit the " + twr_int_to_text(cap) +
                    "-byte PSDU maximum",
                "REQ-PHY-01");
    }
    // The application payload is a separate concept and must fit inside
    // the MAC PSDU after the header/timestamps/FCS are removed.
    if (f.geometry.any_set() && f.application_payload_bytes > f.mac_psdu_bytes)
        rej("frame.application_payload_bytes",
            ConfigReason::FrameLengthOverflow,
            ExchangeStatus::Unsupported,
            "the application payload (" + twr_int_to_text(f.application_payload_bytes) +
                " bytes) does not fit the MAC PSDU (" +
                twr_int_to_text(f.mac_psdu_bytes) + " bytes)",
            "REQ-API-01");

    // The joint PHY capability lookup, now that every axis is known.
    const CapabilityLookup look = d_caps.lookup_phy(c.radio.native_sample_rate_hz,
                                                    c.phy.tx_preamble_code,
                                                    c.phy.preamble_symbols,
                                                    f.sfd_mode,
                                                    f.mac_psdu_bytes,
                                                    f.ranging_bit);
    if (!look.allowed) {
        rej("phy.preamble_symbols",
            ConfigReason::Unsupported,
            ExchangeStatus::Unsupported,
            std::string(
                "the (native_rate, code, sync, sfd, psdu, ranging) combination is not "
                "in the capability whitelist [") +
                capability_status_to_string(look.status) + "]: " + look.reason,
            "REQ-PHY-01");
    }
}

bool to_json_string(const TwrConfig& c, std::string& out, std::string& err)
{
    ConfigJsonWriter w;
    write_config_json(c, w);
    if (!w.ok()) {
        err = w.error();
        return false;
    }
    return json::dump(w.value(), out, err, /*pretty=*/true);
}

void write_config_json(const TwrConfig& c, ConfigJsonWriter& w)
{
    w.current_object();
    w.begin_group("meta");
    w.put("schema_version", c.meta.schema_version);
    w.put("profile_version", c.meta.profile_version);
    w.put("calibration_version", c.meta.calibration_version);
    w.put("label", c.meta.label);
    w.end();

    w.begin_group("session");
    w.put("protocol", protocol_to_string(c.session.protocol));
    w.put("role", role_to_string(c.session.role));
    w.put("local_address", c.session.local_address);
    w.put("peer_address", c.session.peer_address);
    w.put("pan_id", c.session.pan_id);
    w.put("session_id", c.session.session_id);
    w.put("exchange_id", c.session.exchange_id);
    w.put("sequence", c.session.sequence);
    w.put("sequence_modulus", c.session.sequence_modulus);
    w.put("measurement_count", c.session.measurement_count);
    w.put_ns("measurement_interval_ns", c.session.measurement_interval);
    w.put("max_attempts_per_exchange", c.session.max_attempts_per_exchange);
    w.put_ns("retry_backoff_ns", c.session.retry_backoff);
    w.put("max_in_flight_exchanges", c.session.max_in_flight_exchanges);
    w.put("require_pan_match", c.session.require_pan_match);
    w.put("require_address_match", c.session.require_address_match);
    w.end();

    w.begin_group("phy");
    w.put("channel", c.phy.channel);
    w.put("center_frequency_hz", c.phy.center_frequency_hz);
    w.put("tx_preamble_code", c.phy.tx_preamble_code);
    w.put("rx_preamble_code", c.phy.rx_preamble_code);
    w.put("preamble_symbols", c.phy.preamble_symbols);
    w.put("prf_class", prf_class_to_string(c.phy.prf_class));
    w.put("data_rate", data_rate_to_string(c.phy.data_rate));
    // PHR rate, in the PHR-rate spelling ("850k" / "same_as_data"), NOT the
    // payload-rate spelling: the two axes have separate value domains and a
    // document that used "6p8m" here was describing a PHR this build cannot
    // transmit.
    w.put("phr_rate", phr_rate_to_string(c.phy.phr_rate));
    w.end();

    w.begin_group("frame");
    w.put("frame_profile", frame_profile_id_to_string(c.frame.frame_profile));
    w.put("sfd_mode", sfd_mode_to_string(c.frame.sfd_mode));
    w.put("sfd_symbols", c.frame.sfd_symbols);
    w.put_timed("sfd_timeout", c.frame.sfd_timeout);
    w.put("phr_mode", phr_mode_to_string(c.frame.phr_mode));
    w.put("ranging_bit", c.frame.ranging_bit);
    // The geometry CLAIM is written out verbatim, so a reader can see what the
    // operator asserted; the validator is what decides whether it agrees with
    // the codec.
    w.begin_group("geometry");
    w.put("mac_header_bytes", c.frame.geometry.mac_header_bytes);
    w.put("timestamp_bytes", c.frame.geometry.timestamp_bytes);
    w.put("mac_footer_bytes", c.frame.geometry.mac_footer_bytes);
    w.put("mac_fcs_bytes", c.frame.geometry.mac_fcs_bytes);
    w.put("phr_bytes", c.frame.geometry.phr_bytes);
    w.end();
    w.put("mac_psdu_bytes", c.frame.mac_psdu_bytes);
    w.put("mac_psdu_includes_fcs", c.frame.mac_psdu_includes_fcs);
    w.put("fcs_append", fcs_appender_to_string(c.frame.fcs_append));
    w.put("fcs_bytes", c.frame.fcs_bytes);
    w.put("application_payload_bytes", c.frame.application_payload_bytes);
    w.put("sts_mode", sts_mode_to_string(c.frame.sts_mode));
    w.put("sts_length_symbols", c.frame.sts_length_symbols);
    w.end();

    w.begin_group("tx");
    w.put("port", c.tx.port);
    w.put_opt("gain_db", c.tx.gain_db);
    w.put_opt("iq_amplitude", c.tx.iq_amplitude);
    w.put_opt("calibrated_tx_power_dbm", c.tx.calibrated_tx_power_dbm);
    w.put("power_policy", tx_power_policy_to_string(c.tx.power_policy));
    w.put("pulse_shaping", pulse_shaping_to_string(c.tx.pulse_shaping));
    w.put_opt("vendor_power_word", c.tx.vendor_power_word);
    w.put("vendor_power_word_backend", c.tx.vendor_power_word_backend);
    w.end();

    w.begin_group("rx");
    w.put("port", c.rx.port);
    w.put_opt("gain_db", c.rx.gain_db);
    w.put("agc", agc_mode_to_string(c.rx.agc));
    w.put_opt("bandwidth_hz", c.rx.bandwidth_hz);
    w.put("detection_threshold", c.rx.detection_threshold);
    w.put("correlation_threshold", c.rx.correlation_threshold);
    w.put("first_path_threshold", c.rx.first_path_threshold);
    w.put("first_path_index", c.rx.first_path_index);
    w.put("first_path_window", c.rx.first_path_window);
    w.put_opt("vendor_pac_value", c.rx.vendor_pac_value);
    w.put("vendor_pac_backend", c.rx.vendor_pac_backend);
    w.put_opt("vendor_pac_applied_step", c.rx.vendor_pac_applied_step);
    w.end();

    w.begin_group("radio");
    w.put("device_args", c.radio.device_args);
    w.put("tx_channel", c.radio.tx_channel);
    w.put("rx_channel", c.radio.rx_channel);
    w.put("native_sample_rate_hz", c.radio.native_sample_rate_hz);
    w.put("clock_source", c.radio.clock_source);
    w.put("time_source", c.radio.time_source);
    w.put("fpga_image", c.radio.fpga_image);
    w.put("dpdk_config", c.radio.dpdk_config);
    w.begin_array("peers");
    for (const auto& p : c.radio.peers) {
        w.push_object();
        w.put("id", p.id);
        w.put("role", role_to_string(p.role));
        w.put("tx_channel", p.tx_channel);
        w.put("rx_channel", p.rx_channel);
        w.put("native_sample_rate_hz", p.native_sample_rate_hz);
        w.put("occupies_resources", p.occupies_resources);
        w.pop_object();
    }
    w.end();
    w.begin_group("readback");
    w.put("present", c.radio.readback.present);
    w.put("sample_rate_hz", c.radio.readback.sample_rate_hz);
    w.put("center_freq_hz", c.radio.readback.center_freq_hz);
    w.put("tx_channel", c.radio.readback.tx_channel);
    w.put("rx_channel", c.radio.readback.rx_channel);
    w.put("mpm_string", c.radio.readback.mpm_string);
    w.put("fpga_image", c.radio.readback.fpga_image);
    w.put("uhd_version", c.radio.readback.uhd_version);
    w.put("clock_source", c.radio.readback.clock_source);
    w.put("time_source", c.radio.readback.time_source);
    w.end();
    w.put("require_readback", c.radio.require_readback);
    w.end();

    w.begin_group("timing");
    w.put_timed("poll_start", c.timing.poll_start);
    w.put_timed("poll_to_response", c.timing.poll_to_response);
    w.put_timed("response_to_final", c.timing.response_to_final);
    w.put_timed("final_to_report", c.timing.final_to_report);
    w.put_timed("post_tx_rx_enable", c.timing.post_tx_rx_enable);
    w.put_timed("min_tx_lead_time", c.timing.min_tx_lead_time);
    w.end();

    w.begin_group("timeouts");
    w.put_timed("poll_rx_window", c.timeouts.poll_rx_window);
    w.put_timed("response_rx_window", c.timeouts.response_rx_window);
    w.put_timed("final_rx_window", c.timeouts.final_rx_window);
    w.put_timed("report_rx_window", c.timeouts.report_rx_window);
    w.put_timed("rx_timeout", c.timeouts.rx_timeout);
    w.put_timed("exchange_timeout", c.timeouts.exchange_timeout);
    w.put_timed("retry_interval", c.timeouts.retry_interval);
    w.end();

    w.begin_group("calibration");
    w.put_ns("tx_link_delay_ns", c.calibration.tx_link_delay);
    w.put_ns("rx_link_delay_ns", c.calibration.rx_link_delay);
    w.put_ns("antenna_delay_ns", c.calibration.antenna_delay);
    w.put_ns("cable_delay_ns", c.calibration.cable_delay);
    w.put("tx_link_delay_native_ticks", c.calibration.tx_link_delay_native_ticks);
    w.put("rx_link_delay_native_ticks", c.calibration.rx_link_delay_native_ticks);
    w.put("native_sample_rate_hz", c.calibration.native_sample_rate_hz);
    w.put("link_delay_unit", time_unit_to_string(c.calibration.link_delay_unit));
    w.put("first_path_algorithm",
          first_path_algorithm_to_string(c.calibration.first_path_algorithm));
    w.put("cfo_compensation",
          compensation_flag_to_string(c.calibration.cfo_compensation));
    w.put("sfo_compensation",
          compensation_flag_to_string(c.calibration.sfo_compensation));
    w.put("calibration_id", c.calibration.calibration_id);
    w.begin_group("record");
    w.put("calibration_id", c.calibration.record.calibration_id);
    w.put("device_serial", c.calibration.record.device_serial);
    w.put("channel", c.calibration.record.channel);
    w.put("native_sample_rate_hz", c.calibration.record.native_sample_rate_hz);
    w.put("profile_version", c.calibration.record.profile_version);
    w.put_opt("gain_db", c.calibration.record.gain_db);
    w.put("valid_until_monotonic_ns", c.calibration.record.valid_until_monotonic_ns);
    w.end();
    w.put("applied_count", c.calibration.applied_count);
    w.put("calibration_required", c.calibration.calibration_required);
    w.end();

    w.begin_group("diagnostics");
    w.put("cir_capture_enabled", c.diagnostics.cir_capture_enabled);
    w.put("cir_capture_max_bytes", c.diagnostics.cir_capture_max_bytes);
    w.put("cir_capture_stride", c.diagnostics.cir_capture_stride);
    w.put("short_iq_enabled", c.diagnostics.short_iq_enabled);
    w.put("short_iq_max_bytes", c.diagnostics.short_iq_max_bytes);
    w.put("short_iq_stride", c.diagnostics.short_iq_stride);
    w.put("raw_frame_dump", c.diagnostics.raw_frame_dump);
    w.put("raw_frame_max_bytes", c.diagnostics.raw_frame_max_bytes);
    w.put("result_output_path", c.diagnostics.result_output_path);
    w.put("result_queue_capacity", c.diagnostics.result_queue_capacity);
    w.put("event_queue_capacity", c.diagnostics.event_queue_capacity);
    w.put_timed("stats_cadence", c.diagnostics.stats_cadence);
    w.put("io_on_realtime_thread", c.diagnostics.io_on_realtime_thread);
    w.end();
}


ValidationReport from_json_string(const std::string& text, TwrConfig& out)
{
    ConfigJsonReader rd;
    json::Value root;
    std::string err;
    if (!json::parse(text, root, err)) {
        rd.bad("", ConfigReason::MalformedJson, err);
        return rd.report();
    }
    if (root.type != json::Type::Object) {
        rd.bad(
            "", ConfigReason::TypeMismatch, "the config document must be a JSON object");
        return rd.report();
    }
    TwrConfig c;

    {
        const std::string p = "meta";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "schema_version",
                                                 "profile_version",
                                                 "calibration_version",
                                                 "label",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.meta.schema_version = rd.get_str(g, p, "schema_version");
            c.meta.profile_version = rd.get_str(g, p, "profile_version");
            c.meta.calibration_version = rd.get_str(g, p, "calibration_version");
            c.meta.label = rd.get_str(g, p, "label");
        }
    }
    {
        const std::string p = "session";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "protocol",
                                                 "role",
                                                 "local_address",
                                                 "peer_address",
                                                 "pan_id",
                                                 "session_id",
                                                 "exchange_id",
                                                 "sequence",
                                                 "sequence_modulus",
                                                 "measurement_count",
                                                 "measurement_interval_ns",
                                                 "max_attempts_per_exchange",
                                                 "retry_backoff_ns",
                                                 "max_in_flight_exchanges",
                                                 "require_pan_match",
                                                 "require_address_match",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.session.protocol =
                rd.get_enum<Protocol>(g, p, "protocol", protocol_from_string);
            c.session.role = rd.get_enum<Role>(g, p, "role", role_from_string);
            c.session.local_address =
                static_cast<uint16_t>(rd.get_u16(g, p, "local_address"));
            c.session.peer_address =
                static_cast<uint16_t>(rd.get_u16(g, p, "peer_address"));
            c.session.pan_id = static_cast<uint16_t>(rd.get_u16(g, p, "pan_id"));
            c.session.session_id = static_cast<uint32_t>(rd.get_u32(g, p, "session_id"));
            c.session.exchange_id =
                static_cast<uint32_t>(rd.get_u32(g, p, "exchange_id"));
            c.session.sequence = static_cast<uint8_t>(rd.get_u8(g, p, "sequence"));
            c.session.sequence_modulus =
                static_cast<uint16_t>(rd.get_u16(g, p, "sequence_modulus"));
            c.session.measurement_count =
                static_cast<uint32_t>(rd.get_u32(g, p, "measurement_count"));
            c.session.measurement_interval = rd.get_dur(g, p, "measurement_interval_ns");
            c.session.max_attempts_per_exchange =
                static_cast<uint32_t>(rd.get_u32(g, p, "max_attempts_per_exchange"));
            c.session.retry_backoff = rd.get_dur(g, p, "retry_backoff_ns");
            c.session.max_in_flight_exchanges =
                static_cast<uint32_t>(rd.get_u32(g, p, "max_in_flight_exchanges"));
            c.session.require_pan_match = rd.get_bool(g, p, "require_pan_match");
            c.session.require_address_match = rd.get_bool(g, p, "require_address_match");
        }
    }
    {
        const std::string p = "phy";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "channel",
                                                 "center_frequency_hz",
                                                 "tx_preamble_code",
                                                 "rx_preamble_code",
                                                 "preamble_symbols",
                                                 "prf_class",
                                                 "data_rate",
                                                 "phr_rate",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.phy.channel = static_cast<uint8_t>(rd.get_u8(g, p, "channel"));
            c.phy.center_frequency_hz = rd.get_double(g, p, "center_frequency_hz");
            c.phy.tx_preamble_code =
                static_cast<uint8_t>(rd.get_u8(g, p, "tx_preamble_code"));
            c.phy.rx_preamble_code =
                static_cast<uint8_t>(rd.get_u8(g, p, "rx_preamble_code"));
            c.phy.preamble_symbols =
                static_cast<uint16_t>(rd.get_u16(g, p, "preamble_symbols"));
            c.phy.prf_class =
                rd.get_enum<PrfClass>(g, p, "prf_class", prf_class_from_string);
            c.phy.data_rate =
                rd.get_enum<DataRate>(g, p, "data_rate", data_rate_from_string);
            // PHR rate, read in the PHR-rate domain.  A payload-rate name such
            // as "6p8m" is now an UnknownEnumValue here, which is the point:
            // the two axes no longer share a value domain.
            c.phy.phr_rate = rd.get_enum<PhrRate>(g, p, "phr_rate", phr_rate_from_string);
        }
    }
    {
        const std::string p = "frame";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = {
                "frame_profile", "sfd_mode",           "sfd_symbols",
                "sfd_timeout",   "phr_mode",           "ranging_bit",
                "geometry",      "mac_psdu_bytes",     "mac_psdu_includes_fcs",
                "fcs_append",    "fcs_bytes",          "application_payload_bytes",
                "sts_mode",      "sts_length_symbols", nullptr
            };
            rd.reject_unknown_keys(g, p, kKeys);
            c.frame.frame_profile = rd.get_enum<FrameProfileId>(
                g, p, "frame_profile", frame_profile_id_from_string);
            c.frame.sfd_mode =
                rd.get_enum<SfdMode>(g, p, "sfd_mode", sfd_mode_from_string);
            c.frame.sfd_symbols = static_cast<uint16_t>(rd.get_u16(g, p, "sfd_symbols"));
            c.frame.sfd_timeout = rd.get_timed(g, p, "sfd_timeout");
            c.frame.phr_mode =
                rd.get_enum<PhrMode>(g, p, "phr_mode", phr_mode_from_string);
            c.frame.ranging_bit = rd.get_bool(g, p, "ranging_bit");
            const json::Value* geo = g->find("geometry");
            if (!geo)
                rd.bad(
                    p + ".geometry", ConfigReason::MissingKey, "missing required group");
            else if (geo->type != json::Type::Object)
                rd.bad(p + ".geometry", ConfigReason::TypeMismatch, "expected an object");
            else {
                static const char* const kGeo[] = { "mac_header_bytes", "timestamp_bytes",
                                                    "mac_footer_bytes", "mac_fcs_bytes",
                                                    "phr_bytes",        nullptr };
                rd.reject_unknown_keys(geo, p + ".geometry", kGeo);
                c.frame.geometry.mac_header_bytes = static_cast<uint16_t>(
                    rd.get_u16(geo, p + ".geometry", "mac_header_bytes"));
                c.frame.geometry.timestamp_bytes = static_cast<uint16_t>(
                    rd.get_u16(geo, p + ".geometry", "timestamp_bytes"));
                c.frame.geometry.mac_footer_bytes = static_cast<uint16_t>(
                    rd.get_u16(geo, p + ".geometry", "mac_footer_bytes"));
                c.frame.geometry.mac_fcs_bytes = static_cast<uint16_t>(
                    rd.get_u16(geo, p + ".geometry", "mac_fcs_bytes"));
                c.frame.geometry.phr_bytes =
                    static_cast<uint16_t>(rd.get_u16(geo, p + ".geometry", "phr_bytes"));
            }
            c.frame.mac_psdu_bytes =
                static_cast<uint16_t>(rd.get_u16(g, p, "mac_psdu_bytes"));
            c.frame.mac_psdu_includes_fcs = rd.get_bool(g, p, "mac_psdu_includes_fcs");
            c.frame.fcs_append =
                rd.get_enum<FcsAppender>(g, p, "fcs_append", fcs_appender_from_string);
            c.frame.fcs_bytes = static_cast<uint16_t>(rd.get_u16(g, p, "fcs_bytes"));
            c.frame.application_payload_bytes =
                static_cast<uint16_t>(rd.get_u16(g, p, "application_payload_bytes"));
            c.frame.sts_mode =
                rd.get_enum<StsMode>(g, p, "sts_mode", sts_mode_from_string);
            c.frame.sts_length_symbols =
                static_cast<uint16_t>(rd.get_u16(g, p, "sts_length_symbols"));
        }
    }
    {
        const std::string p = "tx";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "port",
                                                 "gain_db",
                                                 "iq_amplitude",
                                                 "calibrated_tx_power_dbm",
                                                 "power_policy",
                                                 "pulse_shaping",
                                                 "vendor_power_word",
                                                 "vendor_power_word_backend",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.tx.port = static_cast<uint8_t>(rd.get_u8(g, p, "port"));
            c.tx.gain_db = rd.get_opt_double(g, p, "gain_db");
            c.tx.iq_amplitude = rd.get_opt_double(g, p, "iq_amplitude");
            c.tx.calibrated_tx_power_dbm =
                rd.get_opt_double(g, p, "calibrated_tx_power_dbm");
            c.tx.power_policy = rd.get_enum<TxPowerPolicy>(
                g, p, "power_policy", tx_power_policy_from_string);
            c.tx.pulse_shaping = rd.get_enum<PulseShaping>(
                g, p, "pulse_shaping", pulse_shaping_from_string);
            c.tx.vendor_power_word = rd.get_opt_u32(g, p, "vendor_power_word");
            c.tx.vendor_power_word_backend =
                rd.get_str(g, p, "vendor_power_word_backend");
        }
    }
    {
        const std::string p = "rx";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "port",
                                                 "gain_db",
                                                 "agc",
                                                 "bandwidth_hz",
                                                 "detection_threshold",
                                                 "correlation_threshold",
                                                 "first_path_threshold",
                                                 "first_path_index",
                                                 "first_path_window",
                                                 "vendor_pac_value",
                                                 "vendor_pac_backend",
                                                 "vendor_pac_applied_step",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.rx.port = static_cast<uint8_t>(rd.get_u8(g, p, "port"));
            c.rx.gain_db = rd.get_opt_double(g, p, "gain_db");
            c.rx.agc = rd.get_enum<AgcMode>(g, p, "agc", agc_mode_from_string);
            c.rx.bandwidth_hz = rd.get_opt_double(g, p, "bandwidth_hz");
            c.rx.detection_threshold = rd.get_double(g, p, "detection_threshold");
            c.rx.correlation_threshold = rd.get_double(g, p, "correlation_threshold");
            c.rx.first_path_threshold = rd.get_double(g, p, "first_path_threshold");
            c.rx.first_path_index =
                static_cast<uint16_t>(rd.get_u16(g, p, "first_path_index"));
            c.rx.first_path_window =
                static_cast<uint16_t>(rd.get_u16(g, p, "first_path_window"));
            c.rx.vendor_pac_value = rd.get_opt_u32(g, p, "vendor_pac_value");
            c.rx.vendor_pac_backend = rd.get_str(g, p, "vendor_pac_backend");
            c.rx.vendor_pac_applied_step =
                rd.get_opt_double(g, p, "vendor_pac_applied_step");
        }
    }
    {
        const std::string p = "radio";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "device_args",
                                                 "tx_channel",
                                                 "rx_channel",
                                                 "native_sample_rate_hz",
                                                 "clock_source",
                                                 "time_source",
                                                 "fpga_image",
                                                 "dpdk_config",
                                                 "peers",
                                                 "readback",
                                                 "require_readback",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.radio.device_args = rd.get_str(g, p, "device_args");
            c.radio.tx_channel = static_cast<uint8_t>(rd.get_u8(g, p, "tx_channel"));
            c.radio.rx_channel = static_cast<uint8_t>(rd.get_u8(g, p, "rx_channel"));
            c.radio.native_sample_rate_hz = rd.get_double(g, p, "native_sample_rate_hz");
            c.radio.clock_source = rd.get_str(g, p, "clock_source");
            c.radio.time_source = rd.get_str(g, p, "time_source");
            c.radio.fpga_image = rd.get_str(g, p, "fpga_image");
            c.radio.dpdk_config = rd.get_str(g, p, "dpdk_config");
            const json::Value* peers = g->find("peers");
            if (!peers)
                rd.bad(p + ".peers", ConfigReason::MissingKey, "missing required key");
            else if (peers->type != json::Type::Array) {
                rd.bad(p + ".peers", ConfigReason::TypeMismatch, "expected an array");
            } else {
                for (size_t i = 0; i < peers->items.size(); ++i) {
                    const std::string pp =
                        p + ".peers[" + twr_int_to_text(static_cast<int64_t>(i)) + "]";
                    const json::Value& e = peers->items[i];
                    if (e.type != json::Type::Object) {
                        rd.bad(pp, ConfigReason::TypeMismatch, "expected an object");
                        continue;
                    }
                    static const char* const kPeer[] = { "id",
                                                         "role",
                                                         "tx_channel",
                                                         "rx_channel",
                                                         "native_sample_rate_hz",
                                                         "occupies_resources",
                                                         nullptr };
                    rd.reject_unknown_keys(&e, pp, kPeer);
                    EndpointBinding b;
                    b.id = rd.get_str(&e, pp, "id");
                    b.role = rd.get_enum<Role>(&e, pp, "role", role_from_string);
                    b.tx_channel = static_cast<uint8_t>(rd.get_u8(&e, pp, "tx_channel"));
                    b.rx_channel = static_cast<uint8_t>(rd.get_u8(&e, pp, "rx_channel"));
                    b.native_sample_rate_hz =
                        rd.get_double(&e, pp, "native_sample_rate_hz");
                    b.occupies_resources = rd.get_bool(&e, pp, "occupies_resources");
                    c.radio.peers.push_back(b);
                }
            }
            const json::Value* rb = g->find("readback");
            if (!rb)
                rd.bad(
                    p + ".readback", ConfigReason::MissingKey, "missing required group");
            else if (rb->type != json::Type::Object)
                rd.bad(p + ".readback", ConfigReason::TypeMismatch, "expected an object");
            else {
                const std::string rp = p + ".readback";
                static const char* const kRb[] = {
                    "present",      "sample_rate_hz", "center_freq_hz", "tx_channel",
                    "rx_channel",   "mpm_string",     "fpga_image",     "uhd_version",
                    "clock_source", "time_source",    nullptr
                };
                rd.reject_unknown_keys(rb, rp, kRb);
                c.radio.readback.present = rd.get_bool(rb, rp, "present");
                c.radio.readback.sample_rate_hz = rd.get_double(rb, rp, "sample_rate_hz");
                c.radio.readback.center_freq_hz = rd.get_double(rb, rp, "center_freq_hz");
                c.radio.readback.tx_channel =
                    static_cast<uint8_t>(rd.get_u8(rb, rp, "tx_channel"));
                c.radio.readback.rx_channel =
                    static_cast<uint8_t>(rd.get_u8(rb, rp, "rx_channel"));
                c.radio.readback.mpm_string = rd.get_str(rb, rp, "mpm_string");
                c.radio.readback.fpga_image = rd.get_str(rb, rp, "fpga_image");
                c.radio.readback.uhd_version = rd.get_str(rb, rp, "uhd_version");
                c.radio.readback.clock_source = rd.get_str(rb, rp, "clock_source");
                c.radio.readback.time_source = rd.get_str(rb, rp, "time_source");
            }
            c.radio.require_readback = rd.get_bool(g, p, "require_readback");
        }
    }
    {
        const std::string p = "timing";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "poll_start",
                                                 "poll_to_response",
                                                 "response_to_final",
                                                 "final_to_report",
                                                 "post_tx_rx_enable",
                                                 "min_tx_lead_time",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.timing.poll_start = rd.get_timed(g, p, "poll_start");
            c.timing.poll_to_response = rd.get_timed(g, p, "poll_to_response");
            c.timing.response_to_final = rd.get_timed(g, p, "response_to_final");
            c.timing.final_to_report = rd.get_timed(g, p, "final_to_report");
            c.timing.post_tx_rx_enable = rd.get_timed(g, p, "post_tx_rx_enable");
            c.timing.min_tx_lead_time = rd.get_timed(g, p, "min_tx_lead_time");
        }
    }
    {
        const std::string p = "timeouts";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "poll_rx_window",  "response_rx_window",
                                                 "final_rx_window", "report_rx_window",
                                                 "rx_timeout",      "exchange_timeout",
                                                 "retry_interval",  nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.timeouts.poll_rx_window = rd.get_timed(g, p, "poll_rx_window");
            c.timeouts.response_rx_window = rd.get_timed(g, p, "response_rx_window");
            c.timeouts.final_rx_window = rd.get_timed(g, p, "final_rx_window");
            c.timeouts.report_rx_window = rd.get_timed(g, p, "report_rx_window");
            c.timeouts.rx_timeout = rd.get_timed(g, p, "rx_timeout");
            c.timeouts.exchange_timeout = rd.get_timed(g, p, "exchange_timeout");
            c.timeouts.retry_interval = rd.get_timed(g, p, "retry_interval");
        }
    }
    {
        const std::string p = "calibration";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = { "tx_link_delay_ns",
                                                 "rx_link_delay_ns",
                                                 "antenna_delay_ns",
                                                 "cable_delay_ns",
                                                 "tx_link_delay_native_ticks",
                                                 "rx_link_delay_native_ticks",
                                                 "native_sample_rate_hz",
                                                 "link_delay_unit",
                                                 "first_path_algorithm",
                                                 "cfo_compensation",
                                                 "sfo_compensation",
                                                 "calibration_id",
                                                 "record",
                                                 "applied_count",
                                                 "calibration_required",
                                                 nullptr };
            rd.reject_unknown_keys(g, p, kKeys);
            c.calibration.tx_link_delay = rd.get_dur(g, p, "tx_link_delay_ns");
            c.calibration.rx_link_delay = rd.get_dur(g, p, "rx_link_delay_ns");
            c.calibration.antenna_delay = rd.get_dur(g, p, "antenna_delay_ns");
            c.calibration.cable_delay = rd.get_dur(g, p, "cable_delay_ns");
            c.calibration.tx_link_delay_native_ticks =
                rd.get_i64(g, p, "tx_link_delay_native_ticks");
            c.calibration.rx_link_delay_native_ticks =
                rd.get_i64(g, p, "rx_link_delay_native_ticks");
            c.calibration.native_sample_rate_hz =
                rd.get_double(g, p, "native_sample_rate_hz");
            c.calibration.link_delay_unit =
                rd.get_enum<TimeUnit>(g, p, "link_delay_unit", time_unit_from_string);
            c.calibration.first_path_algorithm = rd.get_enum<FirstPathAlgorithm>(
                g, p, "first_path_algorithm", first_path_algorithm_from_string);
            c.calibration.cfo_compensation = rd.get_enum<CompensationFlag>(
                g, p, "cfo_compensation", compensation_flag_from_string);
            c.calibration.sfo_compensation = rd.get_enum<CompensationFlag>(
                g, p, "sfo_compensation", compensation_flag_from_string);
            c.calibration.calibration_id = rd.get_str(g, p, "calibration_id");
            const json::Value* rec = g->find("record");
            if (!rec)
                rd.bad(p + ".record", ConfigReason::MissingKey, "missing required group");
            else if (rec->type != json::Type::Object)
                rd.bad(p + ".record", ConfigReason::TypeMismatch, "expected an object");
            else {
                const std::string rp = p + ".record";
                static const char* const kRec[] = { "calibration_id",
                                                    "device_serial",
                                                    "channel",
                                                    "native_sample_rate_hz",
                                                    "profile_version",
                                                    "gain_db",
                                                    "valid_until_monotonic_ns",
                                                    nullptr };
                rd.reject_unknown_keys(rec, rp, kRec);
                c.calibration.record.calibration_id =
                    rd.get_str(rec, rp, "calibration_id");
                c.calibration.record.device_serial = rd.get_str(rec, rp, "device_serial");
                c.calibration.record.channel =
                    static_cast<uint8_t>(rd.get_u8(rec, rp, "channel"));
                c.calibration.record.native_sample_rate_hz =
                    rd.get_double(rec, rp, "native_sample_rate_hz");
                c.calibration.record.profile_version =
                    rd.get_str(rec, rp, "profile_version");
                c.calibration.record.gain_db = rd.get_opt_double(rec, rp, "gain_db");
                c.calibration.record.valid_until_monotonic_ns =
                    rd.get_i64(rec, rp, "valid_until_monotonic_ns");
            }
            c.calibration.applied_count =
                static_cast<uint32_t>(rd.get_u32(g, p, "applied_count"));
            c.calibration.calibration_required =
                rd.get_bool(g, p, "calibration_required");
        }
    }
    {
        const std::string p = "diagnostics";
        const json::Value* g = rd.group(p, &root);
        if (g) {
            static const char* const kKeys[] = {
                "cir_capture_enabled",   "cir_capture_max_bytes",
                "cir_capture_stride",    "short_iq_enabled",
                "short_iq_max_bytes",    "short_iq_stride",
                "raw_frame_dump",        "raw_frame_max_bytes",
                "result_output_path",    "result_queue_capacity",
                "event_queue_capacity",  "stats_cadence",
                "io_on_realtime_thread", nullptr
            };
            rd.reject_unknown_keys(g, p, kKeys);
            c.diagnostics.cir_capture_enabled = rd.get_bool(g, p, "cir_capture_enabled");
            c.diagnostics.cir_capture_max_bytes = static_cast<uint64_t>(
                rd.get_u64(g, p, "cir_capture_max_bytes", 0xFFFFFFFFFFFFFFFFull));
            c.diagnostics.cir_capture_stride =
                static_cast<uint32_t>(rd.get_u32(g, p, "cir_capture_stride"));
            c.diagnostics.short_iq_enabled = rd.get_bool(g, p, "short_iq_enabled");
            c.diagnostics.short_iq_max_bytes = static_cast<uint64_t>(
                rd.get_u64(g, p, "short_iq_max_bytes", 0xFFFFFFFFFFFFFFFFull));
            c.diagnostics.short_iq_stride =
                static_cast<uint32_t>(rd.get_u32(g, p, "short_iq_stride"));
            c.diagnostics.raw_frame_dump = rd.get_bool(g, p, "raw_frame_dump");
            c.diagnostics.raw_frame_max_bytes = static_cast<uint64_t>(
                rd.get_u64(g, p, "raw_frame_max_bytes", 0xFFFFFFFFFFFFFFFFull));
            c.diagnostics.result_output_path = rd.get_str(g, p, "result_output_path");
            c.diagnostics.result_queue_capacity =
                static_cast<uint32_t>(rd.get_u32(g, p, "result_queue_capacity"));
            c.diagnostics.event_queue_capacity =
                static_cast<uint32_t>(rd.get_u32(g, p, "event_queue_capacity"));
            c.diagnostics.stats_cadence = rd.get_timed(g, p, "stats_cadence");
            c.diagnostics.io_on_realtime_thread =
                rd.get_bool(g, p, "io_on_realtime_thread");
        }
    }

    // Even a document that failed to import is handed back so the caller can
    // see which groups were readable; the report is the authority.
    out = c;
    return rd.report();
}

bool to_json_string(const EffectiveConfig& e, std::string& out, std::string& err)
{
    ConfigJsonWriter w;
    w.current_object();
    w.put("ok", e.ok);
    w.put("schema_version", e.schema_version);
    w.put("profile_version", e.profile_version);
    w.put("calibration_version", e.calibration_version);
    w.put("config_hash", e.config_hash);
    write_config_json(e.effective, w);
    w.begin_array("changes");
    for (const auto& ch : e.changes) {
        w.push_object();
        w.put("path", ch.path);
        w.put("requested", ch.requested);
        w.put("effective", ch.effective);
        w.put("note", ch.note);
        w.pop_object();
    }
    w.end();
    w.begin_array("violations");
    for (const auto& v : e.validation.violations) {
        w.push_object();
        w.put("field", v.field);
        w.put("reason", config_reason_to_string(v.reason));
        w.put("status", exchange_status_to_string(v.status));
        w.put("message", v.message);
        w.put("requirement", v.requirement);
        w.pop_object();
    }
    w.end();
    w.begin_group("quantised");
    w.put("tick_rate_hz", e.tick_rate_hz);
    w.put("poll_start_ticks", e.poll_start_ticks);
    w.put("poll_to_response_ticks", e.poll_to_response_ticks);
    w.put("response_to_final_ticks", e.response_to_final_ticks);
    w.put("post_tx_rx_enable_ticks", e.post_tx_rx_enable_ticks);
    w.put_ns("poll_to_response_effective_ns", e.poll_to_response_effective);
    w.put_ns("response_to_final_effective_ns", e.response_to_final_effective);
    w.put_ns("post_tx_rx_enable_effective_ns", e.post_tx_rx_enable_effective);
    w.end();
    w.begin_group("frame_budget");
    w.put("poll_bytes", e.poll_bytes);
    w.put("response_bytes", e.response_bytes);
    w.put("final_bytes", e.final_bytes);
    w.put("max_psdu_bytes", e.max_psdu_bytes);
    w.end();
    if (!w.ok()) {
        err = w.error();
        return false;
    }
    return json::dump(w.value(), out, err, /*pretty=*/true);
}


bool timestamp_normalize(int64_t ticks,
                         int32_t frac_num,
                         uint32_t frac_den,
                         int64_t& out_ticks,
                         int32_t& out_num,
                         uint32_t& out_den)
{
    if (frac_den == 0u) {
        if (frac_num != 0)
            return false; // an exact tick cannot carry a fraction
        return detail::ts_normalize_rational(ticks, 0, 0, out_ticks, out_num, out_den);
    }
    if (frac_den > detail::kFractionDenominatorMax)
        return false;
    return detail::ts_normalize_rational(ticks,
                                         static_cast<int64_t>(frac_num),
                                         static_cast<int64_t>(frac_den),
                                         out_ticks,
                                         out_num,
                                         out_den);
}


bool timestamp_make(int64_t ticks,
                    int32_t frac_num,
                    uint32_t frac_den,
                    const ClockDomain& domain,
                    TimestampMarker marker,
                    TimestampSource source,
                    uint32_t corrections,
                    Timestamp& out)
{
    if (!clock_domain_is_valid(domain))
        return false;
    if (!timestamp_marker_is_known(marker))
        return false;
    if (!timestamp_source_is_known(source) || source == TimestampSource::Unknown)
        return false;
    if (corrections > kCorrectionAll)
        return false;

    int64_t nt = 0;
    int32_t nn = 0;
    uint32_t nd = 0u;
    if (!timestamp_normalize(ticks, frac_num, frac_den, nt, nn, nd))
        return false;
    if (!clock_domain_ticks_in_range(domain, nt))
        return false;

    Timestamp ts;
    ts.ticks = nt;
    ts.frac_num = nn;
    ts.frac_den = nd;
    ts.domain = domain;
    ts.marker = marker;
    ts.source = source;
    ts.applied_corrections = corrections;
    ts.calibration_id.clear();
    ts.valid = true;
    out = ts;
    return true;
}


const char* timestamp_rejection_reason(const Timestamp& ts)
{
    if (!ts.valid)
        return "invalid_flag";
    if (!clock_domain_is_valid(ts.domain))
        return "invalid_clock_domain";
    if (!timestamp_marker_is_known(ts.marker))
        return "unknown_marker";
    if (!timestamp_source_is_known(ts.source) || ts.source == TimestampSource::Unknown)
        return "unknown_source";
    if (!clock_domain_ticks_in_range(ts.domain, ts.ticks))
        return "ticks_out_of_domain_range";
    if (ts.applied_corrections > kCorrectionAll)
        return "unknown_correction_bit";
    if (ts.frac_num == 0) {
        if (ts.frac_den != 0u)
            return "fraction_not_normalized";
    } else if (ts.frac_den == 0u) {
        return "fraction_not_normalized";
    } else if (ts.frac_den > detail::kFractionDenominatorMax) {
        return "fraction_denominator_too_large";
    } else if (ts.frac_num < 0 ||
               static_cast<int64_t>(ts.frac_num) >= static_cast<int64_t>(ts.frac_den)) {
        return "fraction_not_normalized";
    } else if (detail::ts_gcd(static_cast<int64_t>(ts.frac_num),
                              static_cast<int64_t>(ts.frac_den)) != 1) {
        return "fraction_not_reduced";
    }

    if (!timestamp_corrections_satisfied(ts.marker, ts.applied_corrections))
        return "corrections_not_applied";
    if (ts.marker != TimestampMarker::RmarkerRx &&
        ts.marker != TimestampMarker::RmarkerTx)
        return "marker_not_protocol_rmarker";
    if (!timestamp_source_is_measurement(ts.source))
        return "source_not_a_measurement";
    if (ts.calibration_id.empty())
        return "calibration_missing";
    if (ts.source != TimestampSource::HardwareMeasured)
        return "source_not_hardware_measured";
    return "accepted";
}


TimeInterval timestamp_interval(const Timestamp& later, const Timestamp& earlier)
{
    TimeInterval r;
    r.later_marker = later.marker;
    r.earlier_marker = earlier.marker;

    // 1. structural validity, including fraction normalization
    const bool frac_bad = !timestamp_fraction_is_normalized(later) ||
                          !timestamp_fraction_is_normalized(earlier);
    if (!timestamp_is_self_consistent(later) || !timestamp_is_self_consistent(earlier)) {
        r.status = frac_bad ? TimeIntervalStatus::FractionNotNormalized
                            : TimeIntervalStatus::InvalidTimestamp;
        return r;
    }

    // 2-4. same counter?  name, tick rate and wrap width are the identity.
    if (later.domain.name != earlier.domain.name) {
        r.status = TimeIntervalStatus::DomainNameMismatch;
        return r;
    }
    if (std::memcmp(&later.domain.tick_rate_hz,
                    &earlier.domain.tick_rate_hz,
                    sizeof(double)) != 0) {
        r.status = TimeIntervalStatus::TickRateMismatch;
        return r;
    }
    if (later.domain.timestamp_bits != earlier.domain.timestamp_bits) {
        r.status = TimeIntervalStatus::WrapWidthMismatch;
        return r;
    }

    // 5. same epoch?  a reboot / time reset makes older ticks meaningless.
    if (!clock_domain_same_epoch(later.domain, earlier.domain)) {
        r.status = TimeIntervalStatus::EpochMismatch;
        return r;
    }
    r.domain = later.domain;

    // 6. marker semantics: the two markers must denote a difference the
    //    protocol actually defines.  A UhdRxFirstIqSample is never silently
    //    equated with an RmarkerRx.
    if (!timestamp_marker_interval_allowed(later.marker, earlier.marker)) {
        r.status = TimeIntervalStatus::MarkerMismatch;
        return r;
    }

    // 7. both timestamps must be faithful instances of their marker
    //    (REQ-TIME-02): no interval across an unapplied mapping.
    if (!timestamp_corrections_satisfied(later.marker, later.applied_corrections) ||
        !timestamp_corrections_satisfied(earlier.marker, earlier.applied_corrections)) {
        r.status = TimeIntervalStatus::CorrectionsNotApplied;
        return r;
    }

    // 8. modular tick difference with wrap resolution / ambiguity refusal
    const TickDelta d = raw_tick_delta(r.domain, later.ticks, earlier.ticks);
    if (d.status != TimeIntervalStatus::Ok) {
        r.status = d.status;
        return r;
    }
    r.wrapped = d.wrapped;

    // 9. exact sub-tick part: frac_later - frac_earlier, carried into the
    //    whole-tick count with floor semantics.
    int64_t whole = 0;
    if (!detail::ts_sub_i64(static_cast<int64_t>(d.ticks), 0, whole)) {
        r.status = TimeIntervalStatus::DurationOutOfRange;
        return r;
    }
    if (whole > static_cast<int64_t>(std::numeric_limits<int64_t>::max())) {
        r.status = TimeIntervalStatus::DurationOutOfRange;
        return r;
    }

    if (later.frac_den != 0u || earlier.frac_den != 0u) {
        int64_t num = 0;
        int64_t den = 0;
        if (later.frac_den == 0u) {
            num = -static_cast<int64_t>(earlier.frac_num);
            den = static_cast<int64_t>(earlier.frac_den);
        } else if (earlier.frac_den == 0u) {
            num = static_cast<int64_t>(later.frac_num);
            den = static_cast<int64_t>(later.frac_den);
        } else {
            // den = db*de <= 32767^2 < 2^31 and |num| < 2^32, so both are
            // exact in int64.
            const int64_t a = static_cast<int64_t>(later.frac_num);
            const int64_t b = static_cast<int64_t>(later.frac_den);
            const int64_t c = static_cast<int64_t>(earlier.frac_num);
            const int64_t e = static_cast<int64_t>(earlier.frac_den);
            num = a * e - c * b;
            den = b * e;
        }
        if (num < 0) {
            if (!detail::ts_sub_i64(whole, 1, whole)) {
                r.status = TimeIntervalStatus::DurationOutOfRange;
                return r;
            }
            num += den;
        }
        if (num < 0 || den <= 0) {
            r.status = TimeIntervalStatus::FractionNotNormalized;
            return r;
        }
        const int64_t g = detail::ts_gcd(num, den);
        num /= g;
        den /= g;
        if (den > static_cast<int64_t>(detail::kFractionDenominatorMax)) {
            // Cannot store the exact value without losing it -> refuse.
            r.status = TimeIntervalStatus::FractionNotNormalized;
            return r;
        }
        if (num == 0) {
            r.frac_num = 0;
            r.frac_den = 0u;
        } else {
            r.frac_num = static_cast<int32_t>(num);
            r.frac_den = static_cast<uint32_t>(den);
        }
        if (whole < 0) {
            // A negative interval is never a valid TWR interval; REQ-PROTO-03
            // forbids clamping it to zero and calling it a success.
            r.status = TimeIntervalStatus::DurationOutOfRange;
            r.frac_num = 0;
            r.frac_den = 0u;
            return r;
        }
    }

    r.ticks = whole;

    // 9b. N04: re-check the upper bound AFTER the sub-tick fraction has been
    //     merged.  `raw_tick_delta()` gated only the integer difference, so a
    //     positive fraction could push `whole + frac` strictly past P/2 while
    //     the integer gate still saw P/2.  Exactly P/2 stays legal (the caller
    //     supplied the sign, so the magnitude is unique there); only a value
    //     that EXCEEDS the API's stated bound is refused.
    const uint64_t wrap_period = clock_domain_wrap_period(r.domain);
    if (wrap_period != 0u && detail::ts_interval_exceeds_half_period(
                                 r.ticks, r.frac_num, r.frac_den, wrap_period / 2u)) {
        r.status = TimeIntervalStatus::WrapAmbiguous;
        r.ticks = 0;
        r.frac_num = 0;
        r.frac_den = 0u;
        return r;
    }

    // 10. nanosecond projection.  REQ-API-02 wants SI time at the API
    //     boundary; the tick-space fields above stay exact for sub-nanosecond
    //     domains such as DW UUS.  The projection is TRUNCATED toward zero and
    //     its exactness is decided by the SAME quotient/remainder that produces
    //     the integer (N03), never by a rate threshold and never by a second
    //     floating-point path that could disagree with the flag.
    r.is_whole_ticks = (r.frac_den == 0u);
    int64_t ns = 0;
    bool lossless = false;
    if (!detail::ts_ticks_to_nanos_exact(
            r.domain.tick_rate_hz, r.ticks, r.frac_num, r.frac_den, ns, lossless)) {
        // A rejected interval is EMPTY, not partially populated: no tick count
        // and no "whole ticks" claim survive a failed projection, so a caller
        // that ignores `status` sees a zeroed value rather than a tick count it
        // might mistake for a measurement.
        r.status = TimeIntervalStatus::DurationOutOfRange;
        r.ticks = 0;
        r.frac_num = 0;
        r.frac_den = 0u;
        r.is_whole_ticks = false;
        r.duration_is_lossless_ns = false;
        return r;
    }
    r.duration_is_lossless_ns = lossless;
    r.duration = Duration::from_nanos(ns);

    r.status = TimeIntervalStatus::Ok;
    return r;
}


ModularTickDistance modular_tick_distance(const Timestamp& a, const Timestamp& b)
{
    ModularTickDistance r;
    if (!timestamp_is_self_consistent(a) || !timestamp_is_self_consistent(b))
        return r;
    if (!clock_domain_is_comparable(a.domain, b.domain))
        return r;
    const uint64_t p = clock_domain_wrap_period(a.domain);
    if (p == 0u)
        return r; // a monotonic counter has no modular distance

    const uint64_t from = static_cast<uint64_t>(a.ticks);
    const uint64_t to = static_cast<uint64_t>(b.ticks);
    r.wrapped = to < from;
    r.ticks = r.wrapped ? (to + p - from) : (to - from);

    // Exact sub-tick part: b.frac - a.frac, as a common-denominator rational.
    // den == 0 means 0/1.  den_a*den_b <= 32767^2 fits in uint32, and the
    // numerator is bounded by the same, so int64 cannot overflow.
    const uint32_t ad = a.frac_den == 0u ? 1u : a.frac_den;
    const uint32_t bd = b.frac_den == 0u ? 1u : b.frac_den;
    if (ad > kMaxFractionDifferenceDenominator / bd)
        return r; // unreachable for legal Timestamps; refused, never rounded
    const uint32_t den = ad * bd;
    const int64_t an = static_cast<int64_t>(a.frac_num) * static_cast<int64_t>(bd);
    const int64_t bn = static_cast<int64_t>(b.frac_num) * static_cast<int64_t>(ad);
    int64_t num = bn - an; // in (-den, den)
    if (num < 0) {
        // A fraction below `a`'s makes the distance one whole tick SHORTER
        // plus the wrapped-up remainder.  Skipping this borrow is the M0 bug
        // that made 255.25 -> 0.75 look like a full tick instead of 1.5, and
        // the two answers sit on opposite sides of a budget boundary.
        // When the integer part is already 0, the downward borrow would go
        // negative: a whole period elapsed minus (a.frac - b.frac), which is
        // outside every budget by construction.
        if (r.ticks == 0) {
            r.ticks = static_cast<int64_t>(p) - 1;
        } else {
            --r.ticks;
        }
        num += den;
    }
    const int64_t g = detail::ts_gcd(num, den);
    num /= g;
    if (num == 0) {
        r.frac_num = 0;
        r.frac_den = 0u;
    } else {
        r.frac_num = static_cast<int32_t>(num);
        r.frac_den = static_cast<uint32_t>(den / g);
    }
    r.ok = true;
    return r;
}


bool timestamp_from_json_string(const std::string& json, Timestamp& out)
{
    detail::JsonCursor c{ json.c_str(), json.c_str() + json.size() };
    Timestamp ts;
    int64_t frac_num_tmp = 0;
    int64_t frac_den_tmp = 0;
    if (!c.eat('{'))
        return false;

    if (!c.key("ticks") || !c.integer(ts.ticks))
        return false;
    if (!c.eat(','))
        return false;
    if (!c.key("frac_num") || !c.integer(frac_num_tmp))
        return false;
    if (!c.eat(','))
        return false;
    if (!c.key("frac_den") || !c.integer(frac_den_tmp))
        return false;
    if (frac_num_tmp < 0 || frac_num_tmp > std::numeric_limits<int32_t>::max())
        return false;
    if (frac_den_tmp < 0 || frac_den_tmp > std::numeric_limits<uint32_t>::max())
        return false;
    ts.frac_num = static_cast<int32_t>(frac_num_tmp);
    ts.frac_den = static_cast<uint32_t>(frac_den_tmp);
    if (!c.eat(','))
        return false;

    if (!c.key("clock_domain"))
        return false;
    std::string domain_raw;
    if (!c.raw_object(domain_raw))
        return false;
    {
        detail::JsonCursor d{ domain_raw.c_str(),
                              domain_raw.c_str() + domain_raw.size() };
        if (!d.eat('{'))
            return false;
        if (!d.key("name") || !d.string(ts.domain.name))
            return false;
        if (!d.eat(','))
            return false;
        if (!d.key("tick_rate_hz") || !d.real(ts.domain.tick_rate_hz))
            return false;
        if (!d.eat(','))
            return false;
        int64_t epoch = 0;
        if (!d.key("epoch_id") || !d.integer(epoch))
            return false;
        if (epoch < 0)
            return false;
        ts.domain.epoch_id = static_cast<uint64_t>(epoch);
        if (!d.eat(','))
            return false;
        int64_t bits = 0;
        if (!d.key("timestamp_bits") || !d.integer(bits))
            return false;
        if (bits < 0 || bits > 63)
            return false;
        ts.domain.timestamp_bits = static_cast<uint32_t>(bits);
        if (!d.eat('}'))
            return false;
    }
    if (!c.eat(','))
        return false;

    {
        std::string m;
        if (!c.key("marker") || !c.string(m))
            return false;
        if (!timestamp_marker_from_string(m, ts.marker))
            return false;
    }
    if (!c.eat(','))
        return false;
    {
        std::string s;
        if (!c.key("source") || !c.string(s))
            return false;
        if (!timestamp_from_string(s, ts.source))
            return false;
    }
    if (!c.eat(','))
        return false;
    {
        int64_t bits = 0;
        if (!c.key("corrections") || !c.integer(bits))
            return false;
        if (bits < 0 || bits > static_cast<int64_t>(kCorrectionAll))
            return false;
        ts.applied_corrections = static_cast<uint32_t>(bits);
    }
    if (!c.eat(','))
        return false;
    {
        std::string flags;
        if (!c.key("corrections_flags") || !c.string(flags))
            return false;
        if (flags != timestamp_corrections_to_string(ts.applied_corrections))
            return false;
    }
    if (!c.eat(','))
        return false;
    if (!c.key("calibration_id") || !c.string(ts.calibration_id))
        return false;
    if (!c.eat(','))
        return false;
    if (!c.key("valid") || !c.boolean(ts.valid))
        return false;
    if (!c.eat('}'))
        return false;
    if (c.p != c.end)
        return false; // trailing garbage

    if (!clock_domain_is_valid(ts.domain))
        return false;
    out = ts;
    return true;
}


namespace detail {

bool clock_domain_ns_per_tick_exact(double tick_rate_hz,
                                    int64_t& out_num,
                                    int64_t& out_den)
{
    out_num = 0;
    out_den = 0;
    if (!std::isfinite(tick_rate_hz) || tick_rate_hz <= 0.0)
        return false;

    // rate == mant * 2^exp exactly (frexp/ldexp are exact for binary FP).
    int be = 0;
    const double mant = std::ldexp(std::frexp(tick_rate_hz, &be), 53);
    if (!(mant >= 4503599627370496.0 && mant < 9007199254740992.0))
        return false; // not the expected normal form
    uint64_t m = static_cast<uint64_t>(mant);
    int e = be - 53;         // rate == m * 2^e
    while ((m & 1u) == 0u) { // absorb the powers of two into e; m ends up odd
        m >>= 1;
        ++e;
    }

    // ns_per_tick = 1e9 / (m * 2^e) = (2^9 * 1953125) / (m * 2^e)
    const int k = detail::kNsPerSecondPow2 - e;
    int64_t num = 0;
    int64_t den = 0;
    if (k >= 0) {
        if (k > 42) // 1953125 * 2^43 > INT64_MAX
            return false;
        num = detail::kNsPerSecondOdd * (static_cast<int64_t>(1) << k);
        den = static_cast<int64_t>(m);
    } else {
        const int s = -k;
        if (s >= 62 || m > (static_cast<uint64_t>(INT64_MAX) >> s))
            return false;
        num = detail::kNsPerSecondOdd;
        den = static_cast<int64_t>(m) << s;
    }
    const int64_t g = detail::ts_gcd(num, den);
    out_num = num / g;
    out_den = den / g;
    return true;
}


bool ts_ticks_to_nanos_exact(double tick_rate_hz,
                             int64_t ticks,
                             int32_t frac_num,
                             uint32_t frac_den,
                             int64_t& out_ns,
                             bool& out_lossless)
{
    out_ns = 0;
    out_lossless = false;
    if (ticks < 0 || frac_num < 0 || (frac_den == 0u && frac_num != 0))
        return false;

    int64_t np = 0;
    int64_t dp = 0;
    if (!clock_domain_ns_per_tick_exact(tick_rate_hz, np, dp))
        return false;

    using u128 = unsigned __int128;
    const u128 fden = (frac_den == 0u) ? u128{ 1 } : u128{ frac_den };
    const u128 x = static_cast<u128>(static_cast<uint64_t>(ticks)) * fden +
                   static_cast<u128>(static_cast<uint64_t>(frac_num));
    const u128 den = fden * static_cast<u128>(static_cast<uint64_t>(dp));
    // `x` and `den` are bounded by (2^63 * 2^15) and (2^63 * 2^15); the product
    // `x*np` is the only term that can leave 128 bits, and it is guarded.
    if (np != 0) {
        const u128 umax = ~static_cast<u128>(0);
        if (x != 0u && static_cast<u128>(static_cast<uint64_t>(np)) > umax / x)
            return false;
    }
    const u128 numer = x * static_cast<u128>(static_cast<uint64_t>(np));
    const u128 q = numer / den;
    const u128 rem = numer % den;
    out_lossless = (rem == 0u);
    if (q > static_cast<u128>(std::numeric_limits<int64_t>::max()))
        return false;
    out_ns = static_cast<int64_t>(q);
    return true;
}


} // namespace detail


bool encode_into(const Frame& f,
                 const FrameProfile& p,
                 uint8_t* out,
                 size_t out_capacity,
                 size_t& written,
                 std::string& error,
                 FrameError* code)
{
    written = 0;

    if (out == nullptr) {
        set_error(error, code, FrameError::NullPointer, "null output buffer");
        return false;
    }
    if (!frame_profile_validate(p, error, code))
        return false;

    // Version must match the profile; an unknown layout is never guessed at.
    if (f.version != p.version) {
        set_error(
            error, code, FrameError::BadVersion, "frame version != profile version");
        return false;
    }

    const FrameLayout* layout = frame_layout(f.function_code);
    if (layout == nullptr) {
        // Report, or an out-of-range enum value.
        if (f.function_code == FrameType::Report) {
            set_error(error,
                      code,
                      FrameError::NotImplemented,
                      "report frame type is reserved; phase 1 emits the "
                      "three-message exchange only (REQ-PROTO-05)");
        } else {
            set_error(error, code, FrameError::BadFunctionCode, "unknown function code");
        }
        return false;
    }

    if (flags_has_reserved_bits(f.flags)) {
        set_error(error,
                  code,
                  FrameError::ReservedFlags,
                  "flags bits 2..15 are reserved and must be zero");
        return false;
    }
    if (flags_sts_present(f.flags)) {
        set_error(error,
                  code,
                  FrameError::StsUnsupported,
                  "sts_present is set; STS is out of scope for phase 1");
        return false;
    }

    const size_t tsb = timestamp_bytes(p);
    const size_t n_ts = layout->timestamp_count;
    // frame_length_for() is the ONE place a frame's length is computed; this
    // function must not restate the arithmetic (R3: no second geometry).
    const size_t need = frame_length_for(f.function_code, p);

    // A value that does not fit the declared field is refused, never
    // truncated: truncation would silently change a measurement instant.
    const uint64_t tmax = timestamp_max_value(p);
    for (size_t i = 0; i < n_ts; ++i) {
        if (f.timestamps[i] > tmax) {
            set_error(error,
                      code,
                      FrameError::TimestampOutOfRange,
                      "timestamp does not fit the profile field width");
            return false;
        }
    }
    // An unused slot must be zero so a frame's bytes are a function of its
    // declared content alone and a stale instant can never ride along.
    for (size_t i = n_ts; i < kMaxTimestamps; ++i) {
        if (f.timestamps[i] != 0) {
            set_error(error,
                      code,
                      FrameError::ReservedTimestampNonZero,
                      "unused timestamp slot is non-zero");
            return false;
        }
    }

    if (need > out_capacity) {
        set_error(error,
                  code,
                  FrameError::BufferTooSmall,
                  "caller buffer smaller than the frame");
        return false;
    }
    if (need + kFrameFcsBytes > p.max_psdu_bytes) {
        set_error(error,
                  code,
                  FrameError::PsduTooLarge,
                  "frame plus FCS exceeds the profile max PSDU");
        return false;
    }

    std::memset(out, 0, need);
    out[kOffVersion] = f.version;
    out[kOffFunctionCode] = static_cast<uint8_t>(f.function_code);
    put_u16_le(out + kOffSessionId, f.session_id);
    put_u16_le(out + kOffSeq, f.seq);
    put_u16_le(out + kOffPanId, f.pan_id);
    put_u16_le(out + kOffSrcAddr, f.src_addr);
    put_u16_le(out + kOffDstAddr, f.dst_addr);
    put_u16_le(out + kOffFlags, f.flags);
    for (size_t i = 0; i < n_ts; ++i) {
        put_uint_le(out + kOffTimestamps + i * tsb, f.timestamps[i], tsb);
    }

    written = need;
    if (code != nullptr)
        *code = FrameError::None;
    error.clear();
    return true;
}


bool encode(const Frame& f,
            const FrameProfile& p,
            std::vector<uint8_t>& out,
            std::string& error,
            FrameError* code)
{
    out.clear();
    std::string ignored;
    if (!frame_profile_validate(p, ignored, code)) {
        error = ignored;
        return false;
    }
    const FrameType t = f.function_code;
    if (!frame_type_is_implemented(t)) {
        set_error(error,
                  code,
                  t == FrameType::Report ? FrameError::NotImplemented
                                         : FrameError::BadFunctionCode,
                  t == FrameType::Report
                      ? "report frame type is reserved; phase 1 emits the "
                        "three-message exchange only (REQ-PROTO-05)"
                      : "unknown function code");
        return false;
    }
    const size_t need = frame_length_for(t, p);
    out.resize(need);
    size_t written = 0;
    if (!encode_into(f, p, out.data(), out.size(), written, error, code)) {
        out.clear();
        return false;
    }
    out.resize(written);
    return true;
}


bool decode(const uint8_t* data,
            size_t n,
            const FrameProfile& p,
            Frame& out,
            std::string& error,
            FrameError* code)
{
    out = Frame();
    out.version = 0; // invalid on purpose; see the comment above

    if (data == nullptr && n != 0) {
        set_error(error, code, FrameError::NullPointer, "null payload pointer");
        return false;
    }
    if (!frame_profile_validate(p, error, code))
        return false;

    if (n < kFrameHeaderBytes) {
        set_error(error, code, FrameError::Truncated, "shorter than the 14-byte header");
        return false;
    }

    if (data[kOffVersion] != p.version) {
        set_error(error, code, FrameError::BadVersion, "unknown frame version");
        return false;
    }

    FrameType type = FrameType::Poll;
    if (!frame_type_from_wire(data[kOffFunctionCode], type)) {
        set_error(error, code, FrameError::BadFunctionCode, "unknown function code");
        return false;
    }
    const FrameLayout* layout = frame_layout(type);
    if (layout == nullptr) {
        set_error(error,
                  code,
                  FrameError::NotImplemented,
                  "report frame type is reserved; phase 1 decodes the "
                  "three-message exchange only (REQ-PROTO-05)");
        return false;
    }

    const uint16_t flags = get_u16_le(data + kOffFlags);
    if (flags_has_reserved_bits(flags)) {
        set_error(error,
                  code,
                  FrameError::ReservedFlags,
                  "flags bits 2..15 are reserved and must be zero");
        return false;
    }
    if (flags_sts_present(flags)) {
        set_error(error,
                  code,
                  FrameError::StsUnsupported,
                  "sts_present is set; STS is out of scope for phase 1");
        return false;
    }

    const size_t tsb = timestamp_bytes(p);
    // Same single place as encode_into(): a frame's length is computed once.
    const size_t need = frame_length_for(type, p);
    if (n < need) {
        set_error(error,
                  code,
                  FrameError::Truncated,
                  "timestamp block is shorter than the function code requires");
        return false;
    }
    if (n > need) {
        set_error(error,
                  code,
                  FrameError::LengthMismatch,
                  "buffer is longer than the function code requires; a PSDU "
                  "whose FCS was not stripped is the usual cause");
        return false;
    }
    if (need + kFrameFcsBytes > p.max_psdu_bytes) {
        set_error(error,
                  code,
                  FrameError::PsduTooLarge,
                  "frame plus FCS exceeds the profile max PSDU");
        return false;
    }

    Frame f;
    f.version = data[kOffVersion];
    f.function_code = type;
    f.session_id = get_u16_le(data + kOffSessionId);
    f.seq = get_u16_le(data + kOffSeq);
    f.pan_id = get_u16_le(data + kOffPanId);
    f.src_addr = get_u16_le(data + kOffSrcAddr);
    f.dst_addr = get_u16_le(data + kOffDstAddr);
    f.flags = flags;
    for (size_t i = 0; i < kMaxTimestamps; ++i) {
        f.timestamps[i] = 0;
    }
    for (size_t i = 0; i < layout->timestamp_count; ++i) {
        const uint64_t v = get_uint_le(data + kOffTimestamps + i * tsb, tsb);
        if (!timestamp_fits(p, v)) {
            // Unreachable for a whole-byte field, kept as a standing
            // invariant check: a value that cannot be represented
            // unambiguously in this profile must never leave the codec.
            set_error(error,
                      code,
                      FrameError::TimestampOutOfRange,
                      "decoded timestamp is not representable in this profile");
            return false;
        }
        f.timestamps[i] = v;
    }

    out = f;
    if (code != nullptr)
        *code = FrameError::None;
    error.clear();
    return true;
}


bool mac_payload_from_psdu(const uint8_t* psdu,
                           size_t n,
                           const FrameProfile& p,
                           const uint8_t*& payload,
                           size_t& payload_bytes,
                           std::string& error,
                           FrameError* code)
{
    payload = nullptr;
    payload_bytes = 0;
    if (psdu == nullptr && n != 0) {
        set_error(error, code, FrameError::NullPointer, "null psdu pointer");
        return false;
    }
    if (n < kFrameHeaderBytes + kFrameFcsBytes) {
        set_error(error, code, FrameError::Truncated, "psdu shorter than header + FCS");
        return false;
    }
    const size_t m = n - kFrameFcsBytes;
    if (m + kFrameFcsBytes > p.max_psdu_bytes) {
        set_error(
            error, code, FrameError::PsduTooLarge, "psdu exceeds the profile max PSDU");
        return false;
    }
    payload = psdu;
    payload_bytes = m;
    if (code != nullptr)
        *code = FrameError::None;
    error.clear();
    return true;
}


bool frame_profile_validate(const FrameProfile& p, std::string& error, FrameError* code)
{
    if (p.version != kFrameVersion) {
        set_error(error, code, FrameError::BadProfile, "profile version is not 0x01");
        return false;
    }
    if (p.timestamp_bits == 0 || p.timestamp_bits > 64 || (p.timestamp_bits % 8u) != 0) {
        set_error(error,
                  code,
                  FrameError::BadProfile,
                  "timestamp_bits must be a multiple of 8 in [8, 64]");
        return false;
    }
    // A NaN tick rate fails the `> 0.0` test, so no <cmath> is needed.
    if (!(p.timestamp_unit_hz > 0.0) || p.timestamp_unit_hz > 1.0e15) {
        set_error(error,
                  code,
                  FrameError::BadProfile,
                  "timestamp_unit_hz must be finite and in (0, 1e15]");
        return false;
    }
    if (p.max_psdu_bytes < kFrameHeaderBytes || p.max_psdu_bytes > kPhyMaxPsduBytes) {
        set_error(
            error, code, FrameError::BadProfile, "max_psdu_bytes must be in [14, 127]");
        return false;
    }
    if (!p.fcs_appended_by_modulation_layer) {
        set_error(error,
                  code,
                  FrameError::BadProfile,
                  "this codec never emits an FCS; it is appended by the HRP "
                  "modulation layer (UwbHrpPacketSource append_fcs=true)");
        return false;
    }
    if (code != nullptr)
        *code = FrameError::None;
    error.clear();
    return true;
}


bool hex_to_bytes(const std::string& hex, std::vector<uint8_t>& out)
{
    out.clear();
    if (hex.size() % 2 != 0)
        return false;
    out.reserve(hex.size() / 2);
    for (size_t i = 0; i < hex.size(); i += 2) {
        int hi = -1;
        int lo = -1;
        for (int pass = 0; pass < 2; ++pass) {
            const char c = hex[i + static_cast<size_t>(pass)];
            int v = -1;
            if (c >= '0' && c <= '9')
                v = c - '0';
            else if (c >= 'a' && c <= 'f')
                v = c - 'a' + 10;
            else if (c >= 'A' && c <= 'F')
                v = c - 'A' + 10;
            if (v < 0)
                return false;
            if (pass == 0)
                hi = v;
            else
                lo = v;
        }
        out.push_back(static_cast<uint8_t>((hi << 4) | lo));
    }
    return true;
}


std::string preamble_length_unsupported_reason(uint16_t symbols)
{
    const std::string n = std::to_string(symbols);
    if (symbols == 32)
        return "preamble length " + n +
               " SYNC is refused: mod::encode_phr19 advertises preamble-duration index 1 "
               "(= 64) "
               "for it, so the PHR does not self-describe a 32-symbol preamble";
    if (symbols == 128 || symbols == 256 || symbols == 512)
        return "preamble length " + n +
               " SYNC is refused by the CURRENT DEMODULATOR, not by any hardware: "
               "uwb_demod_core.h re-measures only the last max(cfo_min_fit_repetitions, "
               "40) SYNCs "
               "and synthesises the earlier peaks at phase 0, so the zero-phase points "
               "that leak "
               "into the least-squares CFO fit number max(0, reps - 64) -- measured at "
               "128 SYNC as "
               "an injected 20 kHz offset coming back as 3791 Hz; and mod::encode_phr19 "
               "maps " +
               n +
               " to preamble-duration index 1 (= 64), so the PHR does not self-describe "
               "it. "
               "Qorvo's API accepts these lengths (DW3xxx lists 128/256/512 as "
               "non-standard "
               "preamble lengths), so this is a limit of this software path, not a "
               "statement about "
               "what a DW1000 or DW3000 can transmit";
    if (symbols == 1024)
        return "preamble length " + n +
               " SYNC is PHR-legal (preamble-duration index 2) but refused by the "
               "current "
               "demodulator for the same CFO-fit bias as 128/256/512; a software limit, "
               "not a "
               "hardware one";
    if (symbols == 2048)
        return "preamble length " + n +
               " SYNC exceeds the work-grid TX buffer this build sizes for one frame "
               "(mod::kMaxHrpTxSamples); a length this large is a buffer and scheduling "
               "decision, not a hardware capability statement";
    if (symbols == 1 || symbols == 2)
        return "preamble length " + n +
               " SYNC leaves fewer than the 4 measured peaks uwb_demod_core.h stage_cfo "
               "needs";
    if (symbols == 4 || symbols == 8)
        return "preamble length " + n +
               " SYNC is shorter than the 10 initial repetitions uwb_demod_core.h skips "
               "before "
               "the CIR search";
    return "preamble length " + n +
           " is not a SYNC repetition count this build implements: the enumeration is 16 "
           "or 64, "
           "bounded by this repository's demodulator/CFO path (see the PreambleLength "
           "comment). "
           "That is a software limit; the Qorvo parts' own API accepts 64 through 2048";
}


const char* time_reference_event_to_string(TimeReferenceEvent v)
{
    switch (v) {
    case TimeReferenceEvent::PollTransmitRmarker:
        return "poll_tx_rmarker";
    case TimeReferenceEvent::PollReceiveRmarker:
        return "poll_rx_rmarker";
    case TimeReferenceEvent::ResponseTransmitRmarker:
        return "response_tx_rmarker";
    case TimeReferenceEvent::ResponseReceiveRmarker:
        return "response_rx_rmarker";
    case TimeReferenceEvent::FinalTransmitRmarker:
        return "final_tx_rmarker";
    case TimeReferenceEvent::FinalReceiveRmarker:
        return "final_rx_rmarker";
    case TimeReferenceEvent::ReportTransmitRmarker:
        return "report_tx_rmarker";
    case TimeReferenceEvent::ReportReceiveRmarker:
        return "report_rx_rmarker";
    case TimeReferenceEvent::FrameTail:
        return "frame_tail";
    case TimeReferenceEvent::RxEnable:
        return "rx_enable";
    case TimeReferenceEvent::HostMonotonic:
        return "host_monotonic";
    }
    return "invalid";
}


bool data_rate_from_string(const std::string& s, DataRate& out)
{
    if (s == "850k") {
        out = DataRate::R850k;
        return true;
    }
    if (s == "6p8m") {
        out = DataRate::R6p8M;
        return true;
    }
    if (s == "27m") {
        out = DataRate::R27M;
        return true;
    }
    if (s == "7p8m") {
        out = DataRate::R7p8M;
        return true;
    }
    if (s == "27p2m") {
        out = DataRate::R27p2M;
        return true;
    }
    if (s == "6p8m_hprf") {
        out = DataRate::R6p8M_hprf;
        return true;
    }
    return false;
}


bool first_path_algorithm_from_string(const std::string& s, FirstPathAlgorithm& out)
{
    if (s == "leading_edge") {
        out = FirstPathAlgorithm::LeadingEdge;
        return true;
    }
    if (s == "peak") {
        out = FirstPathAlgorithm::Peak;
        return true;
    }
    if (s == "energy_centroid") {
        out = FirstPathAlgorithm::EnergyCentroid;
        return true;
    }
    if (s == "interpolated_peak") {
        out = FirstPathAlgorithm::InterpolatedPeak;
        return true;
    }
    if (s == "strongest_cluster") {
        out = FirstPathAlgorithm::StrongestCluster;
        return true;
    }
    return false;
}


} // namespace twr
} // namespace openuwb
