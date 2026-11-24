#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Design the 65/48 (737.28 -> 998.4 MS/s) 250 MHz channelisation filters.

Task: docs/OpenCode执行方案_65_48重采样250MHz信道滤波器设计与验收_2026-11-24.md

Prototype (virtual) rate
------------------------
    R_IN  = 737.28 MHz
    R_OUT = 998.4  MHz = 65/48 * R_IN
    R_V   = 65 * 737.28 MHz = 47.9232 GHz          (L = 65, M = 48)

Three-band don't-care Parks-McClellan (``scipy.signal.remez``)
--------------------------------------------------------------
    Band0  passband            [0, 250 MHz]                 gain 1
    Band1  in-Nyquist OOB stop [f_oob, 368.64 MHz]          gain 0
    gap    don't-care          (368.64 MHz, 487.28 MHz)     NOT a remez band
    Band2  anti-image stop     [487.28 MHz, R_V/2]          gain 0

The gap is only *checked* after quantisation: it holds no wanted signal and no
first image of the wanted band (the first image starts at R_IN - 250 MHz), so it
is left unconstrained; a large overshoot there is reported, never hidden.

Coefficients are float32, odd length T = 65*Hpad - 1, exactly symmetric, and
DC-normalised so that ``sum(h32) == 65`` within 1e-6 (relative).  Every figure
below is measured on the *quantised* coefficients.

Usage
-----
    python3 tools/dsp/design_resampler_65_48_ch250.py
    python3 tools/dsp/design_resampler_65_48_ch250.py --only ch250_eq32 --dry-run
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np

try:
    import scipy
    from scipy.signal import remez, upfirdn
except Exception as exc:  # pragma: no cover - environment guard
    sys.stderr.write("scipy is required: %s\n" % exc)
    raise

# ---------------------------------------------------------------- constants ---
R_IN = 737.28e6
R_OUT = 998.4e6
L, M = 65, 48
FV = L * R_IN                     # 47.9232 GHz
FV_HALF = FV / 2.0

F_PASS = 250.0e6                  # passband edge
BAND1_HI = R_IN / 2.0             # 368.64 MHz = input Nyquist
BAND2_LO = R_IN - F_PASS          # 487.28 MHz = first image start
DONTCARE = (BAND1_HI, BAND2_LO)

DC_TARGET = float(L)
DC_REL_TOL = 1e-6

# weight / grid-density search space (from the task book)
W_OOB = (0.5, 1.0, 2.0, 4.0, 8.0)
W_IMG = (4.0, 8.0, 16.0, 32.0, 64.0, 128.0)
GRID_DENSITY = (16, 32, 48)

SEARCH_NFFT = 1 << 20             # screening grid
FINAL_NFFT = 1 << 21              # acceptance grid

REMEZ_MAXITER = 400
EXTENSION_HZ = 10e6               # allowed upward f_oob extension
SELFCHECK_SEED = 20261124
SELFCHECK_LEN = 4096

F_PROBE = np.array([276.0e6, 300.0e6, 320.0e6])
RMS_BAND = (276.0e6, 322.0e6)

PROFILES = [
    dict(name="ch250_eq16", hpad=16, f_lo_mhz=330, f_hi_mhz=350,
         ripple_db=0.20, oob_max_db=-35.0, img_max_db=-63.0),
    dict(name="ch250_eq24", hpad=24, f_lo_mhz=300, f_hi_mhz=318,
         ripple_db=0.20, oob_max_db=-38.0, img_max_db=-65.0),
    dict(name="ch250_eq32", hpad=32, f_lo_mhz=282, f_hi_mhz=298,
         ripple_db=0.20, oob_max_db=-40.0, img_max_db=-66.0),
    dict(name="ch250_eq48", hpad=48, f_lo_mhz=268, f_hi_mhz=280,
         ripple_db=0.15, oob_max_db=-42.0, img_max_db=-68.0),
]


# ---------------------------------------------------------------- helpers ---

def dft_at(h64, freqs):
    """Exact DTFT magnitude at arbitrary frequencies (Hz)."""
    n = np.arange(h64.size, dtype=np.float64)
    return np.abs(np.exp(-2j * np.pi * np.outer(np.asarray(freqs, float), n) / FV) @ h64)


def dense_response(h32, nfft):
    h64 = h32.astype(np.float64)
    mag = np.abs(np.fft.rfft(h64, nfft))
    f = np.arange(mag.size, dtype=np.float64) * (FV / nfft)
    return f, mag


def evaluate(h32, f_oob, nfft=FINAL_NFFT):
    """Full post-quantisation metric set for one design."""
    h64 = h32.astype(np.float64)
    f, mag = dense_response(h32, nfft)

    # passband [0, 250 MHz] including exact endpoints
    crit_pb = dft_at(h64, [0.0, F_PASS])
    m_pb = f <= F_PASS
    pb_max = float(max(mag[m_pb].max(), crit_pb.max()))
    pb_min = float(min(mag[m_pb].min(), crit_pb.min()))

    # Band1 [f_oob, 368.64 MHz]
    crit_b1 = dft_at(h64, [f_oob, BAND1_HI])
    m_b1 = (f >= f_oob) & (f <= BAND1_HI)
    b1_max = float(max(mag[m_b1].max(), crit_b1.max()))
    b1_rms = float(np.sqrt(np.mean(mag[m_b1] ** 2)))

    # Band2 images [487.28 MHz, FV/2]
    crit_b2 = dft_at(h64, [BAND2_LO, FV_HALF])
    m_b2 = f >= BAND2_LO
    b2_max = float(max(mag[m_b2].max(), crit_b2.max()))

    # don't-care open interval (368.64, 487.28)
    m_dc = (f > DONTCARE[0]) & (f < DONTCARE[1])
    dc_vals = mag[m_dc]
    dc_idx = int(np.argmax(dc_vals))
    dc_max = float(dc_vals[dc_idx])
    dc_f = float(f[m_dc][dc_idx])

    # probe points + [276, 322] RMS
    probes = dft_at(h64, F_PROBE)
    m_rms = (f >= RMS_BAND[0]) & (f <= RMS_BAND[1])
    rms_276_322 = float(np.sqrt(np.mean(mag[m_rms] ** 2)))

    db = lambda x: float(20.0 * np.log10(x / pb_min))
    return dict(
        nfft=int(nfft),
        f_oob=float(f_oob),
        dc_sum=float(h64.sum()),
        passband_max=pb_max,
        passband_min=pb_min,
        passband_ripple_db=float(20.0 * np.log10(pb_max / pb_min)),
        oob_max_db=db(b1_max),
        oob_rms_db=db(b1_rms),
        at_276mhz_db=db(probes[0]),
        at_300mhz_db=db(probes[1]),
        at_320mhz_db=db(probes[2]),
        dontcare_max_db=db(dc_max),
        dontcare_freq_mhz=dc_f / 1e6,
        img_max_db=db(b2_max),
        rms_276_322_db=db(rms_276_322),
    )


def violations(res, prof):
    """List of hard-spec violations (human readable)."""
    out = []
    if res["passband_ripple_db"] > prof["ripple_db"]:
        out.append("passband_ripple %.4f > %.4f dB"
                   % (res["passband_ripple_db"], prof["ripple_db"]))
    if res["oob_max_db"] > prof["oob_max_db"]:
        out.append("oob_max %.3f > %.3f dB"
                   % (res["oob_max_db"], prof["oob_max_db"]))
    if res["img_max_db"] > prof["img_max_db"]:
        out.append("img_max %.3f > %.3f dB"
                   % (res["img_max_db"], prof["img_max_db"]))
    return out


def meets(res, prof):
    return not violations(res, prof)


def design_taps(T, f_oob, w_oob, w_img, grid_density):
    """Remez design -> symmetric, DC-normalised float32 (sum == 65)."""
    bands = [0.0, F_PASS / FV, f_oob / FV, BAND1_HI / FV, BAND2_LO / FV, 0.5]
    h = remez(int(T), bands, [1.0, 0.0, 0.0],
              weight=[1.0, float(w_oob), float(w_img)], fs=1.0,
              maxiter=REMEZ_MAXITER, grid_density=int(grid_density))
    if not np.all(np.isfinite(h)):
        raise FloatingPointError("remez returned non-finite coefficients")
    h = np.asarray(h, dtype=np.float64)
    h = 0.5 * (h + h[::-1])              # enforce exact linear-phase symmetry
    h = h * (DC_TARGET / h.sum())        # float64 normalisation

    h32 = h.astype("<f4")
    # one float64 correction pass on the quantised values
    s = float(h32.astype(np.float64).sum())
    h32 = (h32.astype(np.float64) * (DC_TARGET / s)).astype("<f4")
    return h32


def try_design(T, f_oob, w_oob, w_img, grid_density, prof, nfft):
    try:
        h32 = design_taps(T, f_oob, w_oob, w_img, grid_density)
    except Exception as exc:            # noqa: BLE001 - remez may not converge
        return dict(ok=False, error=str(exc)[:160], w_oob=w_oob,
                    w_img=w_img, grid_density=grid_density, f_oob=f_oob)
    res = evaluate(h32, f_oob, nfft=nfft)
    res.update(ok=True, w_oob=float(w_oob), w_img=float(w_img),
               grid_density=int(grid_density), f_oob=float(f_oob), h32=h32)
    res["spec_met"] = meets(res, prof)
    res["violations"] = violations(res, prof)
    return res


def select_candidate(cands):
    """Pareto: f_oob is already fixed; then larger [276,322] RMS attenuation,
    then smaller ripple; near-equivalents fall back to the lower don't-care peak.
    """
    best = min(cands, key=lambda c: (c["rms_276_322_db"],
                                     c["passband_ripple_db"]))
    equiv = [c for c in cands
             if c["rms_276_322_db"] <= best["rms_276_322_db"] + 0.2
             and c["passband_ripple_db"] <= best["passband_ripple_db"] + 0.01]
    return min(equiv, key=lambda c: (c["dontcare_max_db"],
                                     c["rms_276_322_db"],
                                     c["passband_ripple_db"],
                                     c["w_oob"], c["w_img"], c["grid_density"]))


def sha256(h32):
    return hashlib.sha256(np.ascontiguousarray(h32, dtype="<f4").tobytes()).hexdigest()


# ------------------------------------------------------------------ search ---

def run_profile(prof):
    name = prof["name"]
    T = L * prof["hpad"] - 1
    t_start = time.time()
    log = []

    f_lo = int(prof["f_lo_mhz"])
    f_hi = int(prof["f_hi_mhz"])
    f_ext_hi = f_hi + int(round(EXTENSION_HZ / 1e6))

    def scan_f(f_mhz, weights=W_OOB, wimg=W_IMG, gd_list=(32,), nfft=SEARCH_NFFT):
        """Return list of meeting candidates at one f_oob over the weight grid."""
        found = []
        f_oob = f_mhz * 1e6
        for gd in gd_list:
            for wo in weights:
                for wi in wimg:
                    r = try_design(T, f_oob, wo, wi, gd, prof, nfft)
                    if r["ok"] and r["spec_met"]:
                        found.append(r)
        return found

    # --- stage 1: find the lowest feasible f_oob (1 MHz step, gd=32) ---------
    f_star = None
    stage1_cands = None
    for f_mhz in range(f_lo, f_ext_hi + 1):
        c = scan_f(f_mhz)
        log.append("stage1 f_oob=%dMHz meeting=%d" % (f_mhz, len(c)))
        if c:
            f_star = f_mhz
            stage1_cands = c
            break

    if f_star is None:
        # no design meets the hard spec even after extension -> best effort
        f_star = f_ext_hi
        stage1_cands = []
        # choose the design with the fewest / smallest violations at f_star
        best = None
        for wo in W_OOB:
            for wi in W_IMG:
                for gd in GRID_DENSITY:
                    r = try_design(T, f_star * 1e6, wo, wi, gd, prof, FINAL_NFFT)
                    if not r["ok"]:
                        continue
                    score = (len(r["violations"]),
                             max([r["passband_ripple_db"] - prof["ripple_db"],
                                  r["oob_max_db"] - prof["oob_max_db"],
                                  r["img_max_db"] - prof["img_max_db"]]))
                    if best is None or score < best[0]:
                        best = (score, r)
        if best is not None:
            stage1_cands = [best[1]]
        chosen = best[1] if best is not None else None
        extended = f_star > f_hi
        return _finish(prof, T, chosen, f_star, extended, log, t_start,
                       search_note="no candidate met hard spec up to f_oob=%dMHz; "
                                   "best-effort closest design" % f_star)

    # --- stage 1b: guard against grid-density effects at f_star-1 ------------
    while f_star > f_lo:
        lower = scan_f(f_star - 1, gd_list=GRID_DENSITY, nfft=FINAL_NFFT)
        if not lower:
            break
        f_star -= 1
        stage1_cands = lower
        log.append("stage1b lowered to f_oob=%dMHz" % f_star)

    # --- stage 2: full weight x grid_density grid at f_star ------------------
    full = scan_f(f_star, gd_list=GRID_DENSITY, nfft=FINAL_NFFT)
    # if the dense acceptance grid disagrees with the screening grid, step up
    tries = 0
    while not full and f_star < f_ext_hi and tries < 12:
        f_star += 1
        tries += 1
        full = scan_f(f_star, gd_list=GRID_DENSITY, nfft=FINAL_NFFT)
        log.append("stage2 stepped up to f_oob=%dMHz meeting=%d"
                   % (f_star, len(full)))

    chosen = select_candidate(full) if full else None
    extended = f_star > f_hi
    note = ("full 3x5x6 grid evaluated at f_oob=%dMHz" % f_star)
    if extended:
        note += ("; f_oob extended by %dMHz beyond the %d-%dMHz band because no "
                 "weight/grid combination met the hard spec inside it"
                 % (f_star - f_hi, f_lo, f_hi))
    return _finish(prof, T, chosen, f_star, extended, log, t_start,
                   search_note=note)


def _finish(prof, T, chosen, f_star, extended, log, t_start, search_note):
    name = prof["name"]
    out = dict(name=name, hpad=int(prof["hpad"]), taps=int(T), f_star=f_star,
               extended=bool(extended), search_note=search_note, log=log,
               elapsed_s=round(time.time() - t_start, 1), chosen=chosen,
               prof=prof)
    return out


# ------------------------------------------------------------- self-check ---

def polyphase_reference(h64, x, L=L, M=M):
    """Hand-written y[m] = sum_k h[(m*M mod L)+L*k] * x[floor(m*M/L)-k]."""
    N = x.size
    T = h64.size
    lout = -(-((N - 1) * L + T) // M)
    y = np.zeros(lout, dtype=np.complex128)
    for m in range(lout):
        a = (m * M) % L
        q = (m * M) // L
        kmin = max(0, q - (N - 1))
        kmax = min((T - 1 - a) // L, q)
        for k in range(kmin, kmax + 1):
            y[m] += h64[a + L * k] * x[q - k]
    return y


def self_check(h32):
    rng = np.random.default_rng(SELFCHECK_SEED)
    x = (rng.standard_normal(SELFCHECK_LEN)
         + 1j * rng.standard_normal(SELFCHECK_LEN)).astype(np.complex128)
    h64 = h32.astype(np.float64)
    ref = upfirdn(h64, x, up=L, down=M)
    man = polyphase_reference(h64, x)
    n = min(ref.size, man.size)
    return dict(max_abs_diff=float(np.max(np.abs(ref[:n] - man[:n]))),
                length_upfirdn=int(ref.size), length_reference=int(man.size),
                length_match=bool(ref.size == man.size),
                input_len=SELFCHECK_LEN, seed=SELFCHECK_SEED)


# ---------------------------------------------------------------- outputs ---

def build_record(prof, res, h32, chk):
    T = int(h32.size)
    group_delay_virtual = (T - 1) / 2.0
    rec = dict(
        name=prof["name"],
        taps_file="taps_%s.txt" % prof["name"],
        hpad=int(prof["hpad"]),
        taps=T,
        f_pass_mhz=F_PASS / 1e6,
        f_oob_mhz=res["f_oob"] / 1e6,
        weight_oob=float(res["w_oob"]),
        weight_img=float(res["w_img"]),
        grid_density=int(res["grid_density"]),
        dc_sum=float(res["dc_sum"]),
        passband_ripple_db=float(res["passband_ripple_db"]),
        oob_max_db=float(res["oob_max_db"]),
        oob_rms_db=float(res["oob_rms_db"]),
        at_276mhz_db=float(res["at_276mhz_db"]),
        at_300mhz_db=float(res["at_300mhz_db"]),
        at_320mhz_db=float(res["at_320mhz_db"]),
        dontcare_max_db=float(res["dontcare_max_db"]),
        dontcare_freq_mhz=float(res["dontcare_freq_mhz"]),
        img_max_db=float(res["img_max_db"]),
        rms_276_322_db=float(res["rms_276_322_db"]),
        group_delay_virtual_samples=float(group_delay_virtual),
        group_delay_input_samples=float(group_delay_virtual / L),
        group_delay_output_samples=float(group_delay_virtual / M),
        self_check_upfirdn_max_abs_diff=float(chk["max_abs_diff"]),
        self_check_length_upfirdn=chk["length_upfirdn"],
        self_check_length_reference=chk["length_reference"],
        self_check_length_match=chk["length_match"],
        self_check_seed=chk["seed"],
        spec_met=bool(res["spec_met"]),
        search_extended=bool(res["f_oob"] / 1e6 > prof["f_hi_mhz"]),
        dontcare_warning=bool(res["dontcare_max_db"] > 1.0),
        sha256=sha256(h32),
        symmetry_max_abs_error=float(np.max(np.abs(h32.astype(np.float64)
                                                   - h32[::-1].astype(np.float64)))),
        specs=dict(passband_ripple_db_max=float(prof["ripple_db"]),
                   oob_max_db_max=float(prof["oob_max_db"]),
                   img_max_db_max=float(prof["img_max_db"]),
                   f_oob_range_mhz=[float(prof["f_lo_mhz"]),
                                    float(prof["f_hi_mhz"])]),
        dc_sum_rel_error=float(abs(res["dc_sum"] - DC_TARGET) / DC_TARGET),
    )
    if not rec["spec_met"]:
        rec["spec_misses"] = list(res["violations"])
    if rec["search_extended"]:
        rec["search_extended_reason"] = (
            "no weight/grid combination met the hard spec inside the "
            "%.0f-%.0f MHz f_oob band; extended upward (limit +%.0f MHz)"
            % (prof["f_lo_mhz"], prof["f_hi_mhz"], EXTENSION_HZ / 1e6))
    if rec["dontcare_warning"]:
        rec["dontcare_warning_note"] = (
            "don't-care gap (368.64, 487.28) MHz peaks at %.2f dB (relative to "
            "the passband gain); the gap is unconstrained (no wanted signal and "
            "no first image of the wanted band) but exceeds the <= +1 dB "
            "recommendation and is reported, not hidden." % rec["dontcare_max_db"])
    # reject non-finite values (JSON must stay strict)
    for k, v in list(rec.items()):
        if isinstance(v, float) and not np.isfinite(v):
            raise ValueError("non-finite metric %s=%r" % (k, v))
    return rec


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.join("testdata", "resampler_65_48"))
    ap.add_argument("--only", default=None,
                    help="comma separated profile names (default: all)")
    ap.add_argument("--jobs", type=int, default=4,
                    help="parallel profiles (default 4)")
    ap.add_argument("--dry-run", action="store_true",
                    help="compute and print, do not write taps/JSON")
    args = ap.parse_args()

    wanted = None
    if args.only:
        wanted = {s.strip() for s in args.only.split(",") if s.strip()}
        unknown = wanted - {p["name"] for p in PROFILES}
        if unknown:
            ap.error("unknown profile(s): %s" % ", ".join(sorted(unknown)))
    todo = [p for p in PROFILES if wanted is None or p["name"] in wanted]

    t0 = time.time()
    if args.jobs > 1 and len(todo) > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=min(args.jobs, len(todo))) as pool:
            results = pool.map(run_profile, todo)
    else:
        results = [run_profile(p) for p in todo]

    records = {}
    payload = {}
    for out in results:
        prof = out["prof"]
        print("[%s] f_oob=%dMHz extended=%s elapsed=%.1fs"
              % (out["name"], out["f_star"], out["extended"], out["elapsed_s"]))
        chosen = out["chosen"]
        if chosen is None:
            print("  !! no design produced")
            continue
        h32 = chosen["h32"]
        chk = self_check(h32)
        rec = build_record(prof, chosen, h32, chk)
        records[prof["name"]] = (rec, h32)
        payload[prof["name"]] = rec
        print("  T=%d Hpad=%d f_oob=%.0fMHz wo=%.1f wi=%.0f gd=%d"
              % (rec["taps"], rec["hpad"], rec["f_oob_mhz"],
                 rec["weight_oob"], rec["weight_img"], rec["grid_density"]))
        print("  dc=%.6f ripple=%.4fdB oob=%.2f/%.2f dB img=%.2f dB "
              "dontcare=%.2f dB spec_met=%s"
              % (rec["dc_sum"], rec["passband_ripple_db"], rec["oob_max_db"],
                 rec["oob_rms_db"], rec["img_max_db"], rec["dontcare_max_db"],
                 rec["spec_met"]))
        print("  self-check max|upfirdn - polyphase| = %.3e (len match=%s)"
              % (rec["self_check_upfirdn_max_abs_diff"],
                 rec["self_check_length_match"]))

    if args.dry_run:
        print("dry-run: nothing written (%.1fs)" % (time.time() - t0))
        return 0

    if len(records) != len(todo):
        sys.stderr.write("refusing to write: not all requested profiles produced\n")
        return 1

    os.makedirs(args.out_dir, exist_ok=True)
    for name, (rec, h32) in records.items():
        path = os.path.join(args.out_dir, rec["taps_file"])
        np.ascontiguousarray(h32, dtype="<f4").tofile(path)
        # prove the on-disk file round-trips to the same float32 array
        on_disk = np.fromfile(path, dtype="<f4")
        if on_disk.size != h32.size or not np.array_equal(on_disk, h32):
            raise RuntimeError("round-trip mismatch for %s" % path)

    if wanted is None:
        # full run -> single combined JSON
        design = dict(
            interpolation=L,
            decimation=M,
            input_rate_hz=R_IN,
            output_rate_hz=R_OUT,
            virtual_rate_hz=FV,
            design_method=("3-band don't-care Parks-McClellan "
                           "(scipy.signal.remez); bands passband/oob-stop/"
                           "image-stop, gap (368.64, 487.28) MHz left "
                           "unconstrained and only checked; odd symmetric "
                           "linear-phase FIR, DC-normalised to sum=65, "
                           "float32 quantised; all metrics measured after "
                           "float32 quantisation relative to the passband gain"),
            bands=dict(
                passband_hz=[0.0, F_PASS],
                band1_oob_stop_hz=dict(lo_hz="per-profile f_oob",
                                       hi_hz=BAND1_HI),
                dontcare_hz=[BAND1_HI, BAND2_LO],
                band2_image_stop_hz=[BAND2_LO, FV_HALF],
            ),
            dc_target=DC_TARGET,
            extension_limit_hz=EXTENSION_HZ,
            profiles=[payload[p["name"]] for p in PROFILES],
            generation=dict(
                script="tools/dsp/design_resampler_65_48_ch250.py",
                scipy=scipy.__version__,
                numpy=np.__version__,
                remez_maxiter=REMEZ_MAXITER,
                search_nfft=SEARCH_NFFT,
                final_nfft=FINAL_NFFT,
                self_check_seed=SELFCHECK_SEED,
                self_check_len=SELFCHECK_LEN,
            ),
        )
        with open(os.path.join(args.out_dir, "design_ch250.json"), "w") as f:
            json.dump(design, f, indent=2, allow_nan=False)
            f.write("\n")
        print("wrote %s + %d taps files (%.1fs)"
              % (os.path.join(args.out_dir, "design_ch250.json"),
                 len(records), time.time() - t0))
    else:
        # partial run -> keep partial JSON separate, never clobber combined
        part = os.path.join(args.out_dir, "design_ch250_partial_%s.json"
                            % "_".join(sorted(wanted)))
        with open(part, "w") as f:
            json.dump(dict(profiles=[payload[p["name"]] for p in PROFILES
                                     if p["name"] in payload]),
                      f, indent=2, allow_nan=False)
            f.write("\n")
        print("wrote %s + %d taps files (%.1fs)"
              % (part, len(records), time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
