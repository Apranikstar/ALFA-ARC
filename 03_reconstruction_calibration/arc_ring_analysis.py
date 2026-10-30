"""T06: per-track ARC observables for a DELPHES parameterisation, ALFA_o1_v01 (ARC xenon 2 bar).
Input: the .npz files written by the HTCondor jobs (condor/T06/t06_extract.py): primary + ARC photon hits.

Per event:
  - TOF from the primary's arrival at the ARC vessel (barrel r = 1910 mm, endcap |z| = 2010 mm); keep the
    single-reflection window 0.4-1.4 ns (ring light) in the cell with the most photons.
  - Photon weights: SiPM PDE (ARC_SiPM_QuantumEfficiency: FBK NUV-HD + 5% loss) from the photon energy.
  - Radiator assignment, truth-seeded: expected ring radii r = f * theta_out from the true species and momentum,
    with cos(theta) = 1/(n beta); aerogel light leaves the aerogel into the gas, sin(theta_out) = n_aero/n_gas sin(theta_aero).
    Ring centre = grid search maximising the PDE-weighted photons within +-2.5 mm of the expected radii;
    each photon goes to the nearer expected ring (within 6 mm); each ring is then fitted (centre fixed, radius =
    PDE-weighted median distance).
  - Cherenkov angle: gas theta_c = r_gas / f; aerogel theta_c = asin(n_gas/n_aero * sin(r_aero / f)).
    f is calibrated per theta from pi+ at 50 GeV (beta ~ 1) on the gas ring, because the mirror radii differ
    per cell (0.22-0.369 m in RadiatorCell_FinalOptimisation_o1_v01.xml) and a fixed theta always hits the
    same cell type.
  - Ring found: >= 3 PDE-weighted photons on that ring.
Per grid point and radiator: mean/rms detected photons, mean ring radius, mean theta_c, per-track
sigma(theta_c) (1.4826 * MAD of the per-event values), ring-finding efficiency.
"""
import glob
import math
import os
import re

import numpy as np
from scipy.optimize import least_squares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_MM_NS = 299.792458
R_BARREL, Z_ENDCAP = 1910.0, 2010.0
GATE_LO, GATE_HI = 0.4, 1.4  # single-reflection TOF window [ns] from the vessel entry (T02/T09)
QE_E = np.array([1.1, 1.3, 1.550, 1.771, 2.101, 2.465, 2.655, 2.743, 2.837, 2.938, 3.002, 3.084, 3.163, 3.378, 3.647, 3.850, 4.147, 4.320, 4.5, 4.7, 4.9])
QE_V = np.array([0.0, 0.0, 0.13, 0.26, 0.4, 0.51, 0.53, 0.56, 0.57, 0.59, 0.6, 0.6, 0.59, 0.52, 0.47, 0.46, 0.41, 0.2, 0.1, 0.0, 0.0])
MASS = {"pip": 0.13957, "kaonp": 0.493677, "proton": 0.938272}
LABEL = {"pip": "pi+", "kaonp": "K+", "proton": "p"}

# Detected-spectrum averaged indices (PDE x flat-in-energy Cherenkov spectrum, 1.55-4.7 eV)
_E = np.linspace(1.55, 4.7, 200)
_W = np.interp(_E, QE_E, QE_V)
XE_FACTOR = (2.0 / 1.01325) * (273.15 / 293.15)
_XE = 1 + XE_FACTOR * (0.00322869 / (46.301 - (_E / 1.239842) ** 2) + 0.00355393 / (59.578 - (_E / 1.239842) ** 2)
                      + 0.0606764 / (112.74 - (_E / 1.239842) ** 2))
_AERO_E = np.array([1.87855, 1.96673, 2.05490, 2.14308, 2.23126, 2.31943, 2.40761, 2.49579, 2.58396, 2.67214, 2.76032, 2.84849, 2.93667])
_AERO_N = np.array([1.01852, 1.01856, 1.01861, 1.01866, 1.01871, 1.01876, 1.01881, 1.01887, 1.01893, 1.01899, 1.01905, 1.01912, 1.01919])
_AE = np.interp(_E, _AERO_E, _AERO_N, right=_AERO_N[-1] + (_E[-1] - _AERO_E[-1]) * 7.5e-5)  # linear extrapolation above 2.94 eV
N_GAS = float(np.sum(_XE * _W) / np.sum(_W))
N_AERO = float(np.sum(_AE * _W) / np.sum(_W))


def theta_out(n, beta, aerogel):
    c = 1.0 / (n * beta)
    if c >= 1:
        return None
    th = math.acos(c)
    if aerogel:
        s = n / N_GAS * math.sin(th)
        return math.asin(s) if s < 1 else None
    return th


def load(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def per_event(d, mass, p_gev, f_mm):
    beta = p_gev / math.sqrt(p_gev ** 2 + mass ** 2)
    out = []
    th_g, th_a = theta_out(N_GAS, beta, False), theta_out(N_AERO, beta, True)
    exp = {"gas": f_mm * th_g if th_g else None, "aerogel": f_mm * th_a if th_a else None}
    for i in range(len(d["ev_pdg"])):
        px, py, pz = d["ev_px"][i], d["ev_py"][i], d["ev_pz"][i]
        P = math.sqrt(px * px + py * py + pz * pz)
        sin_t, cos_t = math.hypot(px, py) / P, abs(pz) / P
        L = min(R_BARREL / sin_t if sin_t > 0 else 1e12, Z_ENDCAP / cos_t if cos_t > 0 else 1e12)
        t0 = d["ev_t0"][i] + L / (beta * C_MM_NS)
        b, e = d["ev_hit_begin"][i], d["ev_hit_end"][i]
        res = {"gas": None, "aerogel": None, "outlier_frac": None, "gated_npe": 0.0, "centre": None}
        if e > b:
            cell = d["hit_cell"][b:e] + 1000000 * d["hit_barrel"][b:e]
            vals, cnt = np.unique(cell, return_counts=True)
            sel = (cell == vals[np.argmax(cnt)])
            tof = d["hit_t_ns"][b:e] - t0
            sel &= (tof >= GATE_LO) & (tof < GATE_HI)
            x, y = d["hit_x_mm"][b:e][sel].astype(float), d["hit_y_mm"][b:e][sel].astype(float)
            w = np.interp(d["hit_e_ev"][b:e][sel], QE_E, QE_V, left=0.0, right=0.0)
            if "hit_isnoise" in d:  # injected dark counts are detected counts: weight 1
                w = np.where(d["hit_isnoise"][b:e][sel] > 0, 1.0, w)
            rings = [(k, r) for k, r in exp.items() if r]
            if len(x) >= 3 and rings:
                best = None
                cx0, cy0 = np.average(x, weights=w + 1e-9), np.average(y, weights=w + 1e-9)
                for step, span in ((1.0, 20.0), (0.25, 1.5)):
                    g = np.arange(-span, span + 1e-9, step)
                    for dx in g:
                        for dy in g:
                            dist = np.hypot(x - cx0 - dx, y - cy0 - dy)
                            score = sum(w[np.abs(dist - r) < 2.5].sum() for _, r in rings)
                            if best is None or score > best[0]:
                                best = (score, cx0 + dx, cy0 + dy)
                    cx0, cy0 = best[1], best[2]
                dist = np.hypot(x - cx0, y - cy0)
                dev = np.array([np.abs(dist - r) for _, r in rings])
                near = np.argmin(dev, axis=0)
                on_ring = dev.min(axis=0) < 6.0
                res["gated_npe"] = float(w.sum())
                res["outlier_frac"] = float(w[~on_ring].sum() / max(w.sum(), 1e-9))
                # joint fit: one common centre + one radius per ring, PDE-weighted, robust (soft L1, 1 mm scale)
                masks = [(near == j) & (dev[j] < 6.0) for j in range(len(rings))]
                used = [j for j, m in enumerate(masks) if w[m].sum() > 0]
                if used:
                    sel_all = np.zeros_like(near, dtype=bool)
                    for j in used:
                        sel_all |= masks[j]
                    idx = np.array([next(j for j in used if masks[j][i]) for i in np.flatnonzero(sel_all)])
                    xs_, ys_, ws_ = x[sel_all], y[sel_all], w[sel_all]
                    p0 = [cx0, cy0] + [rings[j][1] for j in used]

                    def resid(p):
                        dd = np.hypot(xs_ - p[0], ys_ - p[1])
                        return np.sqrt(ws_) * (dd - np.array([p[2 + used.index(j)] for j in idx]))
                    if len(xs_) >= len(p0):
                        fit = least_squares(resid, p0, loss="soft_l1", f_scale=1.0)
                        cx0, cy0 = fit.x[0], fit.x[1]
                        radii = {j: fit.x[2 + k] for k, j in enumerate(used)}
                    else:
                        radii = {j: rings[j][1] for j in used}
                    for j in used:
                        npe = float(w[masks[j]].sum())
                        res[rings[j][0]] = (npe, float(radii[j]))
                    res["centre"] = (float(cx0), float(cy0))
        out.append(res)
    return out, exp


def robust_sigma(v):
    v = np.asarray(v)
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) > 2 else float("nan")


def main():
    files = sorted(glob.glob("npz/*.npz"))
    pts = {}
    for f in files:
        m = re.match(r"(pip|kaonp|proton)_p([\d.]+)_th(\d+)\.npz", os.path.basename(f))
        if m:
            pts[(m.group(1), float(m.group(2)), int(m.group(3)))] = f
    thetas = sorted({k[2] for k in pts})

    # Calibration per theta and radiator from pi+ 50 GeV (beta ~ 1, true angle known): the ring radius depends on
    # where the ring lands in the cell (mirror aberration vs impact point; the mirror radii also differ per cell).
    # Model q = r_fit / r_expected(F0) as a quadratic in the fitted ring centre (cx, cy); a track at any momentum
    # is corrected as r / q(cx, cy), so the angle is r_corr / F0 (gas) or asin(n_gas/n_aero sin(r_corr / F0)).
    F0 = 184.5
    def design(c):
        c = np.atleast_2d(c)
        return np.c_[np.ones(len(c)), c[:, 0], c[:, 1], c[:, 0] ** 2, c[:, 1] ** 2, c[:, 0] * c[:, 1]]
    qcal = {}
    for th in thetas:
        key = ("pip", 50.0, th)
        if key not in pts:
            continue
        ev, exp = per_event(load(pts[key]), MASS["pip"], 50.0, F0)
        for rad in ("gas", "aerogel"):
            sel = [(e["centre"], e[rad][1] / exp[rad]) for e in ev if e[rad] and e[rad][0] >= 3 and e["centre"] and exp[rad]]
            if len(sel) >= 12:
                C = np.array([c for c, _ in sel]); q = np.array([v for _, v in sel])
                coef, *_ = np.linalg.lstsq(design(C), q, rcond=None)
                qcal[(th, rad)] = coef
    f_cal = {th: F0 / float(qcal[(th, "gas")][0]) for th in thetas if (th, "gas") in qcal}  # report only

    rows = []
    for (sp, p, th), f in sorted(pts.items()):
        fmm = F0
        ev, exp = per_event(load(f), MASS[sp], p, fmm)
        for rad in ("aerogel", "gas"):
            okev = [e for e in ev if e[rad] and e[rad][0] >= 3]
            found = [e[rad] for e in okev]
            npe_all = [e[rad][0] if e[rad] else 0.0 for e in ev]
            if found:
                r = np.array([x[1] for x in found])
                if (th, rad) in qcal and all(e["centre"] for e in okev):
                    q = design(np.array([e["centre"] for e in okev])) @ qcal[(th, rad)]
                    # the calibration is a smooth fit on beta~1 pions: outside its validity (e.g. near threshold, where the
                    # centre is poorly constrained) it can explode, so keep it within +-15% and otherwise leave r as is
                    r = np.where((q > 0.85) & (q < 1.15), r / np.where(q > 0, q, 1.0), r)
                thc = r / fmm if rad == "gas" else np.arcsin(np.clip(N_GAS / N_AERO * np.sin(r / fmm), -1, 1))
                mean_thc, sig_thc, mean_r = float(np.mean(thc) * 1e3), robust_sigma(thc * 1e3), float(np.mean(r))
            else:
                mean_thc = sig_thc = mean_r = float("nan")
            rows.append(dict(species=LABEL[sp], p_GeV=p, theta_deg=th, radiator=rad, n_events=len(ev),
                             mean_npe=float(np.mean(npe_all)), rms_npe=float(np.std(npe_all)),
                             mean_ring_mm=mean_r, mean_thetac_mrad=mean_thc, sigma_thetac_mrad=sig_thc,
                             eff_ring_found=len(found) / max(len(ev), 1), f_mm=fmm,
                             expected_ring_mm=exp[rad] if exp[rad] else float("nan")))

    cols = ["species", "p_GeV", "theta_deg", "radiator", "n_events", "mean_npe", "rms_npe", "mean_thetac_mrad",
            "sigma_thetac_mrad", "eff_ring_found", "mean_ring_mm", "expected_ring_mm", "f_mm"]
    with open("arc_table_xenon.csv", "w") as fo:
        fo.write("# ALFA_o1_v01 ARC, xenon 2 bar 293.15 K; n_gas(detected avg) = %.6f, n_aerogel(detected avg) = %.5f; SiPM PDE applied\n" % (N_GAS, N_AERO))
        fo.write(",".join(cols) + "\n")
        for r in rows:
            fo.write(",".join(("%.4g" % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    # plots: N_pe and sigma vs p per radiator (theta = 90), and K/pi separation vs p per radiator and theta
    import collections
    tab = collections.defaultdict(dict)
    for r in rows:
        tab[(r["species"], r["theta_deg"], r["radiator"])][r["p_GeV"]] = r
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    colours = {"pi+": "#1f6f5c", "K+": "#c2571a", "p": "#5a4fb0"}
    th_show = 90 if 90 in thetas else thetas[-1]
    for j, rad in enumerate(("aerogel", "gas")):
        for sp in ("pi+", "K+", "p"):
            d = tab[(sp, th_show, rad)]
            ps = sorted(d)
            ax[0, j].plot(ps, [d[q]["mean_npe"] for q in ps], "o-", color=colours[sp], label=sp)
            ax[1, j].plot(ps, [d[q]["sigma_thetac_mrad"] for q in ps], "o-", color=colours[sp], label=sp)
        ax[0, j].set_title(f"{rad}: detected photons per track (theta = {th_show} deg)")
        ax[1, j].set_title(f"{rad}: per-track sigma(theta_c) (theta = {th_show} deg)")
        for a in (ax[0, j], ax[1, j]):
            a.set_xscale("log"); a.set_xlabel("p [GeV]"); a.legend(fontsize=9); a.grid(alpha=0.3)
        ax[0, j].set_ylabel("photons"); ax[1, j].set_ylabel("mrad")
    fig.tight_layout(); fig.savefig("arc_npe_sigma_vs_p.png", dpi=110)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    summary_lines = []
    for j, rad in enumerate(("aerogel", "gas")):
        for th in thetas:
            dpi, dk = tab[("pi+", th, rad)], tab[("K+", th, rad)]
            ps = sorted(set(dpi) & set(dk))
            sep = []
            for q in ps:
                a, b = dpi[q], dk[q]
                if a["eff_ring_found"] > 0.5 and b["eff_ring_found"] > 0.5:
                    s = abs(a["mean_thetac_mrad"] - b["mean_thetac_mrad"]) / math.sqrt(a["sigma_thetac_mrad"] ** 2 + b["sigma_thetac_mrad"] ** 2)
                else:
                    s = float("nan")
                sep.append(s)
            ax[j].plot(ps, sep, "o-", label=f"theta {th}")
            above = [q for q, s in zip(ps, sep) if s == s and s >= 3]
            summary_lines.append(f"{rad} theta {th}: 3-sigma K/pi up to {max(above) if above else 'none'} GeV (grid points)")
        ax[j].axhline(3, color="k", ls="--", lw=1)
        ax[j].set_xscale("log"); ax[j].set_yscale("log"); ax[j].set_xlabel("p [GeV]"); ax[j].set_ylabel("K/pi separation [sigma]")
        ax[j].set_title(f"{rad}: K/pi separation (per track)"); ax[j].legend(fontsize=8, ncol=2); ax[j].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("arc_kpi_separation.png", dpi=110)

    print("n_gas %.6f n_aero %.5f; f per theta (mm): %s" % (N_GAS, N_AERO, {k: round(v, 1) for k, v in f_cal.items()}))
    for r in rows:
        if r["species"] == "pi+" and r["p_GeV"] == 10 and r["theta_deg"] == 90:
            print("pi+ 10 GeV 90 deg:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print("\n".join(summary_lines))
    print("points", len(pts), "rows", len(rows))


if __name__ == "__main__":
    main()
