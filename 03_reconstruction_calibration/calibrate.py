"""ARC card calibration v3 from the cell-averaged, PDE-sampled full sim (T06b, arc_table_t06b.csv).

Per theta bin (whole cells) and radiator, two-step physical fit over species and momenta:
  per photon:  sigma1^2 * n/(n-1) = s1^2 + (k * 0.3 B Lsin / p)^2          (aerogel: k = 0, bending < 0.1 mrad)
  per track:   sigma_track^2 = (s1^2 + (k 0.3 B Lsin / p)^2) / N_pe + s_c^2  -> s_c
  yield:       N(beta ~ 1) = mean photons of the found rings, pi+ p >= 20 GeV; efficiency = fraction of tracks with a ring (>= 3 photons)
Smooth forms in the cell off-axis variable x = cot(theta) (barrel) or tan(theta) (endcap), fitted per region; the gap bin
(theta 39.3-44.2 deg, between the last barrel cell and the outer edge of the endcap cells) keeps its own measured values.
Symmetry: the response is even in cot(theta) (theta <-> 180 - theta) and in tan(theta) (rotation about the beam axis), so every
quantity is an EVEN polynomial in x: s1^2, s_c^2, Y, E = c0 + c1 x^2 + c2 x^4 + ...; the lowest degree whose chi2/ndf is within 1 of
the best (or below 1.5) is kept, at most 3 terms (x^4; 4 terms, x^6, for the yield); x is held at the outermost bin centre (no extrapolation).
Y is the yield relative to the card's path scaling N_sat (L/L_ref) with L = t/sin (barrel) or t/|cos| (endcap) as in ARCPID,
N_sat fixed by the barrel fit at x = 0. Output calib_t06b.json, plots/calib_t06b.{pdf,png} (LaTeX, theta, 400-point curves)."""
import csv, json, math, os
import numpy as np
from scipy.optimize import least_squares

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.environ["PATH"] = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:" + os.environ["PATH"]
plt.rcParams.update({"lines.linewidth": 1.8, "lines.markersize": 5, "text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Roman"], "font.size": 12, "legend.fontsize": 9})

D = os.environ.get("T06B_DIR", "/eos/user/h/hfatehi/alfa/work/T06b")  # T06c: C4F10 variant
B_T = 3.0
CARD = {"gas": dict(t=0.169, lref=0.15, pix=1.25, chrom=float(os.environ.get("T06B_GAS_CHROM", 0.81))), "aerogel": dict(t=0.0098, lref=0.01, pix=1.23, chrom=2.68)}
COT_SPLIT = 2.2 / 1.9  # ARCPID region rule: barrel if the straight line reaches r = 1.9 m at |z| <= 2.2 m

rows = list(csv.reader(l for l in open(f"{D}/arc_table_t06b.csv") if not l.startswith("#")))
R = [dict(zip(rows[0], r)) for r in rows[1:]]
BINS = sorted({(float(r["theta_lo"]), float(r["theta_hi"]), r["region"]) for r in R})


def fnum(v):
    try:
        return float(v)
    except ValueError:
        return float("nan")


def card_path(th_deg, t):
    th = math.radians(th_deg)
    return t / math.sin(th) if 1 / math.tan(th) <= COT_SPLIT else t / abs(math.cos(th))


def per_bin(rad):
    out = []
    for lo, hi, reg in BINS:
        S = [r for r in R if float(r["theta_lo"]) == lo and r["radiator"] == rad]
        thm = float(np.mean([fnum(r["theta_mean"]) for r in S]))
        x = 1 / math.tan(math.radians(thm)) if reg == "barrel" else math.tan(math.radians(thm))
        Ls = CARD[rad]["t"] if reg != "endcap" else CARD[rad]["t"] * math.tan(math.radians(thm))
        nmin = 10 if reg == "gap" else 30  # the gap has few rings (gas ring efficiency ~0.12): accept smaller samples there
        F = [r for r in S if np.isfinite(fnum(r["sigma_thetac_mrad"])) and np.isfinite(fnum(r["sigma1_mrad"])) and int(r["n_tracks"]) >= nmin and fnum(r["mean_npe_found"]) >= 4]
        res = dict(lo=lo, hi=hi, region=reg, theta=thm, x=x)
        if len(F) >= (2 if reg == "gap" else 3):
            p = np.array([fnum(r["p_GeV"]) for r in F]); st = np.array([fnum(r["sigma_thetac_mrad"]) for r in F])
            n = np.array([fnum(r["mean_npe_found"]) for r in F]); nt = np.array([int(r["n_tracks"]) for r in F])
            s1m = np.array([fnum(r["sigma1_mrad"]) for r in F]) * np.sqrt(n / (n - 1))
            b2 = (0.3 * B_T * Ls / p * 1e3) ** 2
            w1 = (nt * n) / s1m ** 4
            if rad == "gas" and reg == "gap":
                # too few points for s1 and k together: k fixed to the barrel/endcap value (0.155), s1 from the weighted mean
                k, ke = 0.155, 0.0
                c0 = float(np.sum(w1 * (s1m ** 2 - k * k * b2)) / np.sum(w1)); s1 = math.sqrt(max(c0, 0.01))
                chi2 = float(np.sum(w1 * (s1m ** 2 - k * k * b2 - c0) ** 2) / max(len(p) - 1, 1)); s1e = math.sqrt(max(chi2, 1) / np.sum(w1)) / (2 * s1)
            elif rad == "gas":
                A = np.vstack([np.ones_like(b2), b2]).T
                c = np.linalg.lstsq(A * np.sqrt(w1)[:, None], s1m ** 2 * np.sqrt(w1), rcond=None)[0]
                chi2 = float(np.sum((s1m ** 2 - A @ c) ** 2 * w1) / max(len(p) - 2, 1))
                cov = np.linalg.inv((A * w1[:, None]).T @ A) * max(chi2, 1)
                s1 = math.sqrt(max(c[0], 0.01)); k = math.sqrt(max(c[1], 0)); s1e = math.sqrt(cov[0, 0]) / (2 * s1); ke = math.sqrt(cov[1, 1]) / (2 * max(k, 1e-3))
            else:
                c0 = float(np.sum(w1 * s1m ** 2) / np.sum(w1)); s1 = math.sqrt(c0); k = ke = 0.0
                chi2 = float(np.sum(w1 * (s1m ** 2 - c0) ** 2) / max(len(p) - 1, 1)); s1e = math.sqrt(max(chi2, 1) / np.sum(w1)) / (2 * s1)
            sc2 = st ** 2 - (s1 ** 2 + k * k * b2) / n; wsc = nt / st ** 4
            if reg == "gap":
                wsc = nt.astype(float)  # 6-33 rings per row: inverse-variance weights from the noisy sigma_track bias the mean low
            sc2m = float(np.sum(wsc * sc2) / np.sum(wsc))
            sc2e = math.sqrt(1 / np.sum(wsc)) * max(1.0, math.sqrt(np.sum(wsc * (sc2 - sc2m) ** 2) / max(len(p) - 1, 1)))
            sc = math.sqrt(max(sc2m, 0)); sce = sc2e / (2 * max(sc, 0.1))
            res.update(s1=s1, s1e=s1e, k=k, ke=ke, sc=sc, sce=sce, chi2_s1=chi2, npts=len(F))
        Y = [r for r in S if r["species"] == "pi+" and fnum(r["p_GeV"]) >= 20 and int(r["n_events"]) > 0]
        nev = np.array([int(r["n_events"]) for r in Y]); E = np.array([fnum(r["eff_ring_found"]) for r in Y])
        Nf = np.array([fnum(r["mean_npe_found"]) for r in Y]); nt = np.array([int(r["n_tracks"]) for r in Y])
        Eb = float(np.sum(E * nev) / nev.sum()); res.update(E=Eb, Ee=math.sqrt(max(Eb * (1 - Eb), 0.005) / nev.sum()))
        ok = nt > 0
        if ok.any():
            Nb = float(np.sum(Nf[ok] * nt[ok]) / nt[ok].sum()); Nse = math.sqrt(Nb / nt[ok].sum())  # Poisson-like spread per ring
            res.update(N=Nb, Ne=Nse)
        out.append(res)
    return out


def fit(points, key, ekey, kind):
    """even polynomial in x; for the resolutions (kind 'quad') the polynomial is fitted to q^2 and the model returns sqrt."""
    P = [q for q in points if key in q and np.isfinite(q[key])]
    x = np.array([q["x"] for q in P]); y = np.array([q[key] for q in P]); e = np.array([max(q[ekey], 1e-3) for q in P])
    if kind == "quad":
        yy, ee = y * y, 2 * y * e
    else:
        yy, ee = y, e
    cands = []
    for nt in range(1, min(4 if kind == "yield" else 3, max(len(y) - 2, 1)) + 1):
        A = np.vstack([x ** (2 * i) for i in range(nt)]).T
        c, *_ = np.linalg.lstsq(A / ee[:, None], yy / ee, rcond=None)
        ndf = max(len(y) - nt, 1); chi2 = float(np.sum(((A @ c - yy) / ee) ** 2) / ndf)
        cov = np.linalg.pinv((A / ee[:, None]).T @ (A / ee[:, None])) * max(chi2, 1)
        cur = dict(p=[float(v) for v in c], perr=[float(v) for v in np.sqrt(np.abs(np.diag(cov)))], chi2ndf=chi2, ndf=ndf, kind=kind,
                   xmin=float(x.min()), xmax=float(x.max()))
        cands.append(cur)
    # lowest degree that is within 1 of the best chi2/ndf (or already below 1.5)
    top = min(c["chi2ndf"] for c in cands)
    best = next(c for c in cands if c["chi2ndf"] <= max(top + 1.0, 1.5))
    poly = lambda c, x: sum(ci * x ** (2 * i) for i, ci in enumerate(c))
    model = (lambda c, x: np.sqrt(np.maximum(poly(c, x), 0))) if kind == "quad" else (lambda c, x: np.clip(poly(c, x), 0, None) if kind == "yield" else np.clip(poly(c, x), 0, 1))
    return best, model


out = {}
KEYS = (("s1", "s1e", "quad", r"per-photon $s_1$ [mrad] (no bending)"), ("sc", "sce", "quad", r"per-track correlated $\sigma_c$ [mrad]"),
        ("Y", "Ye", "yield", r"yield scale $N / (N_\mathrm{sat} L/L_\mathrm{ref})$"), ("E", "Ee", "eff", r"ring-finding efficiency"))
fig, axes = plt.subplots(2, 4, figsize=(20, 10.5), gridspec_kw=dict(hspace=0.62))
for j, rad in enumerate(("gas", "aerogel")):
    pts = per_bin(rad)
    bar = [q for q in pts if q["region"] == "barrel"]
    # N_sat from the barrel yield extrapolated to x = 0: fit N = N0 * (1 + x^2)^(-g/2) * card path ratio ... done via Y with N_sat = N(x=0)
    N0 = float(np.interp(0.0, [q["x"] for q in sorted(bar, key=lambda q: q["x"])], [q["N"] for q in sorted(bar, key=lambda q: q["x"])]))
    nsat = N0 * CARD[rad]["lref"] / CARD[rad]["t"]
    for q in pts:
        if "N" in q:
            card = nsat * card_path(q["theta"], CARD[rad]["t"]) / CARD[rad]["lref"]
            q["Y"] = q["N"] / card; q["Ye"] = q["Ne"] / card
    res = dict(Nsat=nsat, points=pts, fits={})
    thf = np.linspace(12, 90, 400)
    for i, (k, ek, kind, lab) in enumerate(KEYS):
        ax = axes[j, i]
        for reg, col in (("barrel", "#1f5f8a"), ("endcap", "#2ca02c")):
            P = [q for q in pts if q["region"] == reg]
            fr, model = fit(P, k, ek, kind)
            res["fits"][f"{k}_{reg}"] = fr
            ax.errorbar([q["theta"] for q in P if k in q], [q[k] for q in P if k in q], xerr=[[q["theta"] - q["lo"] for q in P if k in q], [q["hi"] - q["theta"] for q in P if k in q]],
                        yerr=[q[ek] for q in P if k in q], fmt="o", color=col, capsize=0, label=f"full sim {reg} (cell-averaged)")
            tr = thf[(thf >= P[0]["lo"]) & (thf <= P[-1]["hi"])] if P else thf[:0]
            xr = 1 / np.tan(np.radians(tr)) if reg == "barrel" else np.tan(np.radians(tr))
            ax.plot(tr, model(np.array(fr["p"]), np.clip(xr, fr["xmin"], fr["xmax"])), "--", color=col, label=f"fit, $\\chi^2/\\mathrm{{ndf}} = {fr['chi2ndf']:.1f}$")
        G = [q for q in pts if q["region"] == "gap" and k in q]
        if G:
            ax.errorbar([q["theta"] for q in G], [q[k] for q in G], xerr=[[q["theta"] - q["lo"] for q in G], [q["hi"] - q["theta"] for q in G]], yerr=[q[ek] for q in G],
                        fmt="s", mfc="white", color="#d62728", label="gap (own value)")
        ax.set_xlim(10, 92); ax.set_xlabel(r"$\theta$ [deg]"); ax.set_ylabel(lab); ax.set_title(rad); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    out[rad] = res
fig.tight_layout(); os.makedirs(f"{D}/plots", exist_ok=True)
for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
    fig.savefig(f"{D}/plots/calib_t06b.{ext}", bbox_inches="tight", **kw)
json.dump(out, open(f"{D}/calib_t06b.json", "w"), indent=1, default=float)
for rad, r in out.items():
    print(f"{rad}: N_sat {r['Nsat']:.2f}")
    for k, f in r["fits"].items():
        print(f"  {k:10s} p = {[round(v, 3) for v in f['p']]}  +- {[round(v, 3) for v in f['perr']]}  chi2/ndf {f['chi2ndf']:.2f} (ndf {f['ndf']})")
    for q in r["points"]:
        print("   %5.1f-%5.1f %-6s th %5.1f x %.3f | s1 %s sc %s k %s | N %s Y %s E %.3f" % (q["lo"], q["hi"], q["region"], q["theta"], q["x"],
              f"{q['s1']:.2f}+-{q['s1e']:.2f}" if "s1" in q else "-", f"{q['sc']:.2f}+-{q['sce']:.2f}" if "sc" in q else "-", f"{q['k']:.3f}" if "k" in q else "-",
              f"{q['N']:.2f}" if "N" in q else "-", f"{q['Y']:.3f}" if "Y" in q else "-", q["E"]))
