"""[campaign_v2 copy (T17): v4 cards (campaign_v2/cards), full sim = v2 tables (campaign_v2/calib/<gas>), outputs in
campaign_v2/delphes_validation (runs/, plots/<variant>/, validation_<variant>.json), PNG 600 dpi]
T11 validation of the ALFA Delphes cards (ARCPID) against the card formula and the k4geo full simulation
(T06b: theta-uniform, cell-averaged theta bins, PDE sampled).

Usage: python3 t11_validate.py <card-variant: Xe2bar|C4F10>
Inputs: runs/<variant>_{pip,kp,p}_<chunk>.root (1M single-particle gun events per species, log-uniform p 0.1-60 GeV, theta 8.5-171.5 deg),
        runs/<variant>_zqq.root; the card cards/delphes_card_ALFA_<variant>.tcl (ARCPID parameters are read from it).
Plots (PDF + 300 dpi PNG, board style rules) in plots/<variant>/, numbers in validation_<variant>.json.
Separation per p bin: N_sigma = |theta_a(p) - theta_b(p)| / sqrt(s_a^2 + s_b^2), with s = rms of (reco - true angle) of each species
in the bin (the measured per-track spread), and theta_h(p) at the bin centre."""
import glob
import json
import math
import os
import re
import sys

import numpy as np
import ROOT
from scipy.stats import beta as betadist

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt

# LaTeX text rendering (TeX Live on cvmfs; not in the key4hep view)
os.environ["PATH"] = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:" + os.environ["PATH"]
plt.rcParams.update({"lines.linewidth": 1.8, "lines.markersize": 5, "lines.antialiased": True,
                     "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "text.latex.preamble": r"\usepackage{amsmath}", "font.size": 12, "legend.fontsize": 9})
TEXSP = {"pi": r"\pi", "K": r"K", "p": r"p"}
TEXV = {"Xe2bar": r"Xe 2\,bar", "C4F10": r"C$_4$F$_{10}$"}[sys.argv[1]]


def pair_tex(a, b):
    return f"${TEXSP[a]}/{TEXSP[b]}$"
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.gSystem.Load("/afs/cern.ch/work/h/hfatehi/alfa/delphes/libDelphes.so")  # our ARCPID build, not the key4hep stack's Delphes
ROOT.EnableImplicitMT(16)

V = sys.argv[1]
W = "/eos/user/h/hfatehi/alfa/campaign_v2/delphes_validation"
D = "/afs/cern.ch/work/h/hfatehi/alfa/delphes"
OUT = f"{W}/plots/{V}"
os.makedirs(OUT, exist_ok=True)
MASS = {"pi": 0.13957039, "K": 0.493677, "p": 0.93827209}
PDG = {"pi": 211, "K": 321, "p": 2212}
COL = {"pi": "#1f77b4", "K": "#d62728", "p": "#2ca02c", "pi/K": "#1f5f8a", "K/p": "#b5651d"}
TH_SEP = (82.6, 97.4)  # the central barrel z-cell (|z| < 25.7 cm at r = 2 m), i.e. the full-sim bin 82.6-90 deg and its mirror

# ---- card parameters
card = open(f"/eos/user/h/hfatehi/alfa/campaign_v2/cards/delphes_card_ALFA_{V}.tcl").read()
blk = card[card.index("module ARCPID ARCPID"):]
blk = blk[:blk.index("\n}\n")]
P = {k: float(v) for k, v in re.findall(r"^\s*set (\w+) ([-+0-9.eE]+)\s*$", blk, re.M)}
B = float(re.search(r"^set B ([0-9.]+)", card, re.M).group(1))
RAD = {r: dict(n=P[f"{pre}RefractiveIndex"], dnde=P[f"{pre}DnDe"], sE=P[f"{pre}SigmaE"], s1=P[f"{pre}SigmaSingle"],
              nsat=P[f"{pre}Nsat"], lref=P[f"{pre}Lref"], t=P[f"{pre}Thickness"], pre=pre)
       for r, pre in (("gas", "Gas"), ("aerogel", "Aerogel"))}
KB, MINPE = P["BendingCoefficient"], P["MinPe"]


def card_angle(r, p, m):
    c = 1.0 / (RAD[r]["n"] * p / np.sqrt(p * p + m * m))
    return np.where(c < 1, np.arccos(np.minimum(c, 1)), np.nan)


# the card's own TFormulas (the module evaluates them per track in eta), evaluated with ROOT for the analytic expectation
FORM = {}
for name, body in re.findall(r"^\s*set (\w+Formula) \{(.*?)\n\s*\}", blk, re.M | re.S):
    f = ROOT.TFormula(f"f_{name}", " ".join(body.split()).replace("eta", "x"))
    FORM[name] = f


def card_eval(name, eta, default):
    return FORM[name].Eval(abs(eta)) if name in FORM else default


def card_sigma(r, p, m, theta=math.pi / 2):
    """per-track sigma and N_pe of the module at polar angle theta (barrel path), with the card formulas."""
    R = RAD[r]; th = card_angle(r, p, m); pre = R["pre"]
    eta = -math.log(math.tan(theta / 2))
    L = R["t"] / math.sin(theta)
    sin2max = 1 - 1 / R["n"] ** 2
    npe = card_eval(f"{pre}YieldScaleFormula", eta, 1.0) * R["nsat"] * (L / R["lref"]) * np.sin(th) ** 2 / sin2max
    s1 = card_eval(f"{pre}SigmaSingleFormula", eta, R["s1"]); sc = card_eval(f"{pre}SigmaCorrFormula", eta, 0.0)
    chrom = R["dnde"] * R["sE"] / (R["n"] * np.tan(th))
    bend = KB * 0.3 * B * L * math.sin(theta) / p
    s = np.sqrt((s1 ** 2 + chrom ** 2 + bend ** 2) / npe + sc ** 2)
    return np.where(npe >= MINPE, s, np.nan), npe


# ---- data
def load(sp):
    """Flat per-track arrays (uproot; the RDataFrame -> numpy conversion of 1M RVecs is too slow)."""
    import uproot, awkward as ak
    files = sorted(glob.glob(f"{W}/runs/{V}_{sp}_*.root"))
    cols = ["P", "PT", "Eta", "PID", "CherenkovAngleGas", "CherenkovAngleAerogel", "CherenkovRecoAngleGas", "CherenkovRecoAngleAerogel",
            "SigmaThetaCGas", "SigmaThetaCAerogel", "NpeGas", "NpeAerogel", "ARCpid", "ARCprobPi", "ARCprobK", "ARCprobP", "ARCmass",
            "ErrorD0", "ErrorPT", "D0"]
    parts = {c: [] for c in cols}; ntrk = []
    for f in files:
        t = uproot.open(f)["Delphes"]
        arr = t.arrays([f"Track.{c}" for c in cols] + ["Track_size"], library="ak")
        for c in cols:
            parts[c].append(ak.to_numpy(ak.flatten(arr[f"Track.{c}"])).astype(float))
        ntrk.append(ak.to_numpy(arr["Track_size"]).astype(float))
    out = {c: np.concatenate(parts[c]) for c in cols}
    out["ntrk_per_event"] = np.concatenate(ntrk)
    n_ev = len(out["ntrk_per_event"])
    # the gun particle only (drop secondaries from conversions etc.; the gun is the only primary)
    sel = np.abs(out["PID"]) == {"pip": 211, "kp": 321, "p": 2212}[sp]
    out = {k: (v[sel] if k != "ntrk_per_event" else v) for k, v in out.items()}
    out["n_events"] = n_ev; out["files"] = len(files)
    return out


d = {"pi": load("pip"), "K": load("kp"), "p": load("p")}
for s_ in d:  # polar angle [deg]; the plots use theta, not eta
    d[s_]["Theta"] = np.degrees(2 * np.arctan(np.exp(-d[s_]["Eta"])))
TH_ALL = r"$8.5^\circ < \theta < 171.5^\circ$"
counts = {s: dict(events=int(d[s]["n_events"]), tracks=int(len(d[s]["P"])), files=d[s]["files"]) for s in d}
print(json.dumps(counts))

edges = np.logspace(np.log10(0.5), np.log10(60), 61)
cent = np.sqrt(edges[:-1] * edges[1:])
pf = np.logspace(np.log10(0.5), np.log10(60), 400)


def cp_band(k, n, cl=0.6827):
    a = (1 - cl) / 2
    lo = np.where(k > 0, betadist.ppf(a, k, n - k + 1), 0.0)
    hi = np.where(k < n, betadist.ppf(1 - a, k + 1, n - k), 1.0)
    return lo, hi


def save(fig, name):
    fig.tight_layout()
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
        fig.savefig(f"{OUT}/{name}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)


# ---- full sim T06b (Xe 2 bar): per (species, p, theta bin, radiator)
import csv
FULLSIM = {"Xe2bar": f"{os.path.dirname(W)}/calib/Xe2bar/arc_table_t06b.csv", "C4F10": f"{os.path.dirname(W)}/calib/C4F10/arc_table_t06b.csv"}[V]
HAVE_FS = os.path.exists(FULLSIM)  # each card against the full sim of its own gas
FS_TAG = "T06b" if V == "Xe2bar" else "T06c"
rows = list(csv.reader(l for l in open(FULLSIM) if not l.startswith("#")))
R6 = [dict(zip(rows[0], r)) for r in rows[1:]]
SPN = {"pi+": "pi", "K+": "K", "p": "p"}
def fsim(rad, lo):
    """{(species, p): (mean, sigma, n_tracks)} in the full-sim bin starting at theta_lo = lo"""
    out = {}
    for r in R6:
        if r["radiator"] == rad and abs(float(r["theta_lo"]) - lo) < 0.05 and r["sigma_thetac_mrad"] not in ("nan", "") and int(r["n_tracks"]) >= 30:
            out[(SPN[r["species"]], float(r["p_GeV"]))] = (float(r["mean_thetac_mrad"]), float(r["sigma_thetac_mrad"]), int(r["n_tracks"]))
    return out
result = {"variant": V, "card": {k: v for k, v in P.items()}, "B_T": B, "counts": counts, "reach_3sigma_GeV": {}, "fullsim_T06b_reach_central_cell_GeV": {}}


def crossing(x, y):
    """last 3-sigma crossing, log-log interpolation"""
    x, y = np.asarray(x, float), np.asarray(y, float); g = np.isfinite(y) & (y > 0)
    x, y = x[g], y[g]
    above = np.where(y >= 3)[0]
    if not above.size or above[-1] == len(y) - 1:
        return float("nan")
    i = above[-1]
    return float(math.exp(np.interp(math.log(3), [math.log(y[i + 1]), math.log(y[i])], [math.log(x[i + 1]), math.log(x[i])])))


# ---- 1. separation power, barrel slice TH_SEP
fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
PMIN = {"gas": 10.0, "aerogel": 2.0}  # gas: no ring below ~10 GeV for K (threshold 9.6 GeV); aerogel separates from ~2 GeV
for ax, r in zip(axes, ("gas", "aerogel")):
    key = "Gas" if r == "gas" else "Aerogel"
    edges = np.logspace(np.log10(PMIN[r]), np.log10(60), 61); cent = np.sqrt(edges[:-1] * edges[1:])
    pf = np.logspace(np.log10(PMIN[r]), np.log10(60), 400)
    spread = {}
    for s in d:
        x = d[s]; ok = (x["Theta"] > TH_SEP[0]) & (x["Theta"] < TH_SEP[1]) & (x[f"CherenkovRecoAngle{key}"] > 0)
        res = (x[f"CherenkovRecoAngle{key}"] - x[f"CherenkovAngle{key}"])[ok]; pp = x["P"][ok]
        idx = np.digitize(pp, edges) - 1
        sd = np.full(len(cent), np.nan); n = np.zeros(len(cent))
        for i in range(len(cent)):
            v = res[idx == i]; n[i] = len(v)
            if len(v) >= 30:
                sd[i] = v.std()
        spread[s] = (sd, n)
    for a, b in (("pi", "K"), ("K", "p")):
        dth = np.abs(card_angle(r, cent, MASS[a]) - card_angle(r, cent, MASS[b]))
        sa, na = spread[a]; sb, nb = spread[b]
        ns = dth / np.sqrt(sa ** 2 + sb ** 2)
        # error from the rms uncertainty, sigma_s = s / sqrt(2N)
        ens = ns * np.sqrt((sa ** 4 / (2 * na) + sb ** 4 / (2 * nb))) / (sa ** 2 + sb ** 2)
        ok = np.isfinite(ns) & (ns > 0)
        ax.errorbar(cent[ok], ns[ok], yerr=ens[ok], fmt="o", color=COL[f"{a}/{b}"], ms=4, capsize=0, label=f"{pair_tex(a, b)} Delphes ({TEXV})")
        # card formula, 400 points
        sga, _ = card_sigma(r, pf, MASS[a], math.radians(86.4)); sgb, _ = card_sigma(r, pf, MASS[b], math.radians(86.4))
        nf = np.abs(card_angle(r, pf, MASS[a]) - card_angle(r, pf, MASS[b])) / np.sqrt(sga ** 2 + sgb ** 2)
        ax.plot(pf, nf, "--", color=COL[f"{a}/{b}"], label=f"{pair_tex(a, b)} card formula ($\\theta = 86.4^\\circ$)")
        rd, rf = crossing(cent, ns), crossing(pf, nf)
        result["reach_3sigma_GeV"][f"{r} {a}/{b}"] = dict(delphes=round(rd, 1), card_formula=round(rf, 1))
        for xr, ls in ((rd, "-"), (rf, "--")):
            if np.isfinite(xr):
                ax.plot([xr, xr], [2.2, 4.1], ls, color=COL[f"{a}/{b}"], lw=1.4)
        # full sim T06b, central barrel cell (82.6-90 deg), Xe 2 bar: only meaningful for the Xe card
        if HAVE_FS:
            fs = fsim(r, 82.6)
            fp = sorted({k[1] for k in fs if k[0] == a} & {k[1] for k in fs if k[0] == b})
            fx, fy, fe = [], [], []
            for p0 in fp:
                ma, sa0, na0 = fs[(a, p0)]; mb, sb0, nb0 = fs[(b, p0)]
                v = abs(ma - mb) / math.hypot(sa0, sb0)
                ev = math.hypot(sa0 / math.sqrt(na0), sb0 / math.sqrt(nb0)) / math.hypot(sa0, sb0) + v * math.sqrt(sa0 ** 4 / (2 * na0) + sb0 ** 4 / (2 * nb0)) / (sa0 ** 2 + sb0 ** 2)
                fx.append(p0); fy.append(v); fe.append(ev)
            if fx:
                ax.errorbar(fx, fy, yerr=fe, fmt="s", mfc="white", color=COL[f"{a}/{b}"], ms=6, capsize=2, label=f"{pair_tex(a, b)} full sim ({FS_TAG}, central cell)")
            result["fullsim_T06b_reach_central_cell_GeV"][f"{r} {a}/{b}"] = round(crossing(fx, fy), 1)
    ax.axhline(3, color="grey", ls=":", lw=1.2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(PMIN[r], 60); ax.set_ylim(0.3, 40)
    ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(r"separation $N_\sigma$")
    ticks = [t for t in (2, 3, 5, 10, 15, 20, 30, 40, 60) if t >= PMIN[r]]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    R = RAD[r]
    ax.set_title(f"ARC {r} ({TEXV}), central barrel cell ${TH_SEP[0]:.1f}^\\circ < \\theta < {TH_SEP[1]:.1f}^\\circ$: $n = {R['n']}$, $N_\\mathrm{{sat}} = {R['nsat']:.4g}$", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
save(fig, "separation_vs_p")
edges = np.logspace(np.log10(0.5), np.log10(60), 61); cent = np.sqrt(edges[:-1] * edges[1:]); pf = np.logspace(np.log10(0.5), np.log10(60), 400)

# ---- 2. PID fractions vs p (ARCpid), Clopper-Pearson bands, full acceptance
fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
edges = np.logspace(0, np.log10(60), 61); cent = np.sqrt(edges[:-1] * edges[1:])  # below 1 GeV: no ring for any species
conf = {}
for ax, s in zip(axes, ("pi", "K", "p")):
    x = d[s]; pid = np.abs(x["ARCpid"]).astype(int); idx = np.digitize(x["P"], edges) - 1
    N = np.array([np.sum(idx == i) for i in range(len(cent))])
    for q, lab in (("pi", r"as $\pi$"), ("K", r"as $K$"), ("p", r"as $p$"), (None, "no ring")):
        k = np.array([np.sum((idx == i) & (pid == (PDG[q] if q else 0))) for i in range(len(cent))])
        f = np.where(N > 0, k / np.maximum(N, 1), np.nan); lo, hi = cp_band(k, np.maximum(N, 1))
        c = COL.get(q, "#7f7f7f")
        ax.fill_between(cent, lo, hi, color=c, alpha=0.25, lw=0, step=None)
        ax.plot(cent, f, "o", color=c, ms=3.5, label=lab)
    tot = len(pid)
    conf[s] = {lab: round(float(np.sum(pid == v) / tot), 4) for lab, v in (("pi", 211), ("K", 321), ("p", 2212), ("none", 0))}
    ax.set_xscale("log"); ax.set_xlim(1, 60); ax.set_ylim(-0.02, 1.02); ax.set_xlabel(r"$p$ [GeV]")
    ax.set_title(f"true ${TEXSP[s]}^+$ ({TEXV}): \\texttt{{ARCpid}}, {tot:,} tracks, " + TH_ALL, fontsize=11); ax.grid(alpha=0.3, which="both")
TK = [1, 2, 3, 5, 10, 20, 30, 60]
for ax in axes:
    ax.set_xticks(TK); ax.set_xticklabels([str(t) for t in TK]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axes[0].set_ylabel("fraction"); fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=False, fontsize=10)
save(fig, "pid_fractions_vs_p")
result["confusion_all_p"] = conf
# confusion in p ranges
rng = {}
for lo_, hi_ in ((1, 10), (10, 20), (20, 30), (30, 45), (45, 60)):
    rng[f"{lo_}-{hi_} GeV"] = {s: {lab: round(float(np.mean(np.abs(d[s]["ARCpid"][(d[s]["P"] >= lo_) & (d[s]["P"] < hi_)]) == v)), 4)
                                for lab, v in (("pi", 211), ("K", 321), ("p", 2212), ("none", 0))} for s in d}
result["confusion_by_p"] = rng

# ---- 3. ARCprobK vs p for true K+ and pi+ (board request): mean and fraction with ARCprobK > 0.9
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
edges = np.logspace(np.log10(2), np.log10(60), 61); cent = np.sqrt(edges[:-1] * edges[1:])  # no kaon ring below the aerogel threshold
for s in ("K", "pi"):
    x = d[s]; ok = x["ARCprobK"] >= 0; idx = np.digitize(x["P"][ok], edges) - 1; pk = x["ARCprobK"][ok]
    m = np.array([pk[idx == i].mean() if np.sum(idx == i) > 30 else np.nan for i in range(len(cent))])
    e = np.array([pk[idx == i].std() / math.sqrt(np.sum(idx == i)) if np.sum(idx == i) > 30 else np.nan for i in range(len(cent))])
    axes[0].errorbar(cent, m, yerr=e, fmt="o", color=COL[s], ms=4, label=f"true ${TEXSP[s]}^+$")
    idx2 = np.digitize(x["P"], edges) - 1
    N = np.array([np.sum(idx2 == i) for i in range(len(cent))]); k = np.array([np.sum((idx2 == i) & (x["ARCprobK"] > 0.9)) for i in range(len(cent))])
    lo, hi = cp_band(k, np.maximum(N, 1))
    nz = k > 0  # log axis: bins without a selected track have no point (and no band)
    axes[1].fill_between(cent, np.where(nz, lo, np.nan), np.where(nz, hi, np.nan), color=COL[s], alpha=0.25, lw=0)
    axes[1].plot(cent[nz], (k / np.maximum(N, 1))[nz], "o", color=COL[s], ms=4, label=f"true ${TEXSP[s]}^+$: " + ("efficiency" if s == "K" else "mis-ID"))
axes[0].set_ylabel(r"$\langle P_K \rangle$ (tracks with a ring)"); axes[1].set_ylabel(r"fraction with $P_K > 0.9$ (all tracks)")
axes[1].set_yscale("log"); axes[1].set_ylim(1e-4, 1.5)
for ax in axes:
    ax.set_xscale("log"); ax.set_xlim(2, 60); ax.set_xlabel(r"$p$ [GeV]"); ax.grid(alpha=0.3, which="both"); ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
axes[0].set_title(f"$P_K$ (\\texttt{{ARCprobK}}) vs $p$ ({TEXV}), " + TH_ALL); axes[1].set_title(r"kaon selection $P_K > 0.9$ (Clopper--Pearson 68\% bands)")
TK = [2, 3, 5, 10, 20, 30, 60]
for ax in axes:
    ax.set_xticks(TK); ax.set_xticklabels([str(t) for t in TK]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
save(fig, "probK_vs_p")

# ---- 4. N_pe and per-track resolution vs theta, folded about 90 deg, against the cell-averaged full sim (T06b bins, x error = bin width)
#         N_pe: pi+ 20-60 GeV in 1-deg bins, rings with >= 3 photons (the same selection as the full sim: pi+ >= 20 GeV, found rings;
#         a 45 GeV K+ is only ~95% saturated in the gas, 1 - (9.6/45)^2).
#         sigma_track: pi+ 28-32 GeV in the full-sim bins; full sim: pi+ 30 GeV.
BIN_EDGES = sorted({float(r["theta_lo"]) for r in R6} | {float(r["theta_hi"]) for r in R6})
GAP = (39.3, 44.2)
def fs_bins(rad, key, sp="pi+", pmin=None, p0=None):
    out = []
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        S = [r for r in R6 if r["radiator"] == rad and abs(float(r["theta_lo"]) - lo) < 0.05 and r["species"] == sp
             and ((pmin is not None and float(r["p_GeV"]) >= pmin) or (p0 is not None and float(r["p_GeV"]) == p0)) and int(r["n_tracks"]) >= 30]
        if not S:
            continue
        w = np.array([int(r["n_tracks"]) for r in S], float); v = np.array([float(r[key]) for r in S])
        tm = float(np.mean([float(r["theta_mean"]) for r in S]))
        if key == "sigma_thetac_mrad":
            out.append((tm, lo, hi, float(np.average(v, weights=w)), float(np.average(v, weights=w) / math.sqrt(2 * w.sum()))))
        else:
            out.append((tm, lo, hi, float(np.average(v, weights=w)), float(math.sqrt(np.average(v, weights=w) / w.sum()))))
    return out
x = d["pi"]; sel = (x["P"] > 20)
ee = np.linspace(10, 90, 81); ec = 0.5 * (ee[1:] + ee[:-1]); NB = len(ec)
fig, axes = plt.subplots(2, 2, figsize=(14, 9.5))
npe_th = {}
for j, (r, key) in enumerate((("gas", "Gas"), ("aerogel", "Aerogel"))):
    ax = axes[0, j]
    v = x[f"Npe{key}"][sel]; th = np.minimum(x["Theta"][sel], 180 - x["Theta"][sel]); idx = np.digitize(th, ee) - 1
    ok = v >= MINPE
    m = np.array([v[(idx == i) & ok].mean() if np.sum((idx == i) & ok) > 20 else np.nan for i in range(NB)])
    e = np.array([v[(idx == i) & ok].std() / math.sqrt(max(np.sum((idx == i) & ok), 1)) for i in range(NB)])
    ax.axvspan(*GAP, color="#bbbbbb", alpha=0.45, lw=0, label=r"barrel--endcap gap, $39.3^\circ$--$44.2^\circ$")
    ax.errorbar(ec, m, yerr=e, fmt="o", color="#1f5f8a", ms=3.5, label=f"Delphes, $\\pi^+$ 20--60 GeV ({TEXV})")
    if HAVE_FS:
        fb = fs_bins(r, "mean_npe_found", pmin=20)
        ax.errorbar([q[0] for q in fb], [q[3] for q in fb], xerr=[[q[0] - q[1] for q in fb], [q[2] - q[0] for q in fb]], yerr=[q[4] for q in fb],
                    fmt="s", mfc="white", color="#d62728", ms=6, label=f"full sim {FS_TAG}, $\\pi^+ \\geq 20$ GeV (cell-averaged)")
    result.setdefault("npe_vs_theta_pi20", {})[r] = [(round(float(t), 2), round(float(v), 3)) for t, v in zip(ec, m) if np.isfinite(v)]
    if np.isfinite(m).any():
        npe_th[r] = dict(at_86deg=round(float(np.nanmean(m[-7:])), 2), max=round(float(np.nanmax(m)), 2), theta_of_max=round(float(ec[int(np.nanargmax(m))]), 1))
    ax.set_xlim(10, 92); ax.set_ylim(0, None); ax.set_xlabel(r"$\theta$ [deg] (folded)"); ax.set_ylabel(r"photons per ring $N_\mathrm{pe}$ ($\geq 3$)")
    ax.set_title(f"{r}: photon yield vs $\\theta$"); ax.grid(alpha=0.3); ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    # per-track resolution in the full-sim bins
    ax = axes[1, j]
    xp = d["pi"]; sp = (xp["P"] > 28) & (xp["P"] < 32) & (xp[f"CherenkovRecoAngle{key}"] > 0)
    thp = np.minimum(xp["Theta"][sp], 180 - xp["Theta"][sp]); res = (xp[f"CherenkovRecoAngle{key}"] - xp[f"CherenkovAngle{key}"])[sp] * 1e3
    dl = []
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        q = res[(thp >= lo) & (thp < hi)]
        if len(q) > 50:
            q16, q84 = np.percentile(q, [16, 84]); sd = (q84 - q16) / 2
            dl.append((0.5 * (lo + hi), lo, hi, sd, sd / math.sqrt(2 * len(q))))
    ax.axvspan(*GAP, color="#bbbbbb", alpha=0.45, lw=0)
    ax.errorbar([q[0] for q in dl], [q[3] for q in dl], xerr=[[q[0] - q[1] for q in dl], [q[2] - q[0] for q in dl]], yerr=[q[4] for q in dl],
                fmt="o", color="#1f5f8a", ms=4, label=f"Delphes, $\\pi^+$ 28--32 GeV ({TEXV})")
    if HAVE_FS:
        fb = fs_bins(r, "sigma_thetac_mrad", p0=30.0)
        ax.errorbar([q[0] for q in fb], [q[3] for q in fb], xerr=[[q[0] - q[1] for q in fb], [q[2] - q[0] for q in fb]], yerr=[q[4] for q in fb],
                    fmt="s", mfc="white", color="#d62728", ms=6, label=f"full sim {FS_TAG}, $\\pi^+$ 30 GeV (PDE sampled)")
        result.setdefault("sigma_track_pi30_vs_theta", {})[r] = dict(delphes=[(round(q[0], 1), round(q[3], 3)) for q in dl], fullsim=[(round(q[0], 1), round(q[3], 3)) for q in fb])
    ax.set_xlim(10, 92); ax.set_ylim(0, None); ax.set_xlabel(r"$\theta$ [deg] (folded)"); ax.set_ylabel(r"per-track $\sigma_{\theta_c}$ [mrad] (robust)")
    ax.set_title(f"{r}: per-track resolution vs $\\theta$"); ax.grid(alpha=0.3); ax.legend(fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False)
save(fig, "npe_sigma_vs_theta")
result["npe_vs_theta_K45"] = npe_th

# ---- 5. tracking sanity: tracks per event, sigma(d0), sigma(pT)/pT at 45 GeV, 62.5 < theta < 117.5 deg
san = {}
for s in d:
    x = d[s]; sel = (x["P"] > 40) & (x["P"] < 50) & (x["Theta"] > 62.5) & (x["Theta"] < 117.5)
    san[s] = dict(tracks_per_event=round(float(d[s]["ntrk_per_event"].mean()), 3), n_45GeV=int(sel.sum()),
                  ErrorD0_um=round(float(np.median(x["ErrorD0"][sel]) * 1e3), 2), ErrorPT_over_PT=float("%.3g" % np.median(x["ErrorPT"][sel] / x["PT"][sel])),
                  d0_rms_um=round(float(np.std(x["D0"][sel]) * 1e3), 2))
result["tracking_45GeV_theta62-118"] = san
# ---- 6. Z -> qq full chain
z = f"{W}/runs/{V}_zqq.root"
if os.path.exists(z):
    df = ROOT.RDataFrame("Delphes", z)
    zz = df.AsNumpy(["Track_size", "Jet_size", "MissingET.MET", "Jet.BTag", "Track.ARCpid", "Track.PID", "Track.P"])
    ntr = np.asarray(zz["Track_size"], float); nj = np.asarray(zz["Jet_size"], float)
    met = np.concatenate([np.asarray(v) for v in zz["MissingET.MET"]])
    bt = np.concatenate([np.asarray(v) for v in zz["Jet.BTag"]]) if nj.sum() else np.array([])
    pid = np.concatenate([np.asarray(v) for v in zz["Track.ARCpid"]]); tpid = np.concatenate([np.asarray(v) for v in zz["Track.PID"]])
    kk = np.abs(tpid) == 321
    result["zqq"] = dict(events=len(ntr), tracks_per_event=round(float(ntr.mean()), 2), jets_per_event=round(float(nj.mean()), 2),
                         met_mean_GeV=round(float(met.mean()), 2), btag_fraction=round(float((bt > 0).mean()), 3) if bt.size else None,
                         frac_tracks_with_ring=round(float(np.mean(pid != 0)), 3), true_K_tagged_K=round(float(np.mean(np.abs(pid[kk]) == 321)), 3),
                         true_pi_tagged_K=round(float(np.mean(np.abs(pid[np.abs(tpid) == 211]) == 321)), 4))

json.dump(result, open(f"{W}/validation_{V}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in result.items() if k != "card"}, indent=1))
