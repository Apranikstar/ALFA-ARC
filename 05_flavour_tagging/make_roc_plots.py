"""[campaign_v2 copy (T17): only the v2 trainings (v4 cards, Kalman-fixed Delphes); figures T17_r01_roc_<X>_v2,
T17_r02_roc_summary_v2 in campaign_v2/figures; numbers in campaign_v2/tagging/roc_numbers_v2.json]
ROC curves of the ParT jet-flavour taggers (T14) for the XeVsC4F10 deck.
For every signal flavour X in (b, c, s, u, d, g): the discriminant is the tagger's probability P(X); the signal efficiency is the
fraction of true-X jets above a cut, the mistag rate the fraction of true-Y jets above the same cut, one curve per background Y.
x = signal efficiency (linear), y = mistag rate (log). Colour = background flavour (fixed), line style = model.
Inputs: <BASE>/training/<model>/pred.root (weaver --predict-output: score_label_<f>, label_<f>).
Outputs: plots/r01_roc_<X>.{pdf,png} (all models), plots/r02_roc_summary.{pdf,png} (2x3 panels, ARC Xe vs C4F10), roc_numbers.json
Usage: python3 make_roc_plots.py [model=label ...]   default: part_baseline_Xe2bar part_arc_Xe2bar part_arc_C4F10"""
import json, os, sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.environ["PATH"] = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:" + os.environ["PATH"]
plt.rcParams.update({"lines.linewidth": 1.8, "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
                     "font.size": 12, "legend.fontsize": 9})
BASE = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2"
OUT = "/eos/user/h/hfatehi/alfa/campaign_v2/figures"
FL = ["b", "c", "s", "u", "d", "g"]
LAB = {f: f"recojet_is{f.upper()}" for f in FL}
TEX = {"b": "b", "c": "c", "s": "s", "u": "u", "d": "d", "g": "g"}
# fixed colour per flavour (colour-blind-safe Okabe-Ito)
COL = {"b": "#D55E00", "c": "#009E73", "s": "#E69F00", "u": "#0072B2", "d": "#56B4E9", "g": "#CC79A7"}
MODELS = [("part_arc_Xe2bar", r"ParT + ARC, Xe 2\,bar", "-"), ("part_arc_C4F10", r"ParT + ARC, $\mathrm{C_4F_{10}}$", "--"),
          ("part_baseline_Xe2bar", r"ParT, no ARC PID", ":")]
if len(sys.argv) > 1:
    MODELS = [(a.split("=")[0], a.split("=")[1] if "=" in a else a, ls) for a, ls in zip(sys.argv[1:], ["-", "--", ":", "-."])]
EFF = np.linspace(0.0, 1.0, 400)  # signal-efficiency grid for the curves and the numbers


def load(model):
    t = uproot.open(f"{BASE}/training/{model}/pred.root")["Events"]  # weaver --predict-output tree
    a = t.arrays([f"score_{LAB[f]}" for f in FL] + [LAB[f] for f in FL], library="np")
    truth = np.argmax(np.vstack([a[LAB[f]] for f in FL]).T, axis=1)
    score = np.vstack([a[f"score_{LAB[f]}"] for f in FL]).T
    return truth, score


def roc(sig, bkg):
    """mistag rate of bkg at each signal efficiency of EFF, cutting on the score"""
    s = np.sort(sig)[::-1]; b = np.sort(bkg)
    out = np.full_like(EFF, np.nan)
    for i, e in enumerate(EFF):
        if e <= 0 or e > 1:
            continue
        cut = s[min(int(np.ceil(e * len(s))) - 1, len(s) - 1)]
        out[i] = 1.0 - np.searchsorted(b, cut, side="left") / len(b)
    out[out < 5.0 / len(b)] = np.nan  # fewer than 5 background jets passing: not shown (statistics)
    return out


def legend(ax, X, models, ncol=4, fs=8.5):
    """colour = background flavour, line style = model (proxy handles)"""
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=COL[Y], lw=2.2, label=f"${X}$ vs ${Y}$") for Y in FL if Y != X]
    h += [Line2D([], [], color="k", ls=ls, lw=1.8, label=lab) for _, lab, ls in models]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=ncol, frameon=False, fontsize=fs)


data = {m: load(m) for m, _, _ in MODELS if os.path.exists(f"{BASE}/training/{m}/pred.root")}
numbers = {}
for X in FL:
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ix = FL.index(X)
    for m, lab, ls in MODELS:
        if m not in data:
            continue
        truth, score = data[m]
        sig = score[truth == ix, ix]
        for Y in FL:
            if Y == X:
                continue
            bkg = score[truth == FL.index(Y), ix]
            mt = roc(sig, bkg)
            ok = np.isfinite(mt) & (mt > 0)
            ax.plot(EFF[ok], mt[ok], ls, color=COL[Y])
            for e in (0.5, 0.7, 0.9):
                numbers.setdefault(m, {}).setdefault(f"{X}_vs_{Y}", {})[f"eff{int(e * 100)}"] = float(np.interp(e, EFF, np.nan_to_num(mt, nan=1.0)))
    ax.set_yscale("log"); ax.set_xlim(0, 1); ax.set_ylim(1e-4, 1.2)
    ax.set_xlabel(f"${X}$-jet efficiency (signal)"); ax.set_ylabel(r"mistag rate (background jets passing)")
    ax.set_title(f"${X}$-tagging with $P({X})$: ZH, $H \\to jj$, 240\\,GeV, ALFA", fontsize=12)
    ax.grid(alpha=0.3, which="both")
    legend(ax, X, [mm for mm in MODELS if mm[0] in data])
    fig.tight_layout()
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
        fig.savefig(f"{OUT}/T17_r01_roc_{X}_v2.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
# summary: all signal flavours, the two ARC models
fig, axes = plt.subplots(2, 3, figsize=(16, 10.5), gridspec_kw=dict(hspace=0.5))
for ax, X in zip(axes.flat, FL):
    ix = FL.index(X)
    for m, lab, ls in MODELS[:2]:
        if m not in data:
            continue
        truth, score = data[m]; sig = score[truth == ix, ix]
        for Y in FL:
            if Y == X:
                continue
            mt = roc(sig, score[truth == FL.index(Y), ix]); ok = np.isfinite(mt) & (mt > 0)
            ax.plot(EFF[ok], mt[ok], ls, color=COL[Y])
    ax.set_yscale("log"); ax.set_xlim(0, 1); ax.set_ylim(1e-4, 1.2); ax.grid(alpha=0.3, which="both")
    ax.set_xlabel(f"${X}$-jet efficiency"); ax.set_ylabel("mistag rate"); ax.set_title(f"${X}$-tagging", fontsize=12)
    legend(ax, X, [mm for mm in MODELS[:2] if mm[0] in data], ncol=4, fs=8)
fig.suptitle(r"ParT + ARC: solid Xe 2\,bar, dashed $\mathrm{C_4F_{10}}$", fontsize=13, y=0.99)
for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
    fig.savefig(f"{OUT}/T17_r02_roc_summary_v2.{ext}", bbox_inches="tight", **kw)
json.dump(numbers, open("/eos/user/h/hfatehi/alfa/campaign_v2/tagging/roc_numbers_v2.json", "w"), indent=1)
print("models", list(data), "-> r01_roc_<X>, r02_roc_summary, roc_numbers.json")
