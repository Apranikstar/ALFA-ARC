"""T16 plots: ALFA card tracking resolution (KalmanCovT covariance, as Delphes computes it per track) vs the Reco-workshop
DELPHES curves and vs the IDEA card, plus the effect of the Kalman fix. LaTeX, PDF + 600 dpi PNG, legends outside the axes.
Inputs: grid_alfa_new.txt (fixed code, new OT material), grid_alfa_v0.txt (original code, old table), grid_idea.txt,
reco_ref.json (digitised slide 5), gun_alfa_tc_orig.txt / gun_alfa_tc_fix.txt (same 10k tracks, original vs fixed filter)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.environ["PATH"] = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:" + os.environ["PATH"]
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 12,
                     "legend.fontsize": 9, "lines.linewidth": 1.8, "axes.grid": True, "grid.alpha": 0.3})
COL = {1: "#0072B2", 5: "#E69F00", 10: "#E69F00", 46: "#009E73", 45: "#009E73", 182: "#CC79A7"}  # validated 4-slot palette
OUT = "plots"


def grid(fn):
    """{(p, theta): row} from the grid part and {theta: array} from the momentum scans ('S' lines)"""
    g, s = {}, {}
    for l in open(fn):
        if l.startswith("#"): continue
        x = l.split()
        if x[0] == "S":
            s.setdefault(float(x[2]), []).append([float(v) for v in x[1:]])
        else:
            v = [float(t) for t in x]; g[(v[0], v[1])] = v
    return g, {k: np.array(v) for k, v in s.items()}


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig); print("wrote", name)


A, As = grid("grid_alfa_new.txt"); I, Is = grid("grid_idea.txt"); A0, _ = grid("grid_alfa_v0.txt")
ref = json.load(open("reco_ref.json"))
th = np.arange(10, 91)

# 1. vs the Reco-workshop curves: sigma(d0), sigma(z0) vs theta at 1, 5, 46, 182 GeV
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, var, col, lab in [(axs[0], "d0", 4, r"$\sigma(d_0)$ [$\mu$m]"), (axs[1], "z0", 5, r"$\sigma(z_0)$ [$\mu$m]")]:
    for p in (1, 5, 46, 182):
        y = [A[(p, t)][col] for t in th]
        ax.plot(th, y, color=COL[p], ls="-")
        r = np.array(ref[var][str(p)])
        ax.plot(r[:, 0], r[:, 1], color=COL[p], ls="--", lw=1.4)
        ax.text(91.5, y[-1] * {46: 1.12, 182: 0.89}.get(p, 1), rf"{p}\,GeV", color="0.2", va="center", fontsize=9)  # 46/182 GeV curves meet at 90 deg
    ax.set_yscale("log"); ax.set_xlim(10, 90); ax.set_xlabel(r"$\theta$ [deg]"); ax.set_ylabel(lab)
axs[0].set_title(r"ALFA card (3\,T) vs Reco workshop, 26 Oct 2026", fontsize=11)
axs[1].set_title(r"$\sigma(z_0)$; reference digitised up to $80^\circ$", fontsize=11)
h = [plt.Line2D([], [], color=COL[p]) for p in (1, 5, 46, 182)] + [plt.Line2D([], [], color="k"), plt.Line2D([], [], color="k", ls="--", lw=1.4)]
fig.legend(h, [rf"$p = {p}$\,GeV" for p in (1, 5, 46, 182)] + [r"this card (T16 table, fixed Kalman)", r"Reco workshop slide 5 (digitised)"],
           loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), frameon=False)
fig.tight_layout(rect=(0, 0.06, 1, 1)); save(fig, "t16_vs_reco")

# 2. vs momentum at 90 and 45 deg, ALFA vs IDEA
fig, axs = plt.subplots(2, 2, figsize=(11, 7.6))
spec = [(axs[0, 0], 4, r"$\sigma(d_0)$ [$\mu$m]", 1), (axs[0, 1], 5, r"$\sigma(z_0)$ [$\mu$m]", 1),
        (axs[1, 0], 6, r"$\sigma(p_T)/p_T$", 1), (axs[1, 1], 7, r"$\sigma(\theta)$ [mrad]", 1)]
TC = {90.0: "#0072B2", 45.0: "#E69F00"}
for ax, col, lab, _ in spec:
    for t, c in TC.items():
        ax.plot(As[t][:, 0], As[t][:, col], color=c, ls="-")
        ax.plot(Is[t][:, 0], Is[t][:, col], color=c, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 100); ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(lab)
h = [plt.Line2D([], [], color=c) for c in TC.values()] + [plt.Line2D([], [], color="k"), plt.Line2D([], [], color="k", ls=":")]
fig.legend(h, [r"$\theta = 90^\circ$", r"$\theta = 45^\circ$", r"ALFA card (3\,T, all-MAPS)", r"IDEA card (2\,T, drift chamber)"],
           loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.03), frameon=False)
fig.tight_layout(rect=(0, 0.03, 1, 1)); save(fig, "t16_vs_p")

# 3. vs theta at 10 and 45 GeV, ALFA vs IDEA
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, col, lab in [(axs[0], 4, r"$\sigma(d_0)$ [$\mu$m]"), (axs[1], 5, r"$\sigma(z_0)$ [$\mu$m]"), (axs[2], 6, r"$\sigma(p_T)/p_T$")]:
    for p in (10, 45):
        ax.plot(th, [A[(p, t)][col] for t in th], color=COL[p], ls="-")
        ax.plot(th, [I[(p, t)][col] for t in th], color=COL[p], ls=":")
    ax.set_yscale("log"); ax.set_xlim(10, 90); ax.set_xlabel(r"$\theta$ [deg]"); ax.set_ylabel(lab)
h = [plt.Line2D([], [], color=COL[p]) for p in (10, 45)] + [plt.Line2D([], [], color="k"), plt.Line2D([], [], color="k", ls=":")]
fig.legend(h, [r"$p = 10$\,GeV", r"$p = 45$\,GeV", r"ALFA card", r"IDEA card"], loc="lower center", ncol=4,
           bbox_to_anchor=(0.5, -0.07), frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1)); save(fig, "t16_vs_theta")

# 4. the Kalman fix: per-track sigma(pT)/pT, original / fixed filter, same 10k tracks
o = np.loadtxt("gun_alfa_tc_orig.txt"); f = np.loadtxt("gun_alfa_tc_fix.txt")
ok = o[:, 3] == 1
fig, ax = plt.subplots(figsize=(7.5, 4.4))
r = o[ok, 7] / f[ok, 7]
ax.scatter(o[ok, 0], r, s=6, color="#0072B2", lw=0, label=rf"KalmanCovT accepted ({ok.sum()} tracks)")
ax.scatter(o[~ok, 0], np.full((~ok).sum(), 0.3), s=14, marker="v", color="#D55E00", lw=0,
           label=rf"KalmanCovT failed: negative or non-positive-definite covariance ({(~ok).sum()} tracks)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.2, 8); ax.set_xlim(0.1, 60)
ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(r"$\sigma(p_T)/p_T$: original / fixed filter")
ax.set_title(r"ALFA card, 10k $\pi^+$, $|\eta| < 2.6$: same tracks, two covariance updates", fontsize=11)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False, ncol=1)
save(fig, "t16_kalman_fix")
