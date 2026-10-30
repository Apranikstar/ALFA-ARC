"""Material in front of the ECAL: Delphes ALFA cards v4 (TrackCovariance DetectorGeometry) vs the k4geo full-simulation
material scan (T04: TGeo rays to the ARC outer surface, averaged over 60 phi; ALFA_o1_v01 has no vertex detector and no
tracker endcaps yet). Delphes: straight lines from the origin; a barrel layer at r contributes w / (X0 sin theta) if
|r cot theta| lies inside its z range, a disk at z contributes w / (X0 |cos theta|) if |z tan theta| lies inside its r range
(evaluated on 400 points). C4F10: full sim with the gas layer scaled by X0(Xe 2 bar) / X0(C4F10) (as the XeVsC4F10 deck).
Output: campaign_v2/figures/T17_material_delphes_vs_fullsim_v2.{pdf,png} and the numbers at 90 deg."""
import ast, json, math, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.environ["PATH"] = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:" + os.environ["PATH"]
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "font.size": 12,
                     "legend.fontsize": 9, "lines.linewidth": 1.8, "axes.grid": True, "grid.alpha": 0.3})
W = "/eos/user/h/hfatehi/alfa/work"; C2 = "/eos/user/h/hfatehi/alfa/campaign_v2"
X0_XE, X0_C4 = 782.4, 3475.0  # cm (k4geo materials as built by DD4hep)
GROUPS = {"pipe": (r"beam pipe", "#0072B2"), "tracker": (r"vertex + outer tracker", "#E69F00"), "arc": (r"ARC", "#009E73")}


def delphes(card, th):
    s = open(f"{C2}/cards/delphes_card_ALFA_{card}.tcl").read()
    g = re.search(r"set DetectorGeometry \{(.*?)\n\s*\}", s[s.index("module TrackCovariance"):], re.S).group(1)
    out = {k: np.zeros_like(th) for k in GROUPS}
    t = np.radians(th); st, ct = np.sin(t), np.cos(t)
    for line in g.splitlines():
        x = line.split("#")[0].split()
        if len(x) < 13: continue
        typ, lab, a, b, pos, w, x0 = int(x[0]), x[1], float(x[2]), float(x[3]), float(x[4]), float(x[5]), float(x[6])
        grp = "pipe" if lab == "PIPE" else "arc" if lab.startswith("ARC") else "tracker"
        if typ == 1:
            z = pos * ct / st; hit = (z > a) & (z < b); out[grp] += np.where(hit, w / x0 / st * 100, 0)
        else:
            with np.errstate(divide="ignore", invalid="ignore"):
                r = pos * st / ct
            hit = (r > a) & (r < b); out[grp] += np.where(hit, w / x0 / np.abs(ct) * 100, 0)
    return out


src = open(f"{W}/T04/arc_material_scan.py").read()
LAYERS = next(ast.literal_eval(n.value) for n in ast.parse(src).body if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "LAYERS")
M = json.load(open(f"{W}/T04/arc_material_budget.json"))
thf = np.degrees(np.arccos(np.array(M["cos_theta"]))); o = np.argsort(thf); thf = thf[o]


def fullsim(gas):
    out = {k: np.zeros_like(thf) for k in GROUPS}
    for l, c, g in LAYERS:
        x = np.array(M["layers"][l])[o] * 100
        if g == "arc" and gas == "C4F10" and "xenon" in l.lower(): x = x * X0_XE / X0_C4
        grp = "pipe" if l.lower().startswith("beam pipe") else "arc" if g == "arc" else "tracker" if l.startswith("OT") else None
        if grp: out[grp] += x
    return out


th = np.linspace(8, 172, 400)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), sharey=True)
nums = {}
for ax, gas, tex in ((axs[0], "Xe2bar", r"Xe 2\,bar"), (axs[1], "C4F10", r"$\mathrm{C_4F_{10}}$")):
    d, f = delphes(gas, th), fullsim(gas)
    for k, (lab, col) in GROUPS.items():
        ax.plot(th, d[k], color=col, label=f"{lab}: Delphes card v4")
        ax.plot(thf[::3], f[k][::3], "o", color=col, ms=3.5, label=f"{lab}: full sim (k4geo scan)")
    ax.plot(th, sum(d.values()), color="k", label=r"total: Delphes card v4")
    ax.plot(thf[::3], sum(f.values())[::3], "o", color="k", ms=3.5, label=r"total: full sim")
    ax.set_xlim(8, 172); ax.set_xlabel(r"$\theta$ [deg]"); ax.set_title(f"{tex}: material up to the ARC outer surface", fontsize=12)
    i90 = np.argmin(abs(thf - 90)); j90 = np.argmin(abs(th - 90))
    nums[gas] = {k: (round(float(d[k][j90]), 2), round(float(f[k][i90]), 2)) for k in GROUPS}
    nums[gas]["total"] = (round(float(sum(d.values())[j90]), 2), round(float(sum(f.values())[i90]), 2))
axs[0].set_ylabel(r"$x / X_0$ [\%]"); axs[0].set_ylim(0, 16)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.1), frameon=False)
fig.text(0.5, -0.135, r"Full sim has no vertex detector and no tracker endcaps yet (ALFA\_o1\_v01); the ARC layers of the card lie outside the last "
         r"measurement and do not enter the Delphes track fit.", ha="center", fontsize=9)
fig.tight_layout()
for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
    fig.savefig(f"{C2}/figures/T17_material_delphes_vs_fullsim_v2.{ext}", bbox_inches="tight", **kw)
json.dump(nums, open(f"{C2}/figures/T17_material_delphes_vs_fullsim_v2.json", "w"), indent=1)
print("x/X0 [%] at 90 deg (Delphes card, full sim):", json.dumps(nums))
