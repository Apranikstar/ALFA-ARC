"""T11: write cards/delphes_card_ALFA_{Xe2bar,C4F10}.tcl from the bundle card, replacing the ARCPID radiator block with the
k4geo full-simulation calibration v3 (/eos/user/h/hfatehi/alfa/work/T06b/calib_t06b.json).

Calibration v3 is physical:
  - full sim with theta uniform and uniform phi (T06b), analysed in theta bins of whole cells (barrel z-cells of 25.7 cm,
    endcap hexagons of ~0.25 m), so every number is a cell average;
  - each photon kept with probability PDE (a SiPM detects or not), never PDE-weighted;
  - per-photon resolution = pixel (+) chromatic (+) emission; the chromatic term is computed by ARCPID from dn/dE and sigma_E,
    the formula below is the rest (pixel (+) emission); a per-track correlated term (does not average down with N_pe);
  - smooth functions of the cell off-axis variable x = cot(theta) = sinh|eta| (barrel) or tan(theta) = 1/sinh|eta| (endcap),
    and measured constants in the barrel-endcap gap (theta 39.3-44.2 deg: no endcap cells beyond r = 1.7 m).
Usage: python3 examples/ALFA/make_cards.py [Xe2bar] [C4F10]   (default both; only the named cards are written)
Calibration per gas: Xe2bar from work/T06b, C4F10 from work/T06c (the same cell-averaged full sim with C4F10 in the ALFA ARC);
without the T06c calibration the C4F10 card falls back to the Xe calibration x 1.11 on the non-chromatic part."""
import os
import sys
import json
import math

# calibration v4 (T17): reconstruction v2 (pixel centre on the sensor face, T15); v3.1 used work/T06b, work/T06c
CALIB = {"Xe2bar": "/eos/user/h/hfatehi/alfa/campaign_v2/calib/Xe2bar/calib_t06b.json", "C4F10": "/eos/user/h/hfatehi/alfa/campaign_v2/calib/C4F10/calib_t06b.json"}
CAL = json.load(open(CALIB["Xe2bar"]))
ETA_B = -math.log(math.tan(math.radians(44.2) / 2))   # last barrel cell edge
ETA_E = -math.log(math.tan(math.radians(39.3) / 2))   # outer edge of the endcap cells
CHROM0 = {"gas": {"Xe2bar": 0.81, "C4F10": 0.29}, "aerogel": 2.68}  # chromatic term at beta = 1 [mrad], computed by ARCPID itself
C4F10_NONCHROM = 1.11  # standalone cell (no field), PDE sampled: non-chromatic s1 C4F10 1.94 vs Xe 1.75 mrad


def region(rad, key, expr, gap):
    """TFormula: barrel expression in X = sinh|eta| (cot theta), endcap expression in X = 1/sinh|eta| (tan theta), each with X held inside
    the fitted range of that region (no extrapolation), and the measured gap constant in between."""
    xs = {"barrel": "sinh(abs(eta))", "endcap": "(1.0/sinh(abs(eta) + 1e-9))"}
    parts = {}
    for reg in ("barrel", "endcap"):
        f = CAL[rad]["fits"][f"{key}_{reg}"]
        X = f"TMath::Min(TMath::Max({xs[reg]}, {f['xmin']:.4f}), {f['xmax']:.4f})"
        parts[reg] = expr(f["p"]).replace("X", X)
    return ("{\n"
            f"    (abs(eta) <= {ETA_B:.4f}) * ({parts['barrel']}) +\n"
            f"    (abs(eta) > {ETA_B:.4f} && abs(eta) <= {ETA_E:.4f}) * ({gap:.4g}) +\n"
            f"    (abs(eta) > {ETA_E:.4f}) * ({parts['endcap']})\n  }}")


def poly(c):
    """even polynomial c0 + c1 X^2 + c2 X^4 as a TFormula string (X substituted later)"""
    return " + ".join(f"({ci:.6g})" + (f" * pow(X, {2 * i})" if i else "") for i, ci in enumerate(c))


def gap_point(rad, key):
    """measured value in the gap bin; if the gap has too few rings for a fit (gas resolution), the mean of the two neighbouring bins"""
    pts = sorted(CAL[rad]["points"], key=lambda q: q["theta"])
    i = next(k for k, q in enumerate(pts) if q["region"] == "gap")
    v = pts[i].get(key, float("nan"))
    if not math.isfinite(v):
        v = 0.5 * (pts[i - 1][key] + pts[i + 1][key])
    return v


def quad(rad, key, sub2=0.0, scale=1.0):
    """1e-3 * scale * sqrt(poly(X) - sub2) [rad]: the fits are of q^2 in mrad^2."""
    g = gap_point(rad, key)
    gap = 1e-3 * scale * math.sqrt(max(g * g - sub2, 0.0)) if math.isfinite(g) else 0.0
    return region(rad, key, lambda c: f"{1e-3 * scale:.4g} * sqrt(TMath::Max(0.0, {poly(c)} - ({sub2:.5g})))", gap)


def yield_scale(rad):
    return region(rad, "Y", lambda c: f"TMath::Max(0.0, {poly(c)})", gap_point(rad, "Y"))


def efficiency(rad):
    return region(rad, "E", lambda c: f"TMath::Max(0.0, TMath::Min(1.0, {poly(c)}))", gap_point(rad, "E"))


def fits_comment(rad):
    lines = []
    for k, f in CAL[rad]["fits"].items():
        lines.append(f"  ##   {k}: even polynomial in x, c = {', '.join(f'{v:.4g}' for v in f['p'])} (chi2/ndf {f['chi2ndf']:.2f}, x in [{f['xmin']:.3f}, {f['xmax']:.3f}])")
    return "\n".join(lines)


def bending():
    ks = [q["k"] for q in CAL["gas"]["points"] if q.get("k") and q.get("ke") and q["region"] != "gap"]
    kw = [1 / q["ke"] ** 2 for q in CAL["gas"]["points"] if q.get("k") and q.get("ke") and q["region"] != "gap"]
    return sum(k * w for k, w in zip(ks, kw)) / sum(kw)


src = open("cards/delphes_card_ALFA_template.tcl").read()
src = src[src.index("#####################################################################\n", 1) + 70:] if src.startswith("#####################################################################\n# TEMPLATE") else src  # drop the template warning
i = src.index("  ## gas radiator:"); j = src.index("\n", src.index("set AerogelEfficiencyFormula")) + 1
GAS_HDR = {
 "Xe2bar": """  ## gas radiator: xenon at 2 bar, 293.15 K (k4geo Xenon_ARC, ALFA_o1_v01). n and dn/dE at the PDE-weighted mean detected photon
  ## energy <E> = 3.02 eV (T11prep), sigma_E = rms of the detected gas spectrum (T12). Thresholds pi/K/p = 2.71 / 9.59 / 18.2 GeV.
  set GasRefractiveIndex 1.00132
  set GasDnDe 6.17E-5
  set GasSigmaE 0.68
""",
 "C4F10": f"""  ## gas radiator: C4F10 at 1 bar (k4geo C4F10_PFRICH table), for comparison with the July 2026 study. n and dn/dE at the
  ## PDE-weighted mean detected photon energy <E> = 3.05 eV, sigma_E = rms of the detected spectrum (T11prep, standalone ARC cell).
  ## @C4F10_CALIB@
  set GasRefractiveIndex 1.001333
  set GasDnDe 2.25E-5
  set GasSigmaE 0.667
"""}
# ARC passive material in the TrackCovariance table (outside the last measurement: no effect on the track fit, kept right
# for bookkeeping). X0 of the 20 cm ARC layer from the k4geo material scan at 90 deg: Xe 2 bar 5.84 % -> 3.45 m (template),
# C4F10 4.19 % -> 4.78 m (T04 scan, gas layer scaled by X0_Xe / X0_C4F10; T17 material comparison).
ARC_X0 = {"C4F10": 4.78}


def arc_material(card, x0):
    lines = []
    for l in card.split("\n"):
        x = l.split()
        if len(x) >= 13 and x[1] in ("ARCB", "ARCE") and not l.lstrip().startswith("#"):
            x[6] = f"{x0:g}"
            l = "      " + " ".join(x)
        elif "# ARC barrel as passive material" in l:
            l = l.replace("5.8% X0 at 90 deg", "4.2% X0 at 90 deg with C4F10")
        lines.append(l)
    return "\n".join(lines)


for v in (sys.argv[1:] or ["Xe2bar", "C4F10"]):
    own = os.path.exists(CALIB[v])
    CAL = json.load(open(CALIB[v] if own else CALIB["Xe2bar"]))
    K_BEND = bending()
    cg = CHROM0["gas"][v if own else "Xe2bar"]  # remove the chromatic part of the calibrated gas, ARCPID adds the card gas's own
    scale_nc = 1.0 if (own or v == "Xe2bar") else C4F10_NONCHROM
    block = GAS_HDR[v] + f"""  ## per-photon resolution without the chromatic part [rad] (pixel 1.25 (+) emission), cell-averaged full sim (T06b sample, reconstruction v2), smooth in the
  ## cell off-axis variable; even-polynomial fits (s1^2 and s_c^2 in mrad^2, chromatic included in s1 and removed here; Y; E).
  ## Residuals of the smooth forms around the cell-averaged points: s1 <= 2%, s_c 3-9%, Y and E <= 3.6% rms (systematic):
{fits_comment("gas")}
  set GasSigmaSingle {1e-3 * math.sqrt(max(CAL['gas']['fits']['s1_barrel']['p'][0] - cg ** 2, 0)) * scale_nc:.4g}
  set GasSigmaSingleFormula {quad("gas", "s1", sub2=cg * cg, scale=scale_nc)}
  ## per-track correlated term [rad] (does not average down with N_pe): cell-position-dependent bias of the optics
  set GasSigmaCorrFormula {quad("gas", "sc")}
  ## photons at beta = 1 for L = L_ref at normal incidence (cell-averaged); ideal mirror, x 0.90 for a 90% mirror (T07)
  set GasNsat {CAL['gas']['Nsat']:.4g}
  set GasLref 0.15
  ## gas crossed at normal incidence through the k4geo cells (T12 part A): 16.9 cm barrel, 16.7 cm endcap
  set GasThickness 0.169
  ## light collection relative to the path scaling N_sat L / L_ref (photons leave the cell at oblique incidence)
  set GasYieldScaleFormula {yield_scale("gas")}

  ## aerogel radiator: k4geo Aerogel_PFRICH (design value 1.03 not used by the full sim). n and dn/dE at the detected <E> = 3.0 eV,
  ## sigma_E = rms of the detected aerogel spectrum (T12, 90 deg). The chromatic term (2.7 mrad at beta = 1) dominates.
  set AerogelRefractiveIndex 1.0193
  # set AerogelRefractiveIndex 1.03   ;# design value
  set AerogelDnDe 7.9E-4
  set AerogelSigmaE 0.655
  ## per-photon resolution without the chromatic part [rad] (pixel 1.23 (+) emission); fits (mrad, chromatic included):
{fits_comment("aerogel")}
  set AerogelSigmaSingle {1e-3 * math.sqrt(max(CAL['aerogel']['fits']['s1_barrel']['p'][0] - CHROM0['aerogel'] ** 2, 0)):.4g}
  set AerogelSigmaSingleFormula {quad("aerogel", "s1", sub2=CHROM0["aerogel"] ** 2)}
  set AerogelSigmaCorrFormula {quad("aerogel", "sc")}
  set AerogelNsat {CAL['aerogel']['Nsat']:.4g}
  set AerogelLref 0.01
  set AerogelThickness 0.0098
  set AerogelYieldScaleFormula {yield_scale("aerogel")}

  ## per-photon spread from the track curving across the radiator: sigma_bend = BendingCoefficient * 0.3 * B * L * sin(theta) / p,
  ## weighted mean of the per-bin gas fits (T06b sample, reconstruction v2); the reconstruction already corrects for the helix, so k < 1/sqrt(12)
  set BendingCoefficient {K_BEND:.3f}

  ## a ring needs at least MinPe detected photons (as in the T06 reconstruction); the yield fluctuates per track (Poisson)
  set MinPe 3.0
  set SigmaTrack 0.0
  set PoissonFluctuation true
  set UseTrueMomentumForPID false

  ## ring-finding efficiency (ring lost: cell edge, mirror gap; < 3 photons is handled by MinPe), cell-averaged
  set GasEfficiencyFormula {efficiency("gas")}
  set AerogelEfficiencyFormula {efficiency("aerogel")}
"""
    note = ("calibrated on its own cell-averaged full sim (T06c sample: C4F10 in the ALFA ARC, PDE sampled, reconstruction v2)." if own else
            f"no C4F10 full sim yet: Xe calibration, non-chromatic part x {C4F10_NONCHROM} (two standalone cells).")
    block = block.replace("@C4F10_CALIB@", note)
    out = src[:i] + block + src[j:]
    if v != "Xe2bar":  # ARC passive layers in the tracker table: 20 cm with the X0 of the card's gas (template: xenon 2 bar)
        out = arc_material(out, ARC_X0[v])
    out = f"# ALFA Delphes card, ARC gas variant: {v} (generated by examples/ALFA/make_cards.py from delphes_card_ALFA_template.tcl, calibration v4: reconstruction v2, T17)\n" + out
    open(f"cards/delphes_card_ALFA_{v}.tcl", "w").write(out)
    print(v, "written; BendingCoefficient", round(K_BEND, 3))
