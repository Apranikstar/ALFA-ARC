# Results (final versions: reconstruction v2, cards v4)

Full ALFA_o1_v01 simulation, ARC with xenon at 2 bar (main) or C4F10, 1 cm aerogel, SiPM PDE sampled per photon,
cell-averaged (theta uniform, whole-cell bins).

## Per-photon resolution of the gas ring (Xe 2 bar, central barrel cell, pi+ 30-40 GeV, photon truth)

| step | s1 [mrad] |
|---|---|
| photon truth (chromatic dispersion + scattering) | 0.86 |
| + mirror reflection model | 0.91 |
| + 0.8 mm pixels | 1.57 |
| + emission point unknown (mid-depth; tracks bend in 3 T and cross tilted cells) | **2.11** (final reconstruction) |
| pixel centre 1 mm behind the sensor face (earlier reconstruction) | 3.16 |

Per-track resolution at beta = 1: gas 0.70 mrad, aerogel 1.46 mrad (central barrel cell).

## Delphes cards v4 vs the full simulation (3-sigma reach, central barrel cell, GeV)

| | gas pi/K | gas K/p | aerogel pi/K | aerogel K/p |
|---|---|---|---|---|
| Xe 2 bar, Delphes | **24.7** | 42.3 | 9.2 | 15.9 |
| Xe 2 bar, card formula | 25.2 | 43.5 | 9.5 | 16.0 |
| Xe 2 bar, full simulation | 24.3 | 44.4 | 9.2 | 16.3 |
| C4F10, Delphes | **25.0** | 44.0 | 9.6 | 16.0 |
| C4F10, full simulation | 25.2 | 45.5 | 9.1 | 15.6 |

Separation = |theta_a - theta_b| / sqrt(s_a^2 + s_b^2) with s the rms of (reco - true) per species. Delphes reproduces the full
simulation within 2-5%. Z -> qq: true pions tagged as kaons 0.49% (Xe), 0.56% (C4F10). Xenon at 2 bar and C4F10 give the same
reach (same refractive index; the chromatic difference is small against the other terms); xenon costs 1.6% X0 more material.

Card parameters: GasSigmaSingle 1.837 mrad (Xe 2 bar), 1.892 mrad (C4F10), smooth theta dependence (even polynomials in
cot(theta) / tan(theta)), per-track correlated term, yield and ring-finding efficiency vs theta; BendingCoefficient 0.155 / 0.147.

## Material in front of the ECAL (x/X0 at theta = 90 deg)

| | Delphes card | full simulation (k4geo scan) |
|---|---|---|
| beam pipe | 0.68% | 0.68% |
| vertex + outer tracker | 2.04% | 1.78% (no vertex detector in the geometry yet; outer tracker 1.77 vs 1.78) |
| ARC, Xe 2 bar | 5.80% | 5.84% |
| ARC, C4F10 | 4.18% | 4.19% |

## Tracking (Delphes ALFA card, 3 T, all-MAPS, theta = 90 deg)

| p [GeV] | sigma(d0) [um] | sigma(z0) [um] | sigma(pT)/pT |
|---|---|---|---|
| 1 | 11.2 | 11.1 | 0.061% |
| 5 | 3.19 | 2.94 | 0.087% |
| 10 | 2.06 | 1.91 | 0.096% |
| 45 | 1.44 | 1.41 | 0.113% |

Agrees with the ALFA DELPHES curves of the Reco workshop (26 Oct 2026) within 3% at 1 and 182 GeV. IDEA card, same method:
sigma(d0) 14.2 / 4.45 / 2.95 / 1.95 um (ALFA better), sigma(pT)/pT better than ALFA below ~30 GeV (gas tracker).

## Higgs jet-flavour tagging (ParT on ZH, H -> bb/cc/ss/uu/dd/gg)

With the v4 cards: production, ntuples and flat files done (3M jets per card); the three trainings (no ARC PID, ARC Xe 2 bar,
ARC C4F10) are running. ROC curves and mistag rates will be added here.

## What is still idealised

Mirror reflectivity 100% (90% costs ~10% of the photons), track direction from the truth helix (SigmaTrack = 0), no SiPM dark
counts, aerogel n = 1.019 (k4geo) vs 1.03 in some designs, no vertex detector or tracker endcaps in the full-sim geometry yet.
