# 04_delphes: fast simulation of the ALFA ARC

**What it does.** Delphes gets an `ARCPID` module that, for every charged track, predicts what the ARC would measure
(photons per ring, Cherenkov angle with its resolution in each radiator, ring found or not) from formulas fitted to the full
simulation, and computes the pi / K / p probabilities from a chi2 over both radiators. The ALFA cards hold those formulas.

| file / folder | role |
|---|---|
| `delphes-patches/` | Delphes on `delphes/delphes` master `b0465ee`: the ARCPID module, the ARC members of the Track class, the ALFA cards and the card generator, and a fix of the track covariance (`SolTrack::KalmanCovT`: consistent Kalman gain, Joseph-form update). Same commits as github.com/apranikstar/delphes, branch `alfa-arc` |
| `k4simdelphes-patches/` | k4SimDelphes (EDM4hep output) on v00-09: the ARC result of every track written to `EFlowTrack_PID` (type 1023; parameters P(pi), P(K), P(p), angles, resolutions, photons), CMake fixes, `examples/ALFA`. Same as github.com/apranikstar/k4SimDelphes, branch `alfa` |
| `delphes_card_ALFA_Xe2bar.tcl`, `delphes_card_ALFA_C4F10.tcl` | the cards (main: Xe 2 bar) with the ray-tracing ARCFastSim module (patches 0013-0014, `arc/`: cell optics tables, optical tables, dump and card generator); the previous calibrated-formula cards are `*_ARCPIDformulas.tcl`. Tracker: 5-layer MAPS vertex detector, 4 double-layer outer-tracker stations, 3 T; ARC: barrel r = 1.9-2.1 m, endcaps |
| `make_cards.py` | writes both cards from the two calibrations (`03_reconstruction_calibration`): the ARCPID block and the per-gas ARC material |
| `validation/` | 1M pi / K / p per card + Z -> qq through Delphes; `validate.py` compares the separation power, PID fractions, photons and resolution vs theta with the card's own formula and the full simulation |
| `tracking_check/` | the track resolution of the card computed with Delphes' own covariance code on a (p, theta) grid, compared with the ALFA Reco-workshop curves and with the IDEA card |
| `material_delphes_vs_fullsim.py` | x/X0 in the card vs the k4geo material scan (pipe, tracker, ARC; both gases) |

```bash
git clone https://github.com/delphes/delphes && cd delphes && git checkout -b alfa-arc b0465ee && git am ../04_delphes/delphes-patches/*.patch
../env/nightly.sh "./configure && make -j16 HAS_PYTHIA8=true"     # DelphesPythia8 links Delphes statically: rebuild after any change
../env/nightly.sh "python3 examples/ALFA/make_cards.py"           # the cards, from the calibration json
bash ../04_delphes/validation/run_validation.sh Xe2bar C4F10      # then: python3 validation/validate.py Xe2bar
```

**Validation criterion.** The 3-sigma pi/K and K/p reach of Delphes must reproduce the full simulation (same definition:
|theta_a - theta_b| / sqrt(s_a^2 + s_b^2), s = rms of reco - true per species); the current cards agree within 2-5%.
