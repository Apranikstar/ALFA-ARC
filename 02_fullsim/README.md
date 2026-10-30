# 02_fullsim: Geant4 simulation of the ARC inside ALFA

**What it does.** Shoots single pi+, K+ and protons through the full ALFA_o1_v01 detector (ddsim, Geant4 with optical
photons) and keeps every photon detected by the ARC SiPMs. This is the "truth" the Delphes cards are calibrated on.

**Why "cell-averaged".** The polar angle is uniform in 12-90 deg and phi is uniform, so every track hits a random spot of a
random cell. The analysis then bins in whole cells, so each number is an average over the cell, as a real detector sees it.

| file | role |
|---|---|
| `run_fullsim.sh` | one job: `run_fullsim.sh <Xe2bar or C4F10> <particle> <p_GeV> <n_events> <seed> <chunk>`: ddsim, then `extract_arc_photons.py`, copies the `.npz` to EOS, deletes the large ROOT file |
| `fullsim.sub`, `grid.txt` | HTCondor submission: 2 gases x 3 species x 10 momenta (1-50 GeV); 4 x 1000 events below 20 GeV, 40 x 100 events at 20-50 GeV (444 jobs per gas, ~10 CPU-min each) |
| `extract_arc_photons.py` | EDM4hep -> small `.npz`: the primary track and every ARC photon hit (cell, pixel x/y, time, energy) |
| `photon_truth/` | special sample (pi+ 30/40/50 GeV, central barrel cell) that also keeps each photon's true emission point and direction (`--part.keepAllParticles`, `extract_arc_photons_truth.py`); used to find what limits the resolution |

**Run**: `condor_submit fullsim.sub` (paths inside point to the original lxplus/EOS locations; adapt `OUT` and the k4geo path).

**Output**: `<species>_p<GeV>_c<chunk>.npz` per job, ~190 MB per gas in total -> input of `03_reconstruction_calibration`.
