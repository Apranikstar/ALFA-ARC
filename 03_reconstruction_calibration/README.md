# 03_reconstruction_calibration: from photon hits to the numbers in the Delphes card

**Step 1, per photon** (`arc_photon_reco.py`). For every detected photon: take the pixel centre on the sensor face, the
track helix in the 3 T field, assume emission at mid-depth of the radiator, solve the single reflection on the cell's
spherical mirror (refraction at the aerogel surface for aerogel photons) and get the Cherenkov angle theta_c.
Per track and radiator: the mean of the photons near the expected angle, with two 3-sigma rejection passes.

**Step 2, per theta bin** (`make_table.py`). Runs step 1 on all full-sim events (32 independent processes; DD4hep cannot
run inside a multiprocessing pool), samples each photon with the SiPM detection efficiency, and writes one row per
(species, momentum, cell-sized theta bin, radiator): photons per ring, ring-finding efficiency, per-photon resolution s1,
per-track resolution, bias. -> `arc_table_t06b.csv`

**Step 3, smooth calibration** (`calibrate.py`). Fits smooth, physical functions of the cell off-axis variable
(cot theta in the barrel, tan theta in the endcap) to the table: per-photon resolution, per-track correlated term, photon
yield and ring-finding efficiency, separately for barrel and endcap, with the barrel-endcap gap kept at its measured value.
-> `calib_t06b.json` + plots. These fits become the formulas of the Delphes card (`04_delphes/make_cards.py`).

`arc_ring_analysis.py` holds the shared constants: refractive indices of the gas and aerogel averaged over the detected
photon spectrum, particle masses, helpers.

```bash
for k in $(seq 0 31); do ../env/nightly.sh "export T06B_DIR=<npz dir> T06B_OUT=<out>; python3 make_table.py worker $k 32 <out>/parts/part_$k.pkl" & done; wait
../env/nightly.sh "export T06B_DIR=<npz dir> T06B_OUT=<out>; python3 make_table.py merge <out>/parts/part_*.pkl"
../env/nightly.sh "export T06B_DIR=<out>; python3 calibrate.py"
# C4F10: T06B_NGAS=1.001333 for make_table, T06B_GAS_CHROM=0.29 for calibrate
```

**Two details that matter.** (1) The pixel centre from the cellID sits on the mid-plane of the 2 mm SiPM, 1 mm behind the
face where photons are detected; the reconstruction moves it to the face towards the mirror (`face_normal`). Without this the
gas resolution is 3.2 instead of 2.1 mrad. (2) The SiPM efficiency is applied by sampling (a photon is detected or not),
never as a weight: weighting underestimates the per-track spread.
