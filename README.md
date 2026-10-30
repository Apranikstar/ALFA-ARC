# ALFA-ARC

The ARC (Array of RICH Cells) particle identification of the ALFA FCC-ee detector concept, end to end:
**geometry -> full simulation -> photon reconstruction and calibration -> Delphes ARCPID cards -> Higgs jet-flavour tagging**.
Main design: pressurised xenon at 2 bar as the gas radiator (+ 1 cm aerogel); C4F10 as the alternative. Every stage here is the
final version (reconstruction v2 with the pixel centre on the sensor face, Delphes cards v4).

Final numbers: [RESULTS.md](RESULTS.md).

## Layout

| folder | stage |
|---|---|
| `env/` | `nightly.sh`: runs a command in a clean shell with the key4hep nightly and the local k4geo build |
| `01_geometry/` | k4geo patches that add `ALFA_o1_v01` with the ARC barrel and endcaps; the resulting compact files for both gases as plain XML (`ALFA_o1_v01_Xe2bar/`, main; `ALFA_o1_v01_C4F10/`, only the radiator material differs); ddsim steering; see its README |
| `02_fullsim/` | HTCondor full simulation: `T06b` (Xe 2 bar) and `T06c` (C4F10), pi/K/p at 10 momenta, theta uniform 12-90 deg (cell-averaged); `T15` photon-truth sample; extractors to `.npz` |
| `03_reconstruction_calibration/` | per-photon Cherenkov-angle reconstruction, per-track table in cell-sized theta bins, smooth calibration fits |
| `04_delphes/` | Delphes patches (ARCPID module, cards, Kalman covariance fix), k4SimDelphes patches (ARC PID -> EDM4hep), cards v4 and their generator, validation, tracking check, material comparison |
| `05_flavour_tagging/` | ZH -> nunu jj production (k4SimDelphes, EDM4hep), FCCAnalyses (exclusive Durham N = 2 + JetFlavourHelper + ARC features), flattening, ParT (weaver-core) trainings on HTCondor GPUs, ROC curves |

## Software

key4hep nightly (LCG devkey-head, Oct 2026; ROOT 6.40, Geant4 11.4, DD4hep), Pythia 8.316, weaver-core 0.5.4 with
torch 2.11 (LCG_110_cuda) for ParT. Bases of the patches:

- k4geo: `key4hep/k4geo` main `ef839651` + `01_geometry/k4geo-patches/*.patch`
- Delphes: `delphes/delphes` master `b0465ee` + `04_delphes/delphes-patches/*.patch` (also: github.com/apranikstar/delphes, branch `alfa-arc`)
- k4SimDelphes: `key4hep/k4SimDelphes` v00-09 `b85e1d3` + `04_delphes/k4simdelphes-patches/*.patch` (also: github.com/apranikstar/k4SimDelphes, branch `alfa`)

The scripts carry the lxplus/EOS paths of the original run as variables near the top (`/afs/cern.ch/work/h/hfatehi/alfa`,
`/eos/user/h/hfatehi/alfa`, `/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa`); change them to rerun elsewhere.

## Run order

```bash
N=env/nightly.sh            # every step runs in the nightly environment

# 1. geometry
git clone https://github.com/key4hep/k4geo && cd k4geo && git checkout -b alfa-arc ef839651 && git am ../01_geometry/k4geo-patches/*.patch
$N "cd k4geo && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=../install && make -j16 install && ctest -R ALFA_o1_v01"
#    C4F10 variant: cp ../01_geometry/ALFA_o1_v01_C4F10/*.xml FCCee/ALFA/compact/ALFA_o1_v01/
#    view it: python3 utils/dd4hep2root.py -c FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml -o ALFA_o1_v01.root   (ROOT TBrowser)

# 2. full simulation (cell-averaged; ~10 CPU-min per job; the ROOT files are dropped after extraction to .npz)
condor_submit 02_fullsim/T06b/t06b.sub          # Xe 2 bar
condor_submit 02_fullsim/T06c/t06c.sub          # C4F10 (geometry 01_geometry/ALFA_o1_v01_C4F10)

# 3. reconstruction table + calibration (independent processes: DD4hep cannot run inside multiprocessing workers)
cd 03_reconstruction_calibration
for k in $(seq 0 31); do $N "export T06B_DIR=<T06b dir> T06B_OUT=<out>; python3 make_table.py worker $k 32 <out>/parts/part_$k.pkl" & done; wait
$N "export T06B_DIR=<T06b dir> T06B_OUT=<out>; python3 make_table.py merge <out>/parts/part_*.pkl"     # -> arc_table_t06b.csv
$N "export T06B_DIR=<out>; python3 calibrate.py"                                                     # -> calib_t06b.json, plots/
#    C4F10: T06B_DIR=<T06c dir>, T06B_NGAS=1.001333 for the table; T06B_GAS_CHROM=0.29 for the calibration

# 4. Delphes with ARCPID and the cards
git clone https://github.com/delphes/delphes && cd delphes && git checkout -b alfa-arc b0465ee && git am ../04_delphes/delphes-patches/*.patch
$N "cd delphes && ./configure && make -j16 HAS_PYTHIA8=true"    # DelphesPythia8 links Delphes statically: rerun after any change
$N "cd delphes && python3 examples/ALFA/make_cards.py"          # cards v4 from the two calibrations (copies in 04_delphes/)
bash 04_delphes/validation/run_validation.sh Xe2bar C4F10       # 1M pi/K/p per card + Z -> qq
$N "python3 04_delphes/validation/validate.py Xe2bar"           # reach, PID fractions, N_pe and sigma vs theta vs the full sim
#    k4SimDelphes (EDM4hep output with the ARC PID in EFlowTrack_PID): see 04_delphes/k4simdelphes-patches and examples/ALFA/README.md there

# 5. flavour tagging (ZH at 240 GeV, Z -> nunu, H -> bb/cc/ss/uu/dd/gg, 250k events each, both cards with the same seeds)
condor_submit 05_flavour_tagging/condor/prod.sub                # 2 cards x 6 decays x 25 chunks x 10k events -> EDM4hep
$N "export ALFA_CARD=Xe2bar; cd 05_flavour_tagging/fccanalysis && fccanalysis run analysis_jets.py"   # and C4F10
$N "cd 05_flavour_tagging/fccanalysis && python3 flatten_jets.py Xe2bar"                               # one entry per jet, labels
python3 05_flavour_tagging/weaver/patch_weaver_fused.py <weaver site-packages>                        # once (see below)
condor_submit 05_flavour_tagging/condor/train.sub               # ParT: no ARC PID, ARC Xe 2 bar, ARC C4F10 (15 epochs, GPU)
python3 05_flavour_tagging/make_roc_plots.py                    # ROC curves, x signal efficiency (linear), y mistag (log)
```

## Things that were needed (and why)

- **Pixel centre on the sensor face** (`arc_photon_reco.py`, `face_normal`): the cellID pixel centre lies on the mid-plane of the 2 mm SiPM, 1 mm behind the face where photons are detected (on the mirror side; the local frame is flipped between barrel and endcaps). Using the mid-plane inflates the gas per-photon resolution from 2.1 to 3.2 mrad.
- **Kalman covariance in Delphes** (`delphes-patches`, SolTrack::KalmanCovT): the gain was built without the layer's multiple-scattering term and the update lost positive-definiteness with an all-silicon tracker (0.5% of tracks failed, ~5% got a noisy covariance). Now a consistent gain and the Joseph form.
- **ParT training**: weaver-core 0.5.4 reads unpacked jagged buffers across files (`patch_weaver_fused.py`), and ~1 constituent per 200k jets has E < |pz| after rounding (rapidity = NaN): the configs use `pfcand_e_safe = max(E, 1.000001 |p|)`. LCG `setup.sh` fails under `set -u`.
