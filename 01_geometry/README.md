# 01_geometry: ALFA_o1_v01 with the ARC

| folder / file | what it is |
|---|---|
| `k4geo-patches/` | the 5 commits that add ALFA_o1_v01 to k4geo (`git am` on key4hep/k4geo `ef839651`): master file with the ARC barrel and endcaps, ARC compact with the pressurised-xenon radiator and its material table, the SiPM PDE filter plugin, the ddsim steering file, the test. **This is what you build.** |
| `ALFA_o1_v01_Xe2bar/` | the resulting compact files as plain XML, for reading: the **main design**, xenon at 2 bar, 293 K (`Xenon_ARC` in `materials_arc_o1_v01.xml`, Bideau-Mehu dispersion) + 1 cm aerogel. Identical to `FCCee/ALFA/compact/ALFA_o1_v01/` after the patches. |
| `ALFA_o1_v01_C4F10/` | the C4F10 variant: two files on top of the xenon geometry. `ARC_o1_v01_C4F10.xml` = `ARC_o1_v01.xml` with `material="C4F10_PFRICH"` (1 bar) instead of `Xenon_ARC` for the barrel and endcap radiators; `ALFA_o1_v01_C4F10.xml` = the master file including it. Nothing else differs. |
| `SteeringFile_ALFA_o1_v01_ARC.py` | ddsim steering: optical photons in the ARC (DD4hep optical tracker action), tracker and GRAiNITA ECAL as in ALFA_o1_v00, EDM4hep output; `ARC_PDE_IN_SIM=1` applies the SiPM PDE in the simulation (default: in the analysis) |

## Build and view

```bash
git clone https://github.com/key4hep/k4geo && cd k4geo && git checkout -b alfa-arc ef839651 && git am ../01_geometry/k4geo-patches/*.patch
# C4F10 variant: copy the two files next to the xenon ones
cp ../01_geometry/ALFA_o1_v01_C4F10/*.xml FCCee/ALFA/compact/ALFA_o1_v01/
# build and test in the key4hep nightly (env/nightly.sh), then e.g.
ddsim --steeringFile example/SteeringFile_ALFA_o1_v01_ARC.py --compactFile FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml \
      --gun.particle pi+ --gun.energy "40*GeV" --gun.thetaMin "89*deg" --gun.thetaMax "91*deg" -N 10 --outputFile sim.edm4hep.root
python3 utils/dd4hep2root.py -c FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml -o ALFA_o1_v01.root   # open in a ROOT TBrowser
```

ARC envelope: barrel r = 1.90-2.10 m, |z| < 2.20 m; endcaps |z| = 2.00-2.20 m, r = 0.28-1.90 m. Cells from CLD_o3_v01
(`RadiatorCell_FinalOptimisation_o1_v01.xml`), spherical mirrors (R ~ 0.37 m), 80 x 80 mm SiPM arrays with 0.8 mm pixels
(2 mm thick: photons are detected on the face towards the mirror, see 03_reconstruction_calibration). Material at 90 deg:
ARC 5.8% X0 with xenon, 4.2% with C4F10 (10% budget).
