# 05_flavour_tagging: what the ARC brings to Higgs jet-flavour tagging

**Physics sample.** e+e- -> ZH at 240 GeV, Z -> nu nu, H -> bb, cc, ss, uu, dd, gg; 250k events per decay (25 chunks x 10k),
Pythia 8, the same seeds for every detector, so the detectors see identical generated events.

**Chain.**

1. **Production** (`run_production.sh` locally or `condor/prod.sub`): Pythia 8 + Delphes with the ALFA card -> EDM4hep via
   k4SimDelphes (`delphes_edm4hep.sh`, output config `cfg/edm4hep_output_config_ALFA.tcl`), with the ARC result per track.
2. **Jets** (`fccanalysis/analysis_jets.py`, FCCAnalyses): all reconstructed particles clustered into exactly 2 jets
   (exclusive Durham / ee_kt, `ExclusiveJetClusteringHelper`); per constituent the standard features (`jetFlavourHelper.py`:
   kinematics, impact parameters and their significances, track covariance, particle type) plus the ARC PID
   (P(pi), P(K), P(p), gas / aerogel ring seen; `alfa_jet_utils.h`).
3. **Flatten** (`fccanalysis/flatten_jets.py`): one entry per jet with its label (`recojet_isB` ... `recojet_isG`); chunks 0-19
   for training, 20-24 for testing (3M jets per detector).
4. **Train** (`weaver/`, `condor/train.sub`): Particle Transformer (weaver-core 0.5.4), 6 classes, 15 epochs, on HTCondor GPUs:
   `jets_baseline.yaml` (no ARC PID) and `jets_arc.yaml` (with ARC PID).
5. **ROC curves** (`make_roc_plots.py`): for each flavour, the mistag rate of every other flavour vs the signal efficiency
   (x linear, y log), cutting on the tagger output P(X).

**Fixes needed for finite ParT losses.** `weaver/patch_weaver_fused.py` (weaver-core 0.5.4 mixes unpacked jagged buffers
across files) and `pfcand_e_safe = max(E, 1.000001 |p|)` in the configs (about 1 constituent in 200k jets along the beam
axis has E < |pz| after rounding, which makes ParT's rapidity NaN). `weaver/check_inputs.py` and `nan_steps.py` are the probes.
