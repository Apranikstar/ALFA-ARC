"""ALFA Delphes cards with the ray-tracing ARCFastSim module in place of the calibrated ARCPID block.

The module instance keeps the name ARCPID so the execution path and TrackMerger are unchanged. The base card of each
configuration supplies the tracker and the ARC material; only the ARC response block is replaced.
Usage: python3 make_fast_cards.py <A|B|C> <Xe2bar|C4F10> <out.tcl> [--production]
   validation cards (default): TrueHypothesisAngle true (reco angle in the window around the true angle, as the
   full-simulation analysis); --production: false (the best hypothesis), UseTrueTrack false (smeared track in the reco)."""
import re
import sys

X = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v3"
HERE = "/afs/cern.ch/work/h/hfatehi/alfa/campaign_v3_code/arcfast"
BASE = {"A": f"{X}/configA/cards_v42", "B": f"{X}/cards_v52", "C": f"{X}/configC/cards_C2"}
CELLS = {"A": "cells_ALFA_o1_v01.txt", "B": "cells_ALFA_o1_v01.txt", "C": "cells_ALFA_o1_v01_bigSiPM.txt"}
AERO_SCALE = {"A": 1.0 / 1.55285, "B": 1.0, "C": 1.0}  # n - 1 relative to the k4geo table (n = 1.030 at 400 nm)
GAS = {"Xe2bar": ("RINDEX__Xenon_ARC", ""), "C4F10": ("RINDEX__C4F10_PFRICH", "")}
NOTE = {"A": "aerogel n = 1.019 (table scaled), 80 mm SiPM arrays", "B": "aerogel n = 1.030, 80 mm SiPM arrays",
        "C": "aerogel n = 1.030, SiPM arrays 89.6 mm (barrel) / 108.8 mm (endcap)"}

cfg, gas, out = sys.argv[1:4]
production = "--production" in sys.argv
card = open(f"{BASE[cfg]}/delphes_card_ALFA_{gas}.tcl").read()
m = re.search(r"^module ARCPID ARCPID \{.*?^\}\n", card, re.S | re.M)
assert m, "ARCPID block not found"
block = f"""module ARCFastSim ARCPID {{
  ## ARC fast simulation by ray tracing (ARCFastSim): photons generated along the helix in the aerogel and the gas from
  ## the k4geo optical tables, traced through the cell optics of the k4geo geometry (mirror, aerogel slab, SiPM array,
  ## pixels) and reconstructed as in the full-simulation analysis. Configuration {cfg}: {NOTE[cfg]}; gas {gas}.
  set InputArray TrackSmearing/tracks
  set OutputArray tracks
  set Bz $::B

  set CellFile {HERE}/{CELLS[cfg]}
  set OpticsFile {HERE}/arc_optics_tables.txt

  ## radiators: index tables (n - 1 scaled by *IndexScale). Aerogel attenuation as in the Geant4 material of the full
  ## simulation (k4geo Aerogel_PFRICH, checked with a ddsim DEBUG dump): Rayleigh scattering with the Belle II table, no
  ## absorption (deliberately disabled in materials_arc_o1_v01.xml); scattered photons are propagated, not killed.
  set GasIndexTable {GAS[gas][0]}
  set GasIndexScale 1.0
  set GasAbsorptionTable "{GAS[gas][1]}"
  set AerogelIndexTable RINDEX__Aerogel_PFRICH
  set AerogelIndexScale {AERO_SCALE[cfg]:.6f}
  set AerogelAbsorptionTable ""
  set AerogelRayleighTable RAYLEIGH__Aerogel_BelleII
  set AerogelThickness 0.0        ;# mm, 0: from the geometry (10 mm)
  set AerogelSlabInPath true      ;# the photons cross the aerogel again on the way to the sensor

  ## photon detection and readout
  set PDETable ARC_SiPM_QuantumEfficiency
  set PDEScale 1.0
  set MirrorReflectivityTable REFLECTIVITY__ARC_Mirror
  set YieldScale 1.0
  set PixelPitch 0.8              ;# mm
  set SensorSizeScale 1.0         ;# design knob: SiPM array size relative to the geometry
  set SensorShift 0.0             ;# mm along the cell axis
  set MirrorRadiusScale 1.0       ;# design knob: curvature radius, apex fixed
  set DarkCountsPerSensor 0.0

  ## reconstruction (as the full-simulation analysis): emission at mid-radiator, window around the hypothesis, 3-sigma clipping
  set RecoEnergy 3.0              ;# eV, the indices used by the reconstruction
  set RecoRadiusGas 2000.0
  set RecoRadiusAerogel 1925.0
  set RecoZGas 2110.0
  set RecoZAerogel 2025.0
  set Window 0.025
  set MinPe 3.0
  set GasSigmaSingleFloor 1.0E-3
  set AerogelSigmaSingleFloor 1.5E-3
  set GasSigmaCorr 0.0
  set AerogelSigmaCorr 0.0
  set SigmaTrack 0.0
  set UseTrueTrack {"false" if production else "true"}
  set UseTrueMomentumForPID false
  set TrueHypothesisAngle {"false" if production else "true"}
}}
"""
card = card[:m.start()] + block + card[m.end():]
card = card.replace("ARCPID module) are the first-pass values of the ARCPID", "ARCPID module) replaced by the ray-tracing ARCFastSim module", 1)
open(out, "w").write(card)
print("written", out)
