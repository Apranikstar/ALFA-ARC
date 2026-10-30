"""
ddsim steering for ALFA_o1_v01: outer tracker + ARC (optical photons) + GRAiNITA ECAL in one run.
Tracker and ECAL settings follow SteeringFile_ALFA_o1_v00.py; Cerenkov/optical physics and the ARC
optical tracker action follow arcfullsim.py. Output is EDM4hep.

    ddsim --steeringFile SteeringFile_ALFA_o1_v01_ARC.py \
          --compactFile $K4GEO/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml \
          --enableGun --gun.particle pi+ --gun.energy "40*GeV" \
          --gun.thetaMin "90*deg" --gun.thetaMax "90*deg" --gun.distribution uniform -N 10

Set ARC_PDE_IN_SIM=1 in the environment to apply the SiPM photon detection efficiency in the simulation
(ARCPhotonDetectionFilter on the ARC sensors); by default every photon reaching a sensor makes a hit and the
PDE is applied in the analysis.
"""

import os

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import GeV

SIM = DD4hepSimulation()

SIM.compactFile = ["ALFA_o1_v01.xml"]
SIM.runType = "batch"
SIM.numberOfEvents = 10
SIM.outputFile = "alfa_arc_sim.edm4hep.root"
SIM.outputConfig.forceDD4HEP = False

SIM.enableGun = True
SIM.gun.particle = "pi+"
SIM.gun.energy = 40 * GeV  # the DDSim gun sets this as the momentum
SIM.gun.distribution = "uniform"
SIM.gun.multiplicity = 1
SIM.gun.position = "0 0 0*cm"


# Cerenkov and optical photon physics (as in arcfullsim.py)
def setupCerenkov(kernel):
    from DDG4 import PhysicsList

    seq = kernel.physicsList()
    cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    cerenkov.MaxNumPhotonsPerStep = 10
    cerenkov.MaxBetaChangePerStep = 10.0
    cerenkov.TrackSecondariesFirst = False
    cerenkov.VerboseLevel = 0
    cerenkov.enableUI()
    seq.adopt(cerenkov)
    ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
    ph.addParticleConstructor("G4OpticalPhoton")
    ph.VerboseLevel = 0
    # The ARC sensors carry a placeholder RINDEX so photons enter them; see arcfullsim.py
    ph.BoundaryInvokeSD = False
    ph.enableUI()
    seq.adopt(ph)
    return None


SIM.physics.setupUserPhysics(setupCerenkov)
SIM.physics.list = "FTFP_BERT"
SIM.physics.rangecut = 0.7

# ARC: DD4hep optical tracker, no energy filter; optionally the SiPM PDE as a filter
if os.environ.get("ARC_PDE_IN_SIM", "0") == "1":
    SIM.filter.filters["arcpde"] = {"name": "ARCPhotonDetectionFilter/ARCPDE", "parameter": {}}
    SIM.filter.mapDetFilter["ARCBARREL"] = "arcpde"
    SIM.filter.mapDetFilter["ARCENDCAP"] = "arcpde"
else:
    SIM.filter.mapDetFilter["ARCBARREL"] = None
    SIM.filter.mapDetFilter["ARCENDCAP"] = None

# Tracker and ECAL as in SteeringFile_ALFA_o1_v00.py
SIM.filter.tracker = "edep1kev"
SIM.action.tracker = (
    "Geant4TrackerWeightedAction",
    {"HitPositionCombination": 2, "CollectSingleDeposits": False},
)
SIM.action.calo = "GrainitaCaloSDAction"
SIM.action.mapActions = {
    "ARCBARREL": "Geant4OpticalTrackerAction",
    "ARCENDCAP": "Geant4OpticalTrackerAction",
    "Grainita_ECAL_Barrel_v02": (
        "GrainitaCaloSDAction",
        {
            "ReadoutName": "GrainitaEcalBarrelRO",
            "CollectionName": "GrainitaCalorimeterHits",
            "useLightResponseFunction": True,
            "responseFuncSlope": 0.93,
            "responseFuncIntersect": 0.206,
            "responseFuncX0": 0.856,
            "responseFuncAttLength": 3.4,
            "neighborRadius": 3,
            "fiberAttenuationLength": 1e7,
            "outerRadius": 2645,
        },
    ),
}

SIM.field.delta_chord = 0.25
SIM.field.equation = "Mag_UsualEqRhs"
SIM.field.stepper = "ClassicalRK4"

SIM.part.minimalKineticEnergy = 1.0 * GeV
SIM.part.saveProcesses = ["Decay"]
SIM.part.userParticleHandler = ""
SIM.physics.zeroTimePDGs = {17, 11, 13, 15}
