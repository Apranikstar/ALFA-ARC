# T14: FCCAnalyses ntupler for ALFA jet-flavour tagging, e+e- -> ZH (Z -> nu nu, H -> bb/cc/ss/uu/dd/gg) at 240 GeV.
# Based on the user's analysis (JetFlavourHelper as given), with: FCCAnalyses' ExclusiveJetClusteringHelper (exclusive Durham /
# ee_kt clustering of ALL reconstructed particles into exactly 2 jets),
# no ONNX inference (no model yet), and the ARC particle-ID of every constituent (P_pi, P_K, P_p, ring flags).
# Usage: ALFA_CARD=Xe2bar fccanalysis run analysis_jets.py      (or C4F10)
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jetFlavourHelper import JetFlavourHelper  # noqa: E402
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper  # noqa: E402

CARD = os.environ.get("ALFA_CARD", "Xe2bar")
BASE = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2"  # campaign_v2 (T17)

processList = {f"H{d}": {"fraction": 1, "chunks": 25} for d in ("bb", "cc", "ss", "uu", "dd", "gg")}
inputDir = f"{BASE}/edm4hep/{CARD}"
outputDir = f"{BASE}/ntuples/{CARD}"
nCPUS = int(os.environ.get("ALFA_NCPUS", "16"))
if os.environ.get("ALFA_TEST"):  # quick test: one process, output to ALFA_TEST
    processList = {"Hbb": {"fraction": 1}}
    outputDir = os.environ["ALFA_TEST"]
includePaths = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "alfa_jet_utils.h")]

jetFlavourHelper = None
jetClusteringHelper = None

# ARC ParticleID parameters (k4SimDelphes ALFA build): index -> constituent feature name
ARC_PARAMS = {0: "pfcand_arc_probPi", 1: "pfcand_arc_probK", 2: "pfcand_arc_probP", 3: "pfcand_arc_thetaGas", 4: "pfcand_arc_thetaAero"}


class RDFanalysis:
    def analysers(df):
        global jetFlavourHelper, jetClusteringHelper

        coll = {
            "GenParticles": "Particle",
            "PFParticles": "ReconstructedParticles",
            "PFTracks": "EFlowTrack",
            "PFPhotons": "EFlowPhoton",
            "PFNeutralHadrons": "EFlowNeutralHadron",
            "TrackState": "_EFlowTrack_trackStates",
            "TrackerHits": "TrackerHits",
            "CalorimeterHits": "CalorimeterHits",
            "PathLength": "EFlowTrack_L",
            "Bz": "magFieldBz",
        }

        # all reconstructed particles clustered into exactly 2 jets (exclusive Durham / ee_kt, E-scheme, E-ordered)
        jetClusteringHelper = ExclusiveJetClusteringHelper(coll["PFParticles"], 2)
        df = jetClusteringHelper.define(df)
        df = df.Filter("event_njet > 1")

        jetFlavourHelper = JetFlavourHelper(coll, jetClusteringHelper.jets, jetClusteringHelper.constituents)
        df = jetFlavourHelper.define(df)

        # extra jet / event observables (the clustering helper gives jet_p, e, mass, phi, theta, nconst and event_njet)
        df = df.Define("jet_p4", f"JetConstituentsUtils::compute_tlv_jets({jetClusteringHelper.jets})")
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
        df = df.Define("jet_pT", f"JetClusteringUtils::get_pt({jetClusteringHelper.jets})")
        df = df.Define("jet_eta", f"JetClusteringUtils::get_eta({jetClusteringHelper.jets})")

        # ARC particle-ID per constituent (-1 for neutrals and tracks without a ring)
        for k, name in ARC_PARAMS.items():
            df = df.Define(f"_rp_arc{k}", f"ALFA::arc_param_per_rp(EFlowTrack_PID, _EFlowTrack_PID_particle, _EFlowTrack_PID_parameters, ReconstructedParticles, {k})")
            df = df.Define(name, f"ALFA::arc_gather(_rp_arc{k}, _jetc)")
        df = df.Define("pfcand_arc_hasGasRing", "ALFA::arc_has_ring(pfcand_arc_thetaGas)")
        df = df.Define("pfcand_arc_hasAeroRing", "ALFA::arc_has_ring(pfcand_arc_thetaAero)")
        return df

    def output():
        branchList = jetClusteringHelper.outputBranches() + ["jet_pT", "jet_eta", "event_invariant_mass"]
        branchList += [k for k in jetFlavourHelper.definition if k.startswith("pfcand_") or k == "pfmask"]
        branchList += [k for k in jetFlavourHelper.definition if k.startswith("jet_")]
        branchList += list(ARC_PARAMS.values()) + ["pfcand_arc_hasGasRing", "pfcand_arc_hasAeroRing"]
        return branchList
