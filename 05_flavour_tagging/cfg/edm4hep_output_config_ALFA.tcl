# k4SimDelphes EDM4hep output for the ALFA cards (branch names of delphes_card_ALFA_*.tcl; the card has no ScalarHT)
module EDM4HepOutput EDM4HepOutput {
    add ReconstructedParticleCollections EFlowTrack EFlowPhoton EFlowNeutralHadron
    add GenParticleCollections           Particle
    add JetCollections                   Jet
    add MuonCollections                  Muon
    add ElectronCollections              Electron
    add PhotonCollections                Photon
    add MissingETCollections             MissingET
    set RecoParticleCollectionName       ReconstructedParticles
    set RecoMCParticleLinkCollectionName  RecoMCLink
}
