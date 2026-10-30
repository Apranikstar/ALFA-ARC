// T14: extra FCCAnalyses functions for the ALFA jet-flavour ntuples, JIT-included by the analysis script (includePaths).
//  1. The four JetConstituentsUtils functions added by hand in the user's FCCAnalyses (get_px, get_py, get_pz, mask),
//     same signatures and behaviour, so jetFlavourHelper.py runs unchanged on the stock key4hep FCCAnalyses.
//  2. ARC particle identification per jet constituent: k4SimDelphes (ALFA build) writes one edm4hep::ParticleID per track in
//     EFlowTrack_PID, linked to its ReconstructedParticle, with parameters
//       [P(pi), P(K), P(p), theta_c gas, theta_c aerogel, sigma gas, sigma aerogel, Npe gas, Npe aerogel, mass]   (-1 = no ring)
//     arc_param_per_rp() spreads parameter k onto the ReconstructedParticles (default -1: neutrals, no PID), and
//     arc_gather() picks it for each constituent of each jet through the clustering indices.
#ifndef ALFA_JET_UTILS_H
#define ALFA_JET_UTILS_H

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "podio/ObjectID.h"
#include <vector>

namespace FCCAnalyses {
namespace JetConstituentsUtils {
namespace rv = ROOT::VecOps;

inline rv::RVec<FCCAnalysesJetConstituentsData> get_px(const rv::RVec<FCCAnalysesJetConstituents> &jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jc : jcs) {
    FCCAnalysesJetConstituentsData v;
    for (const auto &p : jc) v.push_back(p.momentum.x);
    out.push_back(v);
  }
  return out;
}

inline rv::RVec<FCCAnalysesJetConstituentsData> get_py(const rv::RVec<FCCAnalysesJetConstituents> &jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jc : jcs) {
    FCCAnalysesJetConstituentsData v;
    for (const auto &p : jc) v.push_back(p.momentum.y);
    out.push_back(v);
  }
  return out;
}

inline rv::RVec<FCCAnalysesJetConstituentsData> get_pz(const rv::RVec<FCCAnalysesJetConstituents> &jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  for (const auto &jc : jcs) {
    FCCAnalysesJetConstituentsData v;
    for (const auto &p : jc) v.push_back(p.momentum.z);
    out.push_back(v);
  }
  return out;
}

inline rv::RVec<rv::RVec<int>> mask(const rv::RVec<FCCAnalysesJetConstituentsData> &energies) {
  rv::RVec<rv::RVec<int>> out;
  for (const auto &e_vec : energies) {
    rv::RVec<int> jet_mask;
    for (const auto &e : e_vec) jet_mask.emplace_back(e != 0.0f ? 1 : 0);
    out.emplace_back(jet_mask);
  }
  return out;
}

}  // namespace JetConstituentsUtils
}  // namespace FCCAnalyses

namespace ALFA {
namespace rv = ROOT::VecOps;

// parameter k of the ARC ParticleID of each ReconstructedParticle (default -1)
inline rv::RVec<float> arc_param_per_rp(const rv::RVec<edm4hep::ParticleIDData> &pids, const rv::RVec<podio::ObjectID> &pid_particle,
                                        const rv::RVec<float> &params, const rv::RVec<edm4hep::ReconstructedParticleData> &rps, int k) {
  rv::RVec<float> out(rps.size(), -1.f);
  for (size_t i = 0; i < pids.size() && i < pid_particle.size(); ++i) {
    const int idx = pid_particle[i].index;
    const unsigned b = pids[i].parameters_begin, e = pids[i].parameters_end;
    if (idx < 0 || idx >= (int)out.size() || b + k >= e) continue;
    out[idx] = params[b + k];
  }
  return out;
}

// value per constituent of each jet (constituent indices from JetClusteringUtils::get_constituents)
inline rv::RVec<rv::RVec<float>> arc_gather(const rv::RVec<float> &per_rp, const std::vector<std::vector<int>> &jet_idx) {
  rv::RVec<rv::RVec<float>> out;
  for (const auto &jet : jet_idx) {
    rv::RVec<float> v;
    for (const int i : jet) v.push_back(i >= 0 && i < (int)per_rp.size() ? per_rp[i] : -1.f);
    out.push_back(v);
  }
  return out;
}

// 1 where a ring was measured (angle > 0), else 0
inline rv::RVec<rv::RVec<float>> arc_has_ring(const rv::RVec<rv::RVec<float>> &angle) {
  rv::RVec<rv::RVec<float>> out;
  for (const auto &jet : angle) {
    rv::RVec<float> v;
    for (const float a : jet) v.push_back(a > 0.f ? 1.f : 0.f);
    out.push_back(v);
  }
  return out;
}

}  // namespace ALFA
#endif
