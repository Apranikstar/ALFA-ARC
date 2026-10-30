// T16: track-parameter resolution of a Delphes TrackCovariance geometry, computed with the same SolTrack::KalmanCovT
// call that ObsTrk makes inside Delphes (Res and MS on, pion mass). No smearing statistics: sigma = sqrt(cov diagonal).
//   mode "grid":  sigma(d0), sigma(z0) [um], sigma(pT)/pT, sigma(theta) [mrad] on a (p, theta) grid
//   mode "gun":   N tracks, log-uniform p in [pmin, pmax], uniform |eta| < etamax (as the T11 guns); one line per track
//                 with the KalmanCovT success flag, so the 'negative diagonal' tracks can be located
// Usage: root -l -b -q 'trkres.C+("geom.txt","grid","out.txt")'
#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolTrack.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMath.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdlib>

struct Res { bool ok; int nmeas; double d0, z0, dptpt, dth; };

static Res one(SolGeom *G, double p, double th, double phi, double bz)
{
  TVector3 x(0, 0, 0), mom;
  mom.SetMagThetaPhi(p, th, phi);
  SolTrack trk(x, mom, 1.0, G);
  // method from TRK_METHOD: kalmanT (default, as ObsTrk in Delphes), kalman (partial Kalman, the Delphes fallback),
  // covcalc (global least-squares fit used by SolGridCov)
  static const std::string meth = getenv("TRK_METHOD") ? getenv("TRK_METHOD") : "kalmanT";
  Res r;
  if (meth == "covcalc") { trk.CovCalc(kTRUE, kTRUE); r.ok = true; }
  else if (meth == "kalman") { trk.KalmanCov(kTRUE, kTRUE, 0.13957021); r.ok = trk.fKsuccess; }
  else { trk.KalmanCovT(kTRUE, kTRUE, 0.13957021); r.ok = trk.fKTsuccess; }
  r.nmeas = trk.GetUmeas();
  TMatrixDSym c = trk.Cov();
  double pt = p * TMath::Sin(th), ct = 1. / TMath::Tan(th);
  auto s = [&](int i) { return c(i, i) > 0 ? TMath::Sqrt(c(i, i)) : -1.; };
  r.d0 = s(0) * 1e6;                                      // m -> um
  r.z0 = s(3) * 1e6;
  r.dptpt = s(2) > 0 ? 2 * s(2) * pt / (0.2998 * bz) : -1;  // C = half curvature [1/m]
  r.dth = s(4) > 0 ? s(4) / (1 + ct * ct) * 1e3 : -1;      // mrad
  return r;
}

void trkres(const char *geo, const char *mode, const char *out, int N = 10000, double pmin = 0.1, double pmax = 60,
            double etamax = 2.6, int seed = 1)
{
  std::ifstream f(geo);
  std::string first, rest, line;
  double bz = 0;
  std::getline(f, line);
  std::istringstream(line) >> first >> bz;
  std::stringstream body;
  while (std::getline(f, line)) body << line << "\n";
  SolGeom *G = new SolGeom();
  G->Read(body.str().c_str());
  G->SetBz(bz);
  std::ofstream o(out);
  if (std::string(mode) == "grid") {
    double ps[] = {1, 2, 5, 10, 20, 45, 46, 100, 182};
    o << "# p_GeV theta_deg ok nmeas sd0_um sz0_um dpt_pt dtheta_mrad\n";
    for (double p : ps)
      for (double thd = 10; thd <= 90.01; thd += 1) {
        Res r = one(G, p, thd * TMath::DegToRad(), 0.3, bz);
        o << p << " " << thd << " " << r.ok << " " << r.nmeas << " " << r.d0 << " " << r.z0 << " " << r.dptpt << " " << r.dth << "\n";
      }
    // momentum scan at 90 deg (and 45 deg) for smooth curves
    for (double thd : {90., 45.})
      for (int i = 0; i <= 200; ++i) {
        double p = TMath::Power(10, i / 100.);  // 1..100 GeV
        Res r = one(G, p, thd * TMath::DegToRad(), 0.3, bz);
        o << "S " << p << " " << thd << " " << r.ok << " " << r.nmeas << " " << r.d0 << " " << r.z0 << " " << r.dptpt << " " << r.dth << "\n";
      }
  } else {
    TRandom3 rnd(seed);
    o << "# p_GeV theta_deg phi ok nmeas sd0_um sz0_um dpt_pt dtheta_mrad\n";
    for (int i = 0; i < N; ++i) {
      double p = pmin * TMath::Power(pmax / pmin, rnd.Rndm());
      double eta = (2 * rnd.Rndm() - 1) * etamax, th = 2 * TMath::ATan(TMath::Exp(-eta)), phi = rnd.Rndm() * TMath::TwoPi();
      Res r = one(G, p, th, phi, bz);
      o << p << " " << th * TMath::RadToDeg() << " " << phi << " " << r.ok << " " << r.nmeas << " " << r.d0 << " " << r.z0 << " "
        << r.dptpt << " " << r.dth << "\n";
    }
  }
  std::cout << "wrote " << out << std::endl;
}
