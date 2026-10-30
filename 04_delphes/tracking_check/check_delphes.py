"""Track resolution in Delphes output (pi+ gun): RMS of D0, DZ and of (PT - true PT)/PT, and the mean reported errors, for
40 < p < 50 GeV, 62.5 < theta < 117.5 deg (the T11 sanity selection). Usage: python3 check_delphes.py file.root [...]"""
import sys
import numpy as np
import uproot
for f in sys.argv[1:]:
    t = uproot.open(f)["Delphes"]
    a = t.arrays(["Track.PT", "Track.P", "Track.CtgTheta", "Track.D0", "Track.DZ", "Track.ErrorD0", "Track.ErrorDZ", "Track.ErrorPT",
                  "Particle.PT", "Particle.Status", "Particle.PID"], library="ak")
    import awkward as ak
    p = ak.flatten(a["Track.P"]).to_numpy(); ct = ak.flatten(a["Track.CtgTheta"]).to_numpy()
    th = np.degrees(np.arctan2(1, ct))
    s = (p > 40) & (p < 50) & (th > 62.5) & (th < 117.5)
    d0 = ak.flatten(a["Track.D0"]).to_numpy()[s] * 1e3; dz = ak.flatten(a["Track.DZ"]).to_numpy()[s] * 1e3   # mm -> um
    e0 = ak.flatten(a["Track.ErrorD0"]).to_numpy()[s] * 1e3; ez = ak.flatten(a["Track.ErrorDZ"]).to_numpy()[s] * 1e3
    pt = ak.flatten(a["Track.PT"]).to_numpy()[s]; ept = ak.flatten(a["Track.ErrorPT"]).to_numpy()[s]
    ntr = len(ak.flatten(a["Track.P"])); nev = len(a)
    print(f"{f.split('/')[-1]}: {nev} events, {ntr} tracks; selection {s.sum()}: rms D0 {np.std(d0):.2f} um (mean ErrorD0 {np.mean(e0):.2f}), "
          f"rms DZ {np.std(dz):.2f} um (mean ErrorDZ {np.mean(ez):.2f}), mean ErrorPT/PT {np.mean(ept / pt):.5f}")
