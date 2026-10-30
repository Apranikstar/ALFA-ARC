"""T06 extractor plus photon truth, for T15: every ARC optical-photon hit with its cell / pixel (as T06), the exact hit
position, and the photon's MCParticle (needs ddsim --part.keepAllParticles): creation vertex (= true emission point),
initial momentum (= true emission direction), creation time, and whether its parent is the primary track.
Usage: python3 t15_extract.py sim.edm4hep.root out.npz"""
import math
import sys

import numpy as np
from podio.root_io import Reader

PITCH_MM, OFFSET_MM, MAX_PHOTON_EV = 0.8, -39.6, 20.0


def signed(v, bits):
    return v - (1 << bits) if v >= (1 << (bits - 1)) else v


ev_cols = {k: [] for k in ("pdg", "px", "py", "pz", "mass", "t0", "hit_begin", "hit_end")}
hit_cols = {k: [] for k in ("system", "barrel", "cell", "x_mm", "y_mm", "t_ns", "e_ev", "pos_x", "pos_y", "pos_z",
                            "vtx_x", "vtx_y", "vtx_z", "dir_x", "dir_y", "dir_z", "t_emit", "from_primary")}
for ev in Reader(sys.argv[1]).get("events"):
    mcs = ev.get("MCParticles")
    prim = mcs[0]
    p = prim.getMomentum()
    for k, v in (("pdg", prim.getPDG()), ("px", p.x), ("py", p.y), ("pz", p.z), ("mass", prim.getMass()), ("t0", prim.getTime())):
        ev_cols[k].append(v)
    ev_cols["hit_begin"].append(len(hit_cols["t_ns"]))
    for h in ev.get("ArcCollection"):
        m = h.getMomentum()
        e_ev = math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2) * 1e9
        if e_ev > MAX_PHOTON_EV:
            continue
        c = h.getCellID()
        hit_cols["system"].append(c & 0x1F); hit_cols["barrel"].append((c >> 5) & 0x7); hit_cols["cell"].append((c >> 8) & 0xFFFFFF)
        hit_cols["x_mm"].append(OFFSET_MM + signed((c >> 32) & 0xFFFF, 16) * PITCH_MM)
        hit_cols["y_mm"].append(OFFSET_MM + signed((c >> 48) & 0xFFFF, 16) * PITCH_MM)
        hit_cols["t_ns"].append(h.getTime()); hit_cols["e_ev"].append(e_ev)
        x = h.getPosition(); hit_cols["pos_x"].append(x.x); hit_cols["pos_y"].append(x.y); hit_cols["pos_z"].append(x.z)
        ph = h.getParticle()
        if ph.isAvailable() and ph.getPDG() in (-22, 0, 20022):  # Geant4 optical photon: PDG -22 (DD4hep) or 0
            v = ph.getVertex(); q = ph.getMomentum(); n = math.sqrt(q.x ** 2 + q.y ** 2 + q.z ** 2) or 1.0
            par = ph.getParents()
            fp = 1 if (len(par) and par[0].getObjectID().index == prim.getObjectID().index) else 0
            vals = (v.x, v.y, v.z, q.x / n, q.y / n, q.z / n, ph.getTime(), fp)
        else:
            vals = (np.nan,) * 7 + (-1,)
        for k, val in zip(("vtx_x", "vtx_y", "vtx_z", "dir_x", "dir_y", "dir_z", "t_emit", "from_primary"), vals):
            hit_cols[k].append(val)
    ev_cols["hit_end"].append(len(hit_cols["t_ns"]))

out = {"ev_" + k: np.array(v) for k, v in ev_cols.items()}
I = ("system", "barrel", "cell", "from_primary")
out.update({"hit_" + k: np.array(v, dtype=np.int32 if k in I else np.float64 if k.startswith(("vtx", "pos", "dir")) else np.float32)
            for k, v in hit_cols.items()})
np.savez_compressed(sys.argv[2], **out)
fp = np.array(hit_cols["from_primary"])
print("events", len(ev_cols["pdg"]), "photon hits", len(fp), "with photon truth", int((fp >= 0).sum()), "from primary", int((fp == 1).sum()))
