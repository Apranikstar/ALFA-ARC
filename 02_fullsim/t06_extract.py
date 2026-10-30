"""Extract what T06 needs from one full-ALFA EDM4hep file into a small .npz: the primary per event and every
ARC optical-photon hit (cell, sensor x/y in mm, time, photon energy). Charged hits (E > 20 eV) are dropped.
cellID layout: system:5,barrel:3,cellnumber:24,x:32:-16,y:-16; pixel 0.8 mm, offset -39.6 mm (CartesianGridXY).
Usage: python3 t06_extract.py sim.edm4hep.root out.npz"""
import math
import sys

import numpy as np
from podio.root_io import Reader

PITCH_MM, OFFSET_MM, MAX_PHOTON_EV = 0.8, -39.6, 20.0


def signed(v, bits):
    return v - (1 << bits) if v >= (1 << (bits - 1)) else v


ev_cols = {k: [] for k in ("pdg", "px", "py", "pz", "mass", "t0", "hit_begin", "hit_end")}
hit_cols = {k: [] for k in ("system", "barrel", "cell", "x_mm", "y_mm", "t_ns", "e_ev")}
for ev in Reader(sys.argv[1]).get("events"):
    prim = ev.get("MCParticles")[0]
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
    ev_cols["hit_end"].append(len(hit_cols["t_ns"]))

out = {"ev_" + k: np.array(v) for k, v in ev_cols.items()}
out.update({"hit_" + k: np.array(v, dtype=np.float32 if k in ("x_mm", "y_mm", "t_ns", "e_ev") else np.int32) for k, v in hit_cols.items()})
np.savez_compressed(sys.argv[2], **out)
print("events", len(ev_cols["pdg"]), "photon hits", len(hit_cols["t_ns"]))
