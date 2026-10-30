"""[v2, T15: the pixel centre is put on the sensor entrance face (PIX_DEPTH), not on the sensor mid-plane]
Per-photon Cherenkov-angle reconstruction for the ALFA_o1_v01 ARC (xenon 2 bar), from the T06 .npz files.

For every detected photon:
  - detection point D: the pixel centre from the cellID, put in the world frame with the DD4hep VolumeManager
    (checked against stored SimTrackerHit positions: 0.32 mm mean, 1.06 mm worst);
  - mirror: the sphere that cuts the cell's mirror (centre C, reflecting radius R = inner radius, gas side), from the
    mirror sibling of the sensor node in the TGeo geometry (navigated from the sensor centre);
  - track: a helix from the vertex in a uniform 3 T field along z (checked: -38.9 mrad measured vs -38.9 predicted at
    the outer tracker for 40 GeV pi+), evaluated at the emission point = mid-depth of the gas (barrel r = 2000 mm,
    endcap |z| = 2110 mm) or of the aerogel (barrel r = 1925 mm, endcap |z| = 2025 mm);
  - reflection point M on the sphere: it lies in the plane through C, E, D; solve (a + b).t = 0 for the angle along
    the sphere, with a, b the unit vectors M->E, M->D and t the tangent (vectorised bisection);
  - photon direction k = (M - E)/|M - E| in the gas; for the aerogel hypothesis k is refracted back into the aerogel
    through its exit surface (normal: radial in the barrel, z in the endcap), n_gas -> n_aero;
  - theta_c = angle(k, track direction at E).
Per track and radiator: photons are assigned to the radiator whose hypothesis is nearer the expected angle (truth
seeded with the true species and momentum, as in arc_ring_analysis.py), within 25 mrad; then the PDE-weighted mean
with two passes of 3-sigma rejection around the median.
Units: mm, GeV, mrad."""
import math

import numpy as np
import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kError

B_T = 3.0
R_GAS_MID, R_AERO_MID = 2000.0, 1925.0
Z_GAS_MID, Z_AERO_MID = 2110.0, 2025.0
C_MM_NS = 299.792458
QE_E = np.array([1.1, 1.3, 1.550, 1.771, 2.101, 2.465, 2.655, 2.743, 2.837, 2.938, 3.002, 3.084, 3.163, 3.378, 3.647, 3.850, 4.147, 4.320, 4.5, 4.7, 4.9])
QE_V = np.array([0.0, 0.0, 0.13, 0.26, 0.4, 0.51, 0.53, 0.56, 0.57, 0.59, 0.6, 0.6, 0.59, 0.52, 0.47, 0.46, 0.41, 0.2, 0.1, 0.0, 0.0])
PIX_DEPTH_MM = 1.0  # T15: photons are detected on the sensor face; the cellID pixel centre is on the mid-plane of the 2 mm sensor
COMPACT = "/afs/cern.ch/work/h/hfatehi/alfa/k4geo/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml"


class ArcGeometry:
    def __init__(self, compact=COMPACT):
        ROOT.gSystem.Load("libDDCore")
        self.desc = ROOT.dd4hep.Detector.getInstance()
        self.desc.fromXML(compact)
        self.vm = self.desc.volumeManager()
        if not self.vm.isValid():  # standalone compact files do not load the plugin; ALFA_o1_v01 does
            self.desc.apply("DD4hepVolumeManager", 0, ROOT.nullptr)
            self.vm = self.desc.volumeManager()
        self.mgr = self.desc.manager()
        self.cache = {}

    def cell(self, system, barrel, cell):
        key = (int(system), int(barrel), int(cell))
        if key in self.cache:
            return self.cache[key]
        vid = key[0] | (key[1] << 5) | (key[2] << 8)
        out = None
        try:
            ctx = self.vm.lookupContext(vid)
            o = ctx.localToWorld(ROOT.dd4hep.Position(0, 0, 0))
            ex = ctx.localToWorld(ROOT.dd4hep.Position(1, 0, 0))
            ey = ctx.localToWorld(ROOT.dd4hep.Position(0, 1, 0))
            O = np.array([o.X(), o.Y(), o.Z()]) * 10.0
            # mirror: navigate from the sensor centre to the sensor node, up to its cell, then the mirror daughter.
            # (The endcap VolumeManager entries point to the whole endcap DetElement, so the DetElement tree cannot
            # be used there; the TGeo navigation works for barrel and endcap alike.)
            nav = self.mgr.GetCurrentNavigator()
            node = nav.FindNode(o.X(), o.Y(), o.Z())
            C = R = None
            if node and "sensor" in node.GetName():
                nav.CdUp()
                m_cell = ROOT.TGeoHMatrix(nav.GetCurrentMatrix())
                cellnode = nav.GetCurrentNode()
                for k in range(cellnode.GetNdaughters()):
                    dn = cellnode.GetDaughter(k)
                    if "mirror" in dn.GetName():
                        m_mir = ROOT.TGeoHMatrix(m_cell); m_mir.Multiply(dn.GetMatrix())
                        bn = dn.GetVolume().GetShape().GetBoolNode()
                        t = bn.GetRightMatrix().GetTranslation()
                        loc = np.array([t[0], t[1], t[2]], dtype=float); glob = np.zeros(3)
                        m_mir.LocalToMaster(loc, glob)
                        C, R = glob * 10.0, bn.GetRightShape().GetRmin() * 10.0
                        break
            if C is not None:
                out = dict(C=C, R=R, O=O,
                           EX=(np.array([ex.X(), ex.Y(), ex.Z()]) * 10.0 - O) / 10.0,
                           EY=(np.array([ey.X(), ey.Y(), ey.Z()]) * 10.0 - O) / 10.0)
        except Exception as exc:
            if not getattr(self, "_warned", False):
                print("ArcGeometry: cell lookup failed for", key, "-", exc)
                self._warned = True
            out = None
        self.cache[key] = out
        return out

    @staticmethod
    def to_world(g, x_mm, y_mm):
        """Sensor-local pixel centres (mm, segmentation frame) -> world (mm); the frame is linear, so 3 points suffice."""
        return g["O"][None, :] + np.outer(x_mm, g["EX"]) + np.outer(y_mm, g["EY"])


def face_normal(g):
    """Unit normal of the sensor pointing to the face the photons hit: the side facing the mirror, i.e. away from the mirror
    centre C. EX x EY points the other way in the barrel and this way in the endcaps (checked with photon truth, T15)."""
    nz = np.cross(g["EX"], g["EY"]); nz /= np.linalg.norm(nz)
    return nz if nz @ (g["O"] - g["C"]) > 0 else -nz


def helix_point(p, q, target, b_t=None):
    """Point and unit direction of a helix from the origin with momentum p (GeV, 3-vector) and charge q in B_T along z,
    where it reaches the barrel radius target['r'] or the endcap plane |z| = target['z'] (whichever first)."""
    px, py, pz = p
    pt = math.hypot(px, py)
    phi0 = math.atan2(py, px)
    bt = B_T if b_t is None else b_t
    if bt == 0:  # straight line (standalone ARC, no field)
        pv = np.array(p, dtype=float); u = pv / np.linalg.norm(pv)
        s_r = target["r"] / math.hypot(u[0], u[1]) if math.hypot(u[0], u[1]) > 0 else float("inf")
        s_z = target["z"] / abs(u[2]) if u[2] != 0 else float("inf")
        L = min(s_r, s_z)
        return u * L, u, ("barrel" if s_r <= s_z else "endcap")
    rc = pt / (0.3 * bt) * 1000.0  # mm
    # transverse arc length to radius r: chord 2 rc sin(s/(2 rc)) = r
    s_r = 2 * rc * math.asin(min(1.0, target["r"] / (2 * rc))) if target["r"] < 2 * rc else float("inf")
    s_z = target["z"] / abs(pz) * pt if pz != 0 else float("inf")
    s = min(s_r, s_z)
    psi = -q * s / rc
    x = (rc / -q) * (math.sin(phi0 + psi) - math.sin(phi0))
    y = (rc / -q) * (-math.cos(phi0 + psi) + math.cos(phi0))
    z = s * pz / pt
    d = np.array([pt * math.cos(phi0 + psi), pt * math.sin(phi0 + psi), pz])
    return np.array([x, y, z]), d / np.linalg.norm(d), ("barrel" if s_r <= s_z else "endcap")


def reflect_angle(E, D, C, R):
    """Vectorised over photons D (N x 3): mirror point M and photon direction k from emission point E."""
    e1 = (E - C) / np.linalg.norm(E - C)
    w = (D - C) - np.outer((D - C) @ e1, e1)
    e2 = w / np.linalg.norm(w, axis=1)[:, None]
    beta_d = np.arctan2(np.einsum("ij,ij->i", D - C, e2), (D - C) @ e1)

    def f(al):
        M = C + R * (np.cos(al)[:, None] * e1[None, :] + np.sin(al)[:, None] * e2)
        tvec = -np.sin(al)[:, None] * e1[None, :] + np.cos(al)[:, None] * e2
        a = (E[None, :] - M); a /= np.linalg.norm(a, axis=1)[:, None]
        b = (D - M); b /= np.linalg.norm(b, axis=1)[:, None]
        return np.einsum("ij,ij->i", a + b, tvec), M

    lo, hi = np.zeros(len(D)), beta_d.copy()
    flo, _ = f(lo)
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        fm, _ = f(mid)
        same = np.sign(fm) == np.sign(flo)
        lo = np.where(same, mid, lo); flo = np.where(same, fm, flo); hi = np.where(same, hi, mid)
    _, M = f(0.5 * (lo + hi))
    k = M - E[None, :]
    return M, k / np.linalg.norm(k, axis=1)[:, None]


def refract_into(k, normal, n1, n2):
    """Refract unit vectors k (N x 3) leaving medium n2 into n1 back through a surface with unit normal 'normal':
    returns the direction inside medium n2 (inverse Snell, n1 sin t1 = n2 sin t2)."""
    nrm = normal / np.linalg.norm(normal)
    cos1 = k @ nrm
    sgn = np.sign(cos1); sgn[sgn == 0] = 1
    tang = k - np.outer(cos1, nrm)
    sin1 = np.linalg.norm(tang, axis=1)
    sin2 = np.clip(n1 / n2 * sin1, 0, 1)
    tdir = tang / np.maximum(sin1, 1e-12)[:, None]
    return tdir * sin2[:, None] + np.outer(sgn * np.sqrt(1 - sin2 ** 2), nrm)


def event_angles(d, i, geo, n_gas, n_aero, gate=(0.4, 1.4), b_t=None):
    """Per-photon gas and aerogel hypotheses (mrad) and PDE weights for event i of an npz dict (main cell, TOF-gated)."""
    b, e = d["ev_hit_begin"][i], d["ev_hit_end"][i]
    if e <= b:
        return None
    p = np.array([d["ev_px"][i], d["ev_py"][i], d["ev_pz"][i]])
    q = 1 if d["ev_pdg"][i] > 0 else -1
    m = float(d["ev_mass"][i]); P = np.linalg.norm(p); beta = P / math.sqrt(P * P + m * m)
    cells = d["hit_cell"][b:e] + 1000000 * d["hit_barrel"][b:e]
    vals, cnt = np.unique(cells, return_counts=True)
    main = vals[np.argmax(cnt)]
    sel = cells == main
    # vessel arrival time as in T02/T04 (straight-line path length is adequate for the gate)
    sin_t, cos_t = math.hypot(p[0], p[1]) / P, abs(p[2]) / P
    L = min(1910.0 / sin_t if sin_t > 0 else 1e12, 2010.0 / cos_t if cos_t > 0 else 1e12)
    tof = d["hit_t_ns"][b:e] - (d["ev_t0"][i] + L / (beta * C_MM_NS))
    sel &= (tof >= gate[0]) & (tof < gate[1])
    if sel.sum() < 1:
        return None
    j = np.flatnonzero(sel)
    g = geo.cell(d["hit_system"][b + j[0]], d["hit_barrel"][b + j[0]], d["hit_cell"][b + j[0]])
    if g is None:
        return None
    D = geo.to_world(g, d["hit_x_mm"][b:e][sel].astype(float), d["hit_y_mm"][b:e][sel].astype(float))
    D = D + PIX_DEPTH_MM * face_normal(g)[None, :]
    w = np.interp(d["hit_e_ev"][b:e][sel], QE_E, QE_V, left=0.0, right=0.0)
    if "hit_isnoise" in d:
        w = np.where(d["hit_isnoise"][b:e][sel] > 0, 1.0, w)
    out = {"w": w}
    for rad, tgt in (("gas", {"r": R_GAS_MID, "z": Z_GAS_MID}), ("aerogel", {"r": R_AERO_MID, "z": Z_AERO_MID})):
        E, t, region = helix_point(p, q, tgt, b_t)
        _, k = reflect_angle(E, D, g["C"], g["R"])
        if rad == "aerogel":
            normal = np.array([E[0], E[1], 0.0]) if region == "barrel" else np.array([0.0, 0.0, np.sign(E[2])])
            k = refract_into(k, normal, n_gas, n_aero)
        out[rad] = np.degrees(np.arccos(np.clip(k @ t, -1, 1))) * (math.pi / 180) * 1e3
    return out


def track_angle(theta, w, expected, window=25.0):
    """PDE-weighted mean of the photons near the expected angle, two passes of 3-sigma rejection around the median."""
    m = np.abs(theta - expected) < window
    if w[m].sum() <= 0:
        return None
    for _ in range(2):
        med = np.median(theta[m])
        s = 1.4826 * np.median(np.abs(theta[m] - med)) if m.sum() > 3 else window
        m = m & (np.abs(theta - med) < max(3 * s, 2.0))
        if w[m].sum() <= 0:
            return None
    return float(np.average(theta[m], weights=w[m])), float(w[m].sum()), int(m.sum()), theta[m], w[m]
