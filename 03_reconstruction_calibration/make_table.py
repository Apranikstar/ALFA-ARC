"""[T15 copy: corrected reconstruction arc_photon_reco_v2 (pixel centre on the sensor face); reads the T06b/T06c npz, writes
to T06B_OUT (default work/T15/t06b_v2) and keeps per-track photon angles for the bias and likelihood studies (raw_*.pkl)]
T06b table: cell-averaged ARC response from the theta-uniform full sim (condor/T06b, theta 12-90 deg, uniform phi).
Per-photon reconstruction as T06 (arc_photon_reco.py) with PDE SAMPLING (each photon kept with probability PDE, kept photons
weight 1, as a SiPM does). Tracks are binned in theta so that each bin spans whole cells:
  barrel: one bin per z-cell (pitch 25.7 cm; edges from the geometry at r = 2.0 m): theta 90 / 82.6 / 75.4 / 68.9 / 62.7 / 57.2 / 52.4 / 48.0 / 44.2
  gap:    44.2 - 39.3 deg (last barrel cell edge to the outer edge of the endcap cells at r = 1.7 m, z = 2.1 m)
  endcap: ~one hexagonal cell in r (0.25 m) at z = 2.1 m: theta 39.3 / 36.4 / 31.8 / 26.6 / 20.9 / 14.7 / 12.0
Output arc_table_t06b.csv: one row per (species, p, theta bin, radiator) with the T06 columns plus theta_lo/theta_hi/theta_mean
and region.
Usage: python3 t06b_table.py worker <k> <n> <out.pkl>   (files k, k+n, ... of the size-sorted list; one geometry per process)
       python3 t06b_table.py merge <part_*.pkl>         (writes arc_table_t06b.csv)
Independent processes instead of multiprocessing: DD4hep aborts (malloc corruption) inside pool workers."""
import csv, glob, math, os, pickle, re, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # arc_photon_reco.py, arc_ring_analysis.py
import arc_ring_analysis as ra  # noqa: E402

EDGES = [12.0, 14.7, 20.9, 26.6, 31.8, 36.4, 39.3, 44.2, 48.0, 52.4, 57.2, 62.7, 68.9, 75.4, 82.6, 90.0]
REGION = ["endcap"] * 6 + ["gap"] + ["barrel"] * 8
LABEL = {"pip": "pi+", "kaonp": "K+", "proton": "p"}
# T06B_DIR selects the sample (default T06b, Xe 2 bar; T06c is the C4F10 variant), T06B_NGAS the gas index at the detected spectrum
DIR = os.environ.get("T06B_DIR", "/eos/user/h/hfatehi/alfa/work/T06b")
NPZ = f"{DIR}/npz"
OUT = os.environ["T06B_OUT"]
N_GAS = float(os.environ.get("T06B_NGAS", ra.N_GAS))


GEO = None


def expected(sp, p):
    m = ra.MASS[sp]; beta = p / math.sqrt(p * p + m * m)
    g = ra.theta_out(N_GAS, beta, False); c = 1.0 / (ra.N_AERO * beta)
    return {"gas": g * 1e3 if g else None, "aerogel": math.acos(c) * 1e3 if c < 1 else None}


def run_file(f):
    import arc_photon_reco as pr
    global GEO
    m = re.match(r"(pip|kaonp|proton)_p([\d.]+)_c(\d+)\.npz", os.path.basename(f))
    sp, p, ch = m.group(1), float(m.group(2)), int(m.group(3))
    if GEO is None:
        GEO = pr.ArcGeometry()
    geo = GEO; rng = np.random.default_rng(31 * ch + int(10 * p) + {"pip": 0, "kaonp": 1000, "proton": 2000}[sp])
    d = ra.load(f); exp = expected(sp, p)
    th_ev = np.degrees(np.arctan2(np.hypot(d["ev_px"], d["ev_py"]), d["ev_pz"]))
    acc = {}; raw = []
    for i in range(len(d["ev_pdg"])):
        b = int(np.searchsorted(EDGES, th_ev[i]) - 1)
        if b < 0 or b >= len(EDGES) - 1:
            continue
        a = pr.event_angles(d, i, geo, N_GAS, ra.N_AERO)
        if a is not None:
            keep = rng.random(len(a["w"])) < a["w"]
            a = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in a.items()}
            a["w"] = np.ones(int(keep.sum()))
        if a is not None:  # per-track record for T15 steps 3-4: theta, phi, track position in the sensor frame, photon angles
            pv = np.array([d["ev_px"][i], d["ev_py"][i], d["ev_pz"][i]])
            E, _, _ = pr.helix_point(pv, 1 if d["ev_pdg"][i] > 0 else -1, {"r": pr.R_GAS_MID, "z": pr.Z_GAS_MID})
            hb = int(d["ev_hit_begin"][i]); cells = d["hit_cell"][hb:d["ev_hit_end"][i]]
            vals, cnt = np.unique(cells, return_counts=True); jj = hb + int(np.flatnonzero(cells == vals[np.argmax(cnt)])[0])
            gg = geo.cell(d["hit_system"][jj], d["hit_barrel"][jj], d["hit_cell"][jj])
            loc = ((E - gg["O"]) @ (gg["EX"] / np.linalg.norm(gg["EX"])), (E - gg["O"]) @ (gg["EY"] / np.linalg.norm(gg["EY"]))) if gg else (np.nan, np.nan)
            raw.append(dict(theta=float(th_ev[i]), phi=float(math.atan2(pv[1], pv[0])), loc=loc, gas=a["gas"].astype(np.float32),
                            aerogel=a["aerogel"].astype(np.float32)))
        dg = np.abs(a["gas"] - exp["gas"]) if (a and exp["gas"]) else None
        da = np.abs(a["aerogel"] - exp["aerogel"]) if (a and exp["aerogel"]) else None
        for rad in ("gas", "aerogel"):
            k = (b, rad); A = acc.setdefault(k, dict(track=[], npe=[], resid=[], th=[]))
            A["th"].append(th_ev[i])
            if a is None or exp[rad] is None:
                A["npe"].append(0.0); continue
            other = da if rad == "gas" else dg; mine = dg if rad == "gas" else da
            mask = mine < other if other is not None else np.ones(len(a["w"]), bool)
            r = pr.track_angle(a[rad][mask], a["w"][mask], exp[rad])
            if r is None:
                A["npe"].append(0.0); continue
            mean, npe, nph, th_used, w_used = r
            A["npe"].append(npe)
            if npe >= 3:
                A["track"].append(mean); A["resid"].extend((th_used - mean).tolist())
    return (sp, p), exp, acc, raw


def main():
    if sys.argv[1] == "worker":
        k, n, out = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        files = sorted(glob.glob(f"{NPZ}/*.npz"), key=lambda f: (-os.path.getsize(f), f))[k::n]
        parts = [run_file(f) for f in files]
        pickle.dump([x[:3] for x in parts], open(out, "wb"))
        pickle.dump([(x[0], x[1], x[3]) for x in parts], open(out.replace("part_", "raw_"), "wb"))
        print(f"worker {k}: {len(files)} files")
        return
    parts = [x for f in sys.argv[2:] for x in pickle.load(open(f, "rb"))]
    print("files", len(parts))
    merged = {}
    for key, exp, acc in parts:
        for (b, rad), A in acc.items():
            M = merged.setdefault((key, b, rad), dict(track=[], npe=[], resid=[], th=[], exp=exp[rad]))
            for c in ("track", "npe", "resid", "th"):
                M[c].extend(A[c])
    rows = []
    for ((sp, p), b, rad), M in merged.items():
        t, res = np.array(M["track"]), np.array(M["resid"]); n_ev = len(M["npe"]); ex = M["exp"]
        sig = ra.robust_sigma(t) if len(t) > 2 else float("nan")
        rows.append(dict(species=LABEL[sp], p_GeV=p, theta_lo=EDGES[b], theta_hi=EDGES[b + 1], theta_mean=float(np.mean(M["th"])), region=REGION[b],
                         radiator=rad, n_events=n_ev, n_tracks=len(t), eff_ring_found=len(t) / n_ev if n_ev else float("nan"),
                         mean_npe=float(np.mean(M["npe"])), rms_npe=float(np.std(M["npe"])),
                         mean_npe_found=float(np.mean([v for v in M["npe"] if v >= 3])) if len(t) else float("nan"),
                         mean_thetac_mrad=float(t.mean()) if len(t) else float("nan"), expected_thetac_mrad=ex if ex else float("nan"),
                         bias_mrad=float(t.mean() - ex) if (len(t) and ex) else float("nan"), sigma_thetac_mrad=sig,
                         err_mean_mrad=sig / math.sqrt(len(t)) if len(t) > 2 else float("nan"),
                         sigma1_mrad=ra.robust_sigma(res) if len(res) > 10 else float("nan")))
    rows.sort(key=lambda r: (r["species"], r["p_GeV"], r["theta_lo"], r["radiator"]))
    with open(f"{OUT}/arc_table_t06b.csv", "w") as fo:
        fo.write("# ALFA_o1_v01 ARC xenon 2 bar, full sim theta-uniform 12-90 deg (T06b), cell-sized theta bins; per-photon reco, PDE sampled;"
                 " n_gas %.6f n_aero %.5f\n" % (N_GAS, ra.N_AERO))
        w = csv.DictWriter(fo, fieldnames=list(rows[0])); w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    print("rows", len(rows))


if __name__ == "__main__":
    main()
