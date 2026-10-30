"""Digitise the ALFA Reco-workshop (26 Oct 2026, slide 5) DELPHES vertex-resolution curves from a 400-dpi crop of the slide.
Axis calibration in crop pixels (read off the tick marks): d0 panel x(10 deg)=266, x(90 deg)=1236, y(1e2)=192, y(1e0)=696;
z0 panel x(10 deg)=1466, x(80 deg)=2316, y(1e3)=104, y(1e0)=702. Legend boxes are masked. Output: reco_ref.json {var: {p: [[theta, sigma_um], ...]}}"""
import json
import numpy as np
from PIL import Image
im = np.asarray(Image.open("reco_slide5_crop.png").convert("RGB")).astype(int)
COL = {1: (31, 119, 180), 5: (255, 127, 14), 46: (44, 160, 44), 182: (214, 39, 40)}  # matplotlib tab colours
PAN = {"d0": dict(x0=266, x1=1236, t0=10, t1=90, ya=192, la=2, yb=696, lb=0, xmax=1250, leg=(795, 300)),
       "z0": dict(x0=1466, x1=2316, t0=10, t1=80, ya=104, la=3, yb=702, lb=0, xmax=2399, leg=(2005, 300))}  # leg: legend box (x >, y <)
out = {}
for var, a in PAN.items():
    out[var] = {}
    for p, c in COL.items():
        pts = []
        for x in range(a["x0"], min(a["x1"], a["xmax"]) + 1, 10):
            col = im[:, x, :]
            d = np.abs(col - np.array(c)).sum(1)
            ys = np.where(d < 90)[0]
            ys = ys[(ys > min(a["ya"], a["yb"]) - 40) & (ys < a["yb"] + 5) & (ys > 60)]
            if x > a["leg"][0]: ys = ys[ys > a["leg"][1]]
            if not len(ys): continue
            y = np.median(ys)
            th = a["t0"] + (x - a["x0"]) / (a["x1"] - a["x0"]) * (a["t1"] - a["t0"])
            lg = a["la"] + (y - a["ya"]) / (a["yb"] - a["ya"]) * (a["lb"] - a["la"])
            pts.append([round(th, 2), round(10 ** lg, 3)])
        out[var][p] = pts
json.dump(out, open("reco_ref.json", "w"), indent=0)
for var in out:
    for p in out[var]:
        v = dict((round(t), s) for t, s in out[var][p])
        print(var, p, {t: v.get(t) for t in (10, 20, 30, 45, 60, 80, 90) if t in v})
