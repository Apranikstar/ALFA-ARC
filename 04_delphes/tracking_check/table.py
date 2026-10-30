"""T16 before/after table at 90 deg: sigma(d0) [um], sigma(z0) [um], sigma(pT)/pT for 1, 5, 10, 45 GeV, with IDEA and the
Reco-workshop values (digitised, 46 GeV used for 45; no 10 GeV curve there). Also the material up to the ARC (r < 1.9 m)."""
import json
def g(fn):
    d = {}
    for l in open(fn):
        if l.startswith(("#", "S")): continue
        v = [float(t) for t in l.split()]; d[(v[0], v[1])] = v
    return d
B, N, I = g("grid_alfa_v0.txt"), g("grid_alfa_new.txt"), g("grid_idea.txt")
ref = json.load(open("reco_ref.json"))
def rv(var, p):
    pts = ref[var].get(str(p))
    return min(pts, key=lambda x: abs(x[0] - 90))[1] if pts else None
print("p[GeV] | sd0 before after IDEA Reco | sz0 before after IDEA Reco | dpT/pT before after IDEA")
for p in (1, 5, 10, 45):
    b, n, i = B[(p, 90)], N[(p, 90)], I[(p, 90)]
    pr = 46 if p == 45 else p
    r0, rz = rv("d0", pr), rv("z0", pr)
    f = lambda x: f"{x:.2f}" if x else "-"
    print(f"{p:>4} | {b[4]:.2f} {n[4]:.2f} {i[4]:.2f} {f(r0)} | {b[5]:.2f} {n[5]:.2f} {i[5]:.2f} {f(rz)} | {b[6]:.5f} {n[6]:.5f} {i[6]:.5f}")
for name, fn in [("before", "alfa_xe.txt"), ("after", "alfa_xe_new.txt")]:
    tot = {}
    for l in open(fn).read().splitlines()[1:]:
        x = l.split()
        if x[0] == "1" and float(x[4]) < 1.9: tot[x[1]] = tot.get(x[1], 0) + float(x[5]) / float(x[6]) * 100
    print(name, "x/X0 at 90 deg up to r = 1.9 m [%]:", {k: round(v, 3) for k, v in tot.items()}, "total", round(sum(tot.values()), 2))
