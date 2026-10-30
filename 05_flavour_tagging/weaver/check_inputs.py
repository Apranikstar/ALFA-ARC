"""T14 debug: iterate weaver's training dataset (all files, fetch-step 0.02, as in the trainings) and report
non-finite or extreme input values per batch. Usage: check_inputs.py <card> <yaml> <nbatches>"""
import sys, glob, numpy as np, torch
from weaver.utils.dataset import SimpleIterDataset
B = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa"
card, cfg, nb = sys.argv[1], sys.argv[2], int(sys.argv[3])
files = sorted(glob.glob(f"{B}/flat/{card}/train/*.root"))
ds = SimpleIterDataset({"": files}, cfg, for_training=True, load_range_and_fraction=((0, 1), 1, 1),
                       fetch_by_files=False, fetch_step=0.02, infinity_mode=False, in_memory=False, name="train")
dl = torch.utils.data.DataLoader(ds, batch_size=512, num_workers=4, drop_last=True, pin_memory=False)
worst = {}
for i, (X, y, _) in enumerate(dl):
    for k, v in X.items():
        v = v.numpy(); nf = (~np.isfinite(v)).sum(); mx = np.abs(v[np.isfinite(v)]).max() if nf < v.size else np.inf
        w = worst.setdefault(k, [0, 0.0]); w[0] += nf; w[1] = max(w[1], mx)
        if nf: print(f"batch {i} {k}: {nf} non-finite", flush=True)
    lab = y["_label_"].numpy()
    if i == 0: print("label counts", np.bincount(lab, minlength=6), "shapes", {k: tuple(v.shape) for k, v in X.items()}, flush=True)
    if i + 1 >= nb: break
print("batches", i + 1, {k: (int(a), float(b)) for k, (a, b) in worst.items()})
