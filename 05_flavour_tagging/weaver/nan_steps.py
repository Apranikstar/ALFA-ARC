"""T14 debug: ParT training steps with the real data, report the first non-finite loss/output/grad.
Usage: nan_steps.py <yaml> <amp 0|1> <lr> <nsteps>"""
import sys, glob, math, importlib.util, numpy as np, torch
from weaver.utils.dataset import SimpleIterDataset
B = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa"
cfg, amp, lr, ns = sys.argv[1], sys.argv[2] == "1", float(sys.argv[3]), int(sys.argv[4])
files = sorted(glob.glob(f"{B}/flat/Xe2bar/train/*.root"))
ds = SimpleIterDataset({"": files}, cfg, for_training=True, load_range_and_fraction=((0, 1), 1, 1),
                       fetch_by_files=False, fetch_step=0.02, infinity_mode=False, in_memory=False, name="train")
dl = torch.utils.data.DataLoader(ds, batch_size=512, num_workers=3, drop_last=True)
spec = importlib.util.spec_from_file_location("net", "/eos/user/h/hfatehi/alfa/work/T14/weaver/ParT.py"); net = importlib.util.module_from_spec(spec); spec.loader.exec_module(net)
model, _ = net.get_model(ds.config); model = model.cuda(); model.train()
from weaver.utils.nn.optimizer.ranger import Ranger
opt = Ranger(model.parameters(), lr=lr)
scaler = torch.amp.GradScaler("cuda", enabled=amp)
lossf = torch.nn.CrossEntropyLoss()
for i, (X, y, _) in enumerate(dl):
    ins = [X[k].cuda() for k in ds.config.input_names]; lab = y["_label_"].long().cuda()
    opt.zero_grad()
    with torch.autocast("cuda", enabled=amp):
        out = model(*ins); loss = lossf(out, lab)
    bad_out = (~torch.isfinite(out)).sum().item()
    scaler.scale(loss).backward(); scaler.unscale_(opt)
    gn = torch.nn.utils.clip_grad_norm_(model.parameters(), 1e9).item()
    bad_p = [n for n, p in model.named_parameters() if not torch.isfinite(p).all()]
    if i < 5 or i % 25 == 0 or not math.isfinite(loss.item()) or bad_out or bad_p:
        print(f"step {i} loss {loss.item():.4f} nonfinite_out {bad_out} gradnorm {gn:.3g} scale {scaler.get_scale() if amp else 1} badparams {bad_p[:3]}", flush=True)
    if bad_p or not math.isfinite(loss.item()):
        # where does it start? forward hooks on the bad batch
        print("first non-finite at step", i)
        pairwise_lv_fts = model.mod.pair_embed.pairwise_lv_fts
        v = X["pf_vectors"].double(); m = X["pf_mask"]; f = X["pf_features"]
        pf = pairwise_lv_fts(v.unsqueeze(-1), v.unsqueeze(-2))  # (N, 4, P, P)
        badj = (~torch.isfinite(pf)).flatten(1).any(1).nonzero().flatten().tolist()
        print("jets with non-finite pair features:", badj[:10])
        for j in badj[:3]:
            n = int(m[j].sum()); print(" jet", j, "nconst", n, "label", int(lab[j]))
            print("  vectors px,py,pz,e:", v[j, :, :max(n, 1)].T.tolist()[:4])
            print("  which pair fts bad:", (~torch.isfinite(pf[j])).flatten(1).any(1).tolist())
        print("features non-finite:", (~torch.isfinite(f)).sum().item(), "empty jets:", (m.flatten(1).sum(1) == 0).sum().item(), "1-const jets:", (m.flatten(1).sum(1) == 1).sum().item())
        break
    scaler.step(opt); scaler.update()
    if i + 1 >= ns: break
