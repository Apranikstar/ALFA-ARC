"""T14: flatten the FCCAnalyses event ntuples (2 Durham jets per event) to one entry per jet, with the flavour label of the sample.
Input:  <BASE>/ntuples/<card>/H<f>/chunk<i>.root   (tree "events"; jet_* per jet, pfcand_* per jet and constituent)
Output: <BASE>/flat/<card>/{train,test}/H<f>_<i>.root (tree "tree"; one jet per entry: jet scalars, constituent arrays, labels)
        chunks 0-19 -> train, 20-24 -> test (80/20, split by event so the two jets of an event stay together)
Labels as in the FCCAnalyses weaver stage2: recojet_isB/C/S/U/D/G (one-hot, flavour from "H<f><f>" in the path), plus
jet_npfcand (constituents per jet, as in stage2) and jet_flavour 0..5. Same content as stage2.py, written with uproot (columnar).
Usage: python3 flatten_jets.py <card> [nproc]"""
import glob, os, re, sys
from concurrent.futures import ProcessPoolExecutor

BASE = "/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2"  # campaign_v2 (T17)
FLAV = ["b", "c", "s", "u", "d", "g"]
N_TRAIN_CHUNKS = 20


def flatten(args):
    import awkward as ak
    import numpy as np
    import uproot
    card, path = args
    f = re.search(r"/H(\w\w)/chunk_?(\d+)\.root$", path)
    pair, chunk = f.group(1), int(f.group(2))
    flav = pair[0]  # "bb" -> "b", as stage2 matches "H{f}{f}"
    split = "train" if chunk < N_TRAIN_CHUNKS else "test"
    out = f"{BASE}/flat/{card}/{split}/H{pair}_{chunk}.root"
    t = uproot.open(path)["events"]
    keys = [k for k in t.keys() if k.startswith(("jet_", "pfcand_", "pfmask"))]
    a = t.arrays(keys)
    njet = ak.num(a["jet_e"])
    data = {}
    for k in keys:
        x = a[k]
        data[k] = ak.flatten(x, axis=1)  # [event][jet](...) -> [jet](...)
    n = len(data["jet_e"])
    ev = np.repeat(np.arange(len(njet)) + 100000 * chunk, ak.to_numpy(njet))
    for i, fl in enumerate(FLAV):
        data[f"recojet_is{fl.upper()}"] = np.full(n, 1 if fl == flav else 0, dtype=np.int32)
    data["jet_flavour"] = np.full(n, FLAV.index(flav), dtype=np.int32)
    data["jet_npfcand"] = ak.to_numpy(ak.num(data["pfcand_e"], axis=1)).astype(np.int32)
    data["jet_index"] = ak.to_numpy(ak.flatten(ak.local_index(a["jet_e"], axis=1)))
    data["event_index"] = ev.astype(np.int64)
    # constituent arrays as float32 (weaver reads them as jagged arrays); scalars as float32/int32
    for k in list(data):
        x = data[k]
        if k.startswith(("pfcand_", "pfmask")):
            data[k] = ak.values_astype(x, np.float32)
        elif k.startswith("jet_") and k not in ("jet_flavour", "jet_index", "jet_npfcand"):
            data[k] = ak.to_numpy(x).astype(np.float32)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = f"/tmp/hfatehi_flat_{card}_{pair}_{chunk}.root"
    with uproot.recreate(tmp, compression=uproot.LZ4(4)) as fo:
        fo["tree"] = {k: v for k, v in data.items()}
    os.system(f"xrdcp -s -f {tmp} root://eospublic.cern.ch/{out} && rm -f {tmp}")
    return out, n


if __name__ == "__main__":
    card = sys.argv[1]; nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    files = sorted(glob.glob(f"{BASE}/ntuples/{card}/H*/chunk*.root"))
    print(card, len(files), "input files")
    tot = 0
    with ProcessPoolExecutor(nproc) as ex:
        for out, n in ex.map(flatten, [(card, f) for f in files]):
            tot += n
    print(card, "jets written", tot)
