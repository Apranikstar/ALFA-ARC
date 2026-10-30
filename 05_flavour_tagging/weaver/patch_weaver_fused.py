"""Patch weaver-core 0.5.4 utils/data/tools.py: _get_content_and_offsets must pack the jagged array.
Without it, after concatenating/shuffling several files the raw content buffer does not start at offsets[0]=0 /
is longer than offsets[-1], and _fused_pad_and_stack either fails (broadcast error) or reads misaligned or
uninitialised entries in its numba kernel -> NaN loss. Usage: python3 patch_weaver_fused.py <site-packages dir>"""
import sys, os
p = os.path.join(sys.argv[1], "weaver/utils/data/tools.py")
s = open(p).read()
old = """        if isinstance(layout, (ak.contents.IndexedArray, ak.contents.IndexedOptionArray)):
            layout = ak.to_packed(a).layout
        # unwrap remaining"""
new = """        # ALFA T14 patch: always pack, so that offsets start at 0 and the content holds exactly offsets[-1] entries
        layout = ak.to_packed(a).layout
        # unwrap remaining"""
if new in s:
    print("already patched", p)
else:
    assert s.count(old) == 1, "pattern not found"
    open(p, "w").write(s.replace(old, new)); print("patched", p)
for d, _, fs in os.walk(os.path.join(sys.argv[1], "weaver")):
    for f in fs:
        if f.endswith(".nbi") or f.endswith(".nbc"):
            os.remove(os.path.join(d, f))
