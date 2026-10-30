"""Extract the TrackCovariance DetectorGeometry block (and Bz) of a Delphes card into a text file for trkres.C.
Usage: python3 geom.py <card.tcl> <out.txt>"""
import re, sys
s = open(sys.argv[1]).read()
m = re.search(r"module TrackCovariance \w+ \{(.*?)\n\}", s, re.S)
blk = m.group(1)
g = re.search(r"set DetectorGeometry \{(.*?)\n\s*\}", blk, re.S).group(1)
bz = re.search(r"set Bz\s+(\S+)", blk).group(1)
if bz.startswith("$::"):
    bz = re.search(rf"^set {bz[3:]}\s+(\S+)", s, re.M).group(1)
lines = [l.split("#")[0].strip() for l in g.splitlines()]
lines = [l for l in lines if l]
open(sys.argv[2], "w").write(f"BZ {bz}\n" + "\n".join(lines) + "\n")
print(sys.argv[2], "Bz", bz, len(lines), "layers")
