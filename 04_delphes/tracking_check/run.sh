#!/bin/bash
# run trkres.C in the nightly env against the local Delphes build: run.sh <geom.txt> <mode> <out> [N pmin pmax etamax seed]
D=/afs/cern.ch/work/h/hfatehi/alfa/delphes
W=/eos/user/h/hfatehi/alfa/work/T16
args="\"$1\",\"$2\",\"$3\""; shift 3; for a in "$@"; do args="$args,$a"; done
/afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "cd $W; root -l -b -q -e 'gSystem->Load(\"$D/libDelphes.so\"); gSystem->AddIncludePath(\"-I$D -I$D/external\");' 'trkres.C+($args)'"
