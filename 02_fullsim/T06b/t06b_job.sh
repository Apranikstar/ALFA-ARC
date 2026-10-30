#!/bin/bash
# T06b: cell-averaged full sim. ALFA_o1_v01 (ARC xenon 2 bar), fixed species / momentum, theta uniform in [12, 90] deg, uniform phi,
# so every barrel z-cell (pitch 25.7 cm) and endcap radius is sampled; the analysis bins in theta.
# Simulates in the job's scratch dir, extracts the ARC photons to .npz (T06 extractor), copies .npz + log to EOS, drops the ROOT file.
# Usage: t06b_job.sh <particle> <p_GeV> <n_events> <seed> <chunk>
set -u
PART=$1; P=$2; N=$3; SEED=$4; CH=$5
TAG=${PART/+/p}_p${P}_c${CH}
K=/afs/cern.ch/work/h/hfatehi/alfa/k4geo
X=/afs/cern.ch/work/h/hfatehi/alfa/condor/T06/t06_extract.py
OUT=root://eosuser.cern.ch//eos/user/h/hfatehi/alfa/work/T06b/npz
HERE=$(pwd)
echo "host $(hostname) dir $HERE tag $TAG"
NIGHTLY_DAY=Thu /afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "cd $HERE; echo stack \$KEY4HEP_STACK; \
  ddsim --steeringFile $K/example/SteeringFile_ALFA_o1_v01_ARC.py --compactFile $K/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml \
    --gun.particle $PART --gun.energy '${P}*GeV' --gun.thetaMin '12*deg' --gun.thetaMax '90*deg' --gun.distribution uniform \
    --random.seed $SEED -N $N --outputFile sim.edm4hep.root > sim_$TAG.log 2>&1; echo ddsim exit \$?; \
  python3 $X sim.edm4hep.root $TAG.npz"
status=$?
grep -E "Finished run" sim_$TAG.log
xrdcp -f $TAG.npz $OUT/$TAG.npz && xrdcp -f sim_$TAG.log $OUT/sim_$TAG.log
rm -f sim.edm4hep.root
[ -s $TAG.npz ] || exit 1
exit $status
