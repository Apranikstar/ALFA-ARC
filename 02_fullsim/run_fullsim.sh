#!/bin/bash
# campaign_v3 (aerogel n = 1.030, k4geo alfa-arc 2c138b3b): cell-averaged full sim, ALFA_o1_v01 with the ARC, fixed species /
# momentum, theta uniform in [12, 90] deg, uniform phi; extract ARC photons to .npz, copy to EOS, drop the ROOT file.
# Usage: fullsim_job.sh <Xe2bar|C4F10> <particle> <p_GeV> <n_events> <seed> <chunk>
set -u
GAS=$1; PART=$2; P=$3; N=$4; SEED=$5; CH=$6
TAG=${PART/+/p}_p${P}_c${CH}
K=/afs/cern.ch/work/h/hfatehi/alfa/k4geo
X=/afs/cern.ch/work/h/hfatehi/alfa/ALFA-ARC/02_fullsim/extract_arc_photons.py
case $GAS in Xe2bar) GEO=$K/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml;; C4F10) GEO=$K/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01_C4F10.xml;; esac   # C4F10: files copied there from 01_geometry/ALFA_o1_v01_C4F10
OUT=root://eosexperiment.cern.ch//eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v3/npz/$GAS
HERE=$(pwd)
echo "host $(hostname) dir $HERE gas $GAS tag $TAG"
NIGHTLY_DAY=Thu /afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "cd $HERE; \
  ddsim --steeringFile $K/example/SteeringFile_ALFA_o1_v01_ARC.py --compactFile $GEO \
    --gun.particle $PART --gun.energy '${P}*GeV' --gun.thetaMin '12*deg' --gun.thetaMax '90*deg' --gun.distribution uniform \
    --random.seed $SEED -N $N --outputFile sim.edm4hep.root > sim_$TAG.log 2>&1; echo ddsim exit \$?; \
  python3 $X sim.edm4hep.root $TAG.npz"
status=$?
grep -E "Finished run" sim_$TAG.log
xrdcp -f $TAG.npz $OUT/$TAG.npz && xrdcp -f sim_$TAG.log $OUT/sim_$TAG.log
rm -f sim.edm4hep.root
[ -s $TAG.npz ] || exit 1
exit $status
