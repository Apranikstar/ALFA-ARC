#!/bin/bash
# T15: full ALFA (ARC xenon 2 bar), pi+ at fixed p, theta uniform over the central barrel cell (82.6-97.4 deg), with every
# particle kept (--part.keepAllParticles) so each ARC hit links to its optical photon: true emission point and direction.
# ~145 MB per event of EDM4hep (ECAL showers kept too): extract to .npz in the scratch dir, copy the .npz, drop the ROOT file.
# Usage: t15_job.sh <p_GeV> <n_events> <seed> <chunk>
set -u
P=$1; N=$2; SEED=$3; CH=$4
TAG=pip_p${P}_c${CH}
K=/afs/cern.ch/work/h/hfatehi/alfa/k4geo
X=/afs/cern.ch/work/h/hfatehi/alfa/condor/T15/t15_extract.py
OUT=root://eosuser.cern.ch//eos/user/h/hfatehi/alfa/work/T15/npz
HERE=$(pwd)
echo "host $(hostname) dir $HERE tag $TAG"
NIGHTLY_DAY=Thu /afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "cd $HERE; \
  ddsim --steeringFile $K/example/SteeringFile_ALFA_o1_v01_ARC.py --compactFile $K/FCCee/ALFA/compact/ALFA_o1_v01/ALFA_o1_v01.xml \
    --gun.particle pi+ --gun.energy '${P}*GeV' --gun.thetaMin '82.6*deg' --gun.thetaMax '97.4*deg' --gun.distribution uniform \
    --part.keepAllParticles True --random.seed $SEED -N $N --outputFile sim.edm4hep.root > sim_$TAG.log 2>&1; echo ddsim exit \$?; \
  python3 $X sim.edm4hep.root $TAG.npz"
status=$?
grep -E "Finished run" sim_$TAG.log
xrdcp -f $TAG.npz $OUT/$TAG.npz && xrdcp -f sim_$TAG.log $OUT/sim_$TAG.log
rm -f sim.edm4hep.root
[ -s $TAG.npz ] || exit 1
exit $status
