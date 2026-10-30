#!/bin/bash
# T14: one ParT training on a condor GPU node. Usage: train_condor.sh <name> <card: Xe2bar|C4F10> <config: baseline|arc> <epochs>
set -u
NAME=$1; CARD=$2; CFG=$3; EPOCHS=$4
B=/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2
W=/eos/user/h/hfatehi/alfa/campaign_v2/tagging/weaver
P=/afs/cern.ch/work/h/hfatehi/alfa/weaver_py
OUT=$B/training/$NAME; mkdir -p $OUT
echo "host $(hostname) gpu $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null) start $(date)"
# the LCG setup script uses unset variables: no set -u while sourcing it (it silently killed the first submission)
set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_110_cuda/x86_64-el9-gcc13-opt/setup.sh >/dev/null; set -u
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 PYTHONUSERBASE=$P PYTHONPATH=$P/lib/python3.13/site-packages:${PYTHONPATH:-}
cp $W/jets_$CFG.yaml $W/ParT.py $OUT/
cd $OUT
${WEAVER_DRYRUN:+echo} $P/bin/weaver \
  --data-train "$B/flat/$CARD/train/*.root" \
  --data-test "$B/flat/$CARD/test/*.root" \
  --data-config $OUT/jets_$CFG.yaml --network-config $OUT/ParT.py \
  --model-prefix $OUT/net --num-workers 6 --fetch-step 0.02 \
  --batch-size 512 --start-lr 1e-3 --num-epochs $EPOCHS --optimizer ranger --use-amp --gpus 0 \
  --predict-output $OUT/pred.root --log $OUT/train.log
st=$?
echo "weaver exit $st end $(date)"
exit $st
