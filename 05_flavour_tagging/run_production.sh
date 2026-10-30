#!/bin/bash
# [campaign_v2 copy (T17): v4 cards (Delphes alfa-arc 974fc61, T16 Kalman fix), output under /eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2]
# T14 production: e+e- -> ZH (Z -> nu nu, H -> bb/cc/ss/uu/dd/gg) at 240 GeV, 250k events per decay, Delphes ALFA card -> EDM4hep.
# Both cards (Xe2bar, C4F10) use the same Pythia seed per (decay, chunk): identical generated events, detector-only differences.
# Usage: run_production.sh [nparallel=60] [cards="Xe2bar C4F10"]
J=/eos/user/h/hfatehi/alfa/campaign_v2/tagging
B=/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2
CARDS=${2:-Xe2bar C4F10}; NPAR=${1:-60}; NCH=25; NEV=10000
declare -A DEC=([bb]="25:onIfMatch = 5 -5" [cc]="25:onIfMatch = 4 -4" [ss]="25:onIfMatch = 3 -3" [uu]="25:onIfMatch = 2 -2" [dd]="25:onIfMatch = 1 -1" [gg]="25:onIfMatch = 21 21")
declare -A DID=([bb]=1 [cc]=2 [ss]=3 [uu]=4 [dd]=5 [gg]=6)
mkdir -p $J/cfg/chunks $J/logs
: > $J/jobs.txt
for card in $CARDS; do
  for d in bb cc ss uu dd gg; do
    mkdir -p $B/edm4hep/$card/H$d
    for c in $(seq 0 $((NCH-1))); do
      seed=$(( 100000 * ${DID[$d]} + c + 1 ))
      cm=$J/cfg/chunks/zh_H${d}_$c.cmnd
      [ -f $cm ] || sed -e "s/@DECAY@/$d/; s/@NEV@/$NEV/; s/@HDECAY@/${DEC[$d]}/; s/@SEED@/$seed/" /eos/user/h/hfatehi/alfa/campaign_v2/tagging/cfg/zh_template.cmnd > $cm
      out=$B/edm4hep/$card/H$d/chunk_$c.root
      [ -s $out ] && continue
      echo "$card $d $c $cm $out" >> $J/jobs.txt
    done
  done
done
echo "jobs: $(wc -l < $J/jobs.txt)"
# one card file per chunk with its own Delphes RandomSeed (same for both gases)
run_one() {
  card=$1; d=$2; c=$3; cm=$4; out=$5
  J=/eos/user/h/hfatehi/alfa/campaign_v2/tagging
  seed=$(grep '^Random:seed' $cm | awk '{print $3}')
  cd=$J/cfg/chunks/card_${card}_H${d}_$c.tcl
  sed "s/^set RandomSeed .*/set RandomSeed $seed/" /afs/cern.ch/work/h/hfatehi/alfa/delphes/cards/delphes_card_ALFA_$card.tcl > $cd
  tmp=/tmp/hfatehi_t14_${card}_${d}_$c.root; rm -f $tmp
  nice -n 5 $J/delphes_edm4hep.sh $cd $cm $tmp > $J/logs/${card}_H${d}_$c.log 2>&1
  st=$?
  [ $st -eq 0 ] && [ -s $tmp ] && xrdcp -s -f $tmp root://eosexperiment.cern.ch/$out && rm -f $tmp
  echo "$card H$d $c exit $st" >> $J/logs/done.txt
}
export -f run_one
: > $J/logs/done.txt
cat $J/jobs.txt | xargs -P $NPAR -L 1 bash -c 'run_one "$@"' _
echo ALL DONE >> $J/logs/done.txt
