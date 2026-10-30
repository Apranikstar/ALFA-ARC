#!/bin/bash
# T17 ZH production job (campaign_v2): one (card, decay, chunk) of 10k events, v4 card, same seeds as T14 (identical
# generated events, detector-only differences). Usage: prod_job.sh <Xe2bar|C4F10> <bb|cc|ss|uu|dd|gg> <chunk>
set -u
card=$1; d=$2; c=$3
T=/eos/user/h/hfatehi/alfa/campaign_v2/tagging
B2=/eos/experiment/fcc/ee/analyses/case-studies/higgs/electron_yukawa/alfa/campaign_v2
declare -A DEC=([bb]="25:onIfMatch = 5 -5" [cc]="25:onIfMatch = 4 -4" [ss]="25:onIfMatch = 3 -3" [uu]="25:onIfMatch = 2 -2" [dd]="25:onIfMatch = 1 -1" [gg]="25:onIfMatch = 21 21")
declare -A DID=([bb]=1 [cc]=2 [ss]=3 [uu]=4 [dd]=5 [gg]=6)
seed=$(( 100000 * ${DID[$d]} + c + 1 ))
HERE=$(pwd)
sed -e "s/@DECAY@/$d/; s/@NEV@/10000/; s/@HDECAY@/${DEC[$d]}/; s/@SEED@/$seed/" $T/cfg/zh_template.cmnd > zh.cmnd
sed "s/^set RandomSeed .*/set RandomSeed $seed/" /afs/cern.ch/work/h/hfatehi/alfa/delphes/cards/delphes_card_ALFA_$card.tcl > card.tcl
cp zh.cmnd $T/cfg/chunks/zh_H${d}_$c.cmnd 2>/dev/null; cp card.tcl $T/cfg/chunks/card_${card}_H${d}_$c.tcl 2>/dev/null
bash $T/delphes_edm4hep.sh $HERE/card.tcl $HERE/zh.cmnd $HERE/out.root > log.txt 2>&1
st=$?
tail -3 log.txt
[ $st -eq 0 ] && [ -s out.root ] && xrdcp -s -f out.root root://eosexperiment.cern.ch/$B2/edm4hep/$card/H$d/chunk_$c.root || exit 1
rm -f out.root
