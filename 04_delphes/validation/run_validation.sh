#!/bin/bash
# [campaign_v2 copy (T17): v4 cards from campaign_v2/cards, outputs in campaign_v2/delphes_validation]
# T11: Delphes production with the two ALFA card variants. Guns: 10 chunks x 100k events per species (1M), Z->qq 2000.
# Each chunk gets its own Pythia seed and Delphes RandomSeed (card copy in cfg/).
# Usage: run_t11.sh [variants]   (default: Xe2bar C4F10; e.g. run_t11.sh C4F10 reruns only the C4F10 card and leaves the Xe files alone)
CARDS=${@:-Xe2bar C4F10}
D=/afs/cern.ch/work/h/hfatehi/alfa/delphes
W=/eos/user/h/hfatehi/alfa/campaign_v2/delphes_validation
CARDS_DIR=/eos/user/h/hfatehi/alfa/campaign_v2/cards
mkdir -p $W/cfg $W/runs $W/logs; cd $W
jobs=()
for card in $CARDS; do
  for sp in pip kp p; do
    for c in 0 1 2 3 4 5 6 7 8 9; do
      ci=1; [ $card = C4F10 ] && ci=2; si=1; [ $sp = kp ] && si=2; [ $sp = p ] && si=3
      seed=$(( 1000 * si + 10 * c + ci ))
      sed "s/^set RandomSeed .*/set RandomSeed $seed/" $CARDS_DIR/delphes_card_ALFA_$card.tcl > cfg/card_${card}_${sp}_$c.tcl
      sed -e 's/numberOfEvents = 20000/numberOfEvents = 100000/' -e 's/numberCount = 1000/numberCount = 20000/' $D/examples/ALFA/gun_$sp.cmnd > cfg/gun_${sp}_$c.cmnd
      printf 'Random:setSeed = on\nRandom:seed = %d\n' $seed >> cfg/gun_${sp}_$c.cmnd
      jobs+=("$card $sp $c")
    done
  done
  sed "s/^set RandomSeed .*/set RandomSeed 777/" $CARDS_DIR/delphes_card_ALFA_$card.tcl > cfg/card_${card}_zqq.tcl
done
run() { # card tag cmnd
  rm -f runs/$1_$2.root
  nice -n 10 /afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "cd $D; ./DelphesPythia8 $W/cfg/card_$1_$4.tcl $W/cfg/$3 $W/runs/$1_$2.root" > logs/$1_$2.log 2>&1
  echo "$1 $2 exit $?" >> logs/done.txt
}
: > logs/done.txt
for card in $CARDS; do run $card zqq zqq.cmnd zqq & done
for j in "${jobs[@]}"; do set -- $j; run $1 $2_$3 gun_$2_$3.cmnd $2_$3 & done
wait
echo ALL DONE >> logs/done.txt
