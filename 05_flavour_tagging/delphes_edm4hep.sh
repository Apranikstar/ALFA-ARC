#!/bin/bash
# DelphesPythia8_EDM4HEP of the local k4SimDelphes (built against the ARCPID Delphes), in the key4hep nightly.
# Usage: delphes_edm4hep.sh <card> <cmnd> <out.edm4hep.root>
D=/afs/cern.ch/work/h/hfatehi/alfa/delphes; K=/afs/cern.ch/work/h/hfatehi/alfa/k4SimDelphes/install
NIGHTLY_DAY=Thu /afs/cern.ch/work/h/hfatehi/alfa/nightly.sh "export LD_LIBRARY_PATH=$D:$K/lib64:\$LD_LIBRARY_PATH LD_PRELOAD=$K/lib64/libDelphesEDM4HepConverter.so:$D/libDelphes.so; cd $D; \
  $K/bin/DelphesPythia8_EDM4HEP $1 /eos/user/h/hfatehi/alfa/campaign_v2/tagging/cfg/edm4hep_output_config_ALFA.tcl $2 $3"
