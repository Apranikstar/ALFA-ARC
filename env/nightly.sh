#!/bin/bash
# Run a command in a clean shell with only the key4hep nightly + the local k4geo build.
# The login profile auto-loads the 2026-04-08 release, which blocks the nightly setup.
# Usage: ~/alfa/nightly.sh '<command string>'
# NIGHTLY_DAY=Thu ~/alfa/nightly.sh '...' pins the LCG nightly to that day's build (default: latest);
# batch jobs pin it so a midnight rollover cannot mix library versions with the local k4geo build.
exec env -i HOME="$HOME" USER="$USER" TERM=dumb KRB5CCNAME="$KRB5CCNAME" NIGHTLY_DAY="$NIGHTLY_DAY" bash --noprofile --norc -c '
if [ -n "$NIGHTLY_DAY" ]; then
  source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r "$NIGHTLY_DAY" >/dev/null 2>&1
else
  source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh >/dev/null 2>&1
fi
cd /afs/cern.ch/work/h/hfatehi/alfa/k4geo && k4_local_repo >/dev/null 2>&1
cd - >/dev/null 2>&1
'"$1"
