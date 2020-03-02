#!/bin/bash

for TRJ in trajectory_*
do
  echo "Launching job "$TRJ
  echo "export OMP_NUM_THREADS=1" > tmp_$TRJ.bash
  echo "export PYTHONPATH=${PYTHONPATH}" >> tmp_$TRJ.bash
  echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> tmp_$TRJ.bash
  echo "cd "$TRJ >> tmp_$TRJ.bash
  echo "/usr/bin/time -v -o time.dat python MainKratos.py > outMainKratos" >> tmp_$TRJ.bash
  echo "cd .." >> tmp_$TRJ.bash
  echo "rm tmp_$TRJ.bash" >> tmp_$TRJ.bash
done

for TRJ in trajectory_*
do
  bash tmp_$TRJ.bash
  sleep 0.05
done
