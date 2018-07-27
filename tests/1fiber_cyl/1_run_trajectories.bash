#!/bin/bash
export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in {07..20}
do
    TRJ="trajectory_"$i
   	cd $TRJ
	echo $TRJ
	time python3 MainKratos.py > outMainKratos 2>&1 &
	cd ..
done
