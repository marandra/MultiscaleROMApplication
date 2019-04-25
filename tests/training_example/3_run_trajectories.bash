#!/bin/bash
export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in {14..20}
do
    TRJ="trajectory_"$i
   	cd $TRJ
	echo $TRJ
	python3 MainKratos.py > outMainKratos 2>&1 &
	cd ..
done
