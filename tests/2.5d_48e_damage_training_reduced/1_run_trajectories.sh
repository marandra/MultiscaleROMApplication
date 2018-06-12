export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in 0 1 2 3 4 5 6
do
    TRJ="trajectory_"$i
   	cd $TRJ
	echo $TRJ
	python3 MainKratos.py > outMainKratos 2>&1 &
	cd ..
done
