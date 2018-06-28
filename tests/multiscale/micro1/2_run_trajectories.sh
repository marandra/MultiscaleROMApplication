export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in 06 07 08 09 10 11
do
    TRJ="trajectory_"$i
   	cd $TRJ
	echo $TRJ
	python3 MainKratos.py > outMainKratos 2>&1 &
	cd ..
done
