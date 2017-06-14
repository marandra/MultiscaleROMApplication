export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in trajectory_?
do
       	cd $i
	echo $i
	python3 MainKratos.py > outMainKratos 2>&1 &
	cd ..
done
