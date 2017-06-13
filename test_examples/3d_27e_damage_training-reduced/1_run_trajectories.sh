export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS="$OMP_NUM_THREADS
for i in trajectory_?
do
       	cd $i
	echo $i
	time python3 MainKratos.py >& outMainKratos &
	cd ..
done
