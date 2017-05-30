for i in trajectory_?
do
       	cd $i
	echo $i
	python3 MainKratos.py
	cd ..
done
