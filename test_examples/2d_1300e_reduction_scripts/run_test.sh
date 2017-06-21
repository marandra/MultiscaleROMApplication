export OMP_NUM_THREADS=1
cp -r ~/Cases/validation_files/trajectory_? .
python3 ../../reduction-stage/compute_reduced_bases.py reduced_bases.cfg

cp ~/Cases/validation_files/energy_bases.dat.matlab .
python3 ../../reduction-stage/compute_roq.py reduced_bases_roq.cfg
