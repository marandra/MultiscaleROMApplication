import numpy as np


def compute_permutations(n):
    u = np.empty((2**n, n))
    for d in range(2**n):
        bs = "{0:0{1}b}".format(d, n)
        u[d, :] = [float(x) for x in bs]
    return u

###################
if __name__ == "__main__":
    factor = 0.001
    n = 6
    U = compute_permutations(n)[1:,:]
    U *= factor
    for i, l in enumerate(U):
        string = 'm4 -DM4VAR_INITIALSTRAIN="{}" t_ProjectParameters.m4 > ProjectParameters_{:02d}.json'.format([x for x in l], i)
        print(string)

print('''
for i in {00..62}
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp materials.py $TRAJ
	cp t_model.mdpa $TRAJ/model.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
''')
