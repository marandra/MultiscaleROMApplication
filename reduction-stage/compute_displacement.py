import numpy as np
import sys
import pprint

# arguments:
# 1 strain bases
# 2 computed C (per timesteps)
# 3 mdpa

strain = np.load(sys.argv[1])
weight = np.loadtxt(sys.argv[2])

# generate epsilon_m #TODO read it
nr_steps = 250
end_time = 1.
f_eps_macro = [0.001, 0., 0., 0., 0., 0.]
time = np.linspace(end_time/nr_steps, end_time, nr_steps)
eps_macro = np.outer(time, f_eps_macro)

# read connectivity 
connect = []
flag_elements = False
with open(sys.argv[3], "r") as fi:
    for line in fi.readlines():
        if not flag_elements:
            if "Begin Elements" not in line:
                continue
            else:
                flag_elements = True
        else:
            if "End Elements" in line:
                break
            else:
                    connect.append([int(x) for x in line.split()[2:]])

nr_nodes = np.amax(connect)
nr_elements = len(connect)
nr_dofs = int(np.shape(strain)[0] / nr_elements)
print("nr_nodes: {}".format(nr_nodes))
print("nr_elements: {}".format(nr_elements))
print("nr_dofs: {}".format(nr_dofs))
K = np.zeros((nr_nodes, nr_nodes))
rhs = np.zeros(nr_nodes)
#
## temp loop
#for sm in eps_macro:
#
    for e, nodes in enumerate(connect):
        for i in range(nr_dofs) pass
            B = strain
            rhs[nodes[i]] = B[i]
         
#    sys.exit() 
