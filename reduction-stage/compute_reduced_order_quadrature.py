import sys
import configparser
import glob
import numpy as np
from compute_reduced_bases import ComputeROQ

def read_modes(file_name):
    c= np.loadtxt(file_name)
    #print(c.shape)
    return(c)

# get parameters
fname = sys.argv[1]
conf = configparser.ConfigParser()
conf.read(fname)
file_format = conf['Parameters']['file_format']
nr_elements = int(conf['Parameters']['nr_elements'])
nr_integration_points = int(conf['Parameters']['nr_integration_points'])
#nr_energy_reduced_modes = int(conf['Parameters']['nr_energy_reduced_modes'])
gaussweights_file_name = conf['Parameters']['gaussweights_file_name']

# get files
trajectory_paths = glob.glob('*_?')

gaussweights_file    = glob.glob(trajectory_paths[0] + '/' + gaussweights_file_name)
#print(gaussweights_file)

GaussWeights = np.fromfile(gaussweights_file[0], dtype=np.float32)

# Read the energy modes from the txt file
energy_basis_file    = './Output/EnergyBasis.dat'
energy_modes = read_modes(energy_basis_file)
#print(energy_modes.shape)

print("************************")
print("REDUCED ORDER QUADRATURE")
print("************************")

#num_ener_modes=100
nr_energy_reduced_modes = 100

#Ue_red=energy_modes[:,0:num_ener_modes]
if nr_energy_reduced_modes > energy_modes.shape[1]:
    sys.exit("Error: number of energy modes greater than the total number of computed energy modes")
else:
    num_ener_modes = nr_energy_reduced_modes

energy_modes_red=energy_modes[:,0:num_ener_modes]

factorLEQ=1.0
tol = 1e-10
nGP=energy_modes.shape[1] #In case of use the same number of points as energy modes.
#nGP=num_ener_modes #In case of use the same number of points as energy modes.

# TODO: Define a criteron to choose a subset of modes, or define a value (number of modes) as an input
[w,z] = ComputeROQ(energy_modes_red, GaussWeights, factorLEQ, nGP, tol)

# Print matrix with new weigths
roq_weigths=np.empty([nr_elements, nr_integration_points])
for i in range(nr_elements):
    for j in range(nr_integration_points):
        i_elem=nr_integration_points*(i)+j
        d=np.where(z==i_elem)[0]
        if not d:
            roq_weigths[i][j]=-1.0
        else:
            roq_weigths[i][j]=w[d[0]]

file_roq_weights='./Output/weights.dat'
with open(file_roq_weights,'wb') as ofile_roq_weights:
    np.savetxt(ofile_roq_weights, roq_weigths, fmt='%.13f')

#print("OK")