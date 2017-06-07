import sys
import configparser
import glob
import numpy as np
import logging
from compute_reduced_bases import ComputeROQ

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
##handler = logging.FileHandler(sys.argv[1].rsplit('.', 1)[0] + '.log')
##handler.setLevel(logging.INFO)
##logger.addHandler(handler)


# get parameters
conf = configparser.ConfigParser()
conf.read(sys.argv[1])
nr_elements = int(conf['Parameters']['nr_elements'])
nr_integration_points = int(conf['Parameters']['nr_integration_points'])
nr_energy_reduced_modes = int(conf['Parameters']['nr_energy_reduced_modes'])
integration_weights_filename = conf['Parameters']['integration_weights_filename']
energy_bases_filename = conf['Parameters']['energy_bases_filename']
trajectory_filename = conf['Parameters']['trajectory_filename']
roq_weights_filename = conf['Parameters']['roq_weights_filename']

trajectory_paths = sorted(glob.glob("{}_?".format(trajectory_filename)))
integration_weights_files = glob.glob("{}_?/{}".format(trajectory_filename, integration_weights_filename))

logger.info("REDUCED ORDER QUADRATURE")

integration_weights = np.loadtxt(integration_weights_files[0])
energy_modes = np.loadtxt(energy_bases_filename)
if nr_energy_reduced_modes > energy_modes.shape[1]:
    sys.exit("Error: number of energy modes greater than the total number of computed energy modes")
energy_modes_red = energy_modes[:, 0:nr_energy_reduced_modes]

factorLEQ = 1.0
tol = 1e-10
nGP = energy_modes.shape[1] #In case of use the same number of points as energy modes.

# TODO: Define a criteron to choose a
#  subset of modes, or define a value (number of modes) as an input
[w,z] = ComputeROQ(energy_modes_red, integration_weights, factorLEQ, nGP, tol)

# Print matrix with new weigths
roq_weigths = np.empty([nr_elements, nr_integration_points])
for i in range(nr_elements):
    for j in range(nr_integration_points):
        i_elem = nr_integration_points * i + j
        d = np.where(z==i_elem)[0]
        if not d:
            roq_weigths[i][j] = -1.0
        else:
            roq_weigths[i][j] = w[d[0]]

with open(roq_weights_filename,'wb') as of:
    np.savetxt(of, roq_weigths, fmt='%.13f')

