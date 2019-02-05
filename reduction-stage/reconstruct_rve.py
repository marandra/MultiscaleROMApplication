import argparse
import logging
import numpy
#import meshio

#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Reconstructs RVE displacement"
                                             " field from strain mode weights")
parser.add_argument('mesh_file', help="RVE mesh file (kratos mdpa format)")
parser.add_argument('LHS_matrix', help="modes_displacement file (BtB matrix")
parser.add_argument('RHS_matrix', help="modes_displacement file (BtPhi matrix")
parser.add_argument('weights_file', help="mode weights and macro strain file."
                         " Format: first column correspond macro strain,"
                         " following columsn are mode weights sorted by mode."
                         " Each row is a timestep.")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
args = parser.parse_args()

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S',level=verbosity_level)
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    logger.info("Reconstruct RVE fluctuant displacement field")

    logger.info("loading modes_displacement files")
    lhs = numpy.loadtxt(args.LHS_matrix)
    rhs = numpy.loadtxt(args.RHS_matrix)

    logger.info("loading modes_weights files")
    macro_strain = np.loadtxt(args.weights_file)[:, 0]
    weights = np.loadtxt(args.weights_file)[:, 1:]
    nr_timesteps = numpy.shape(weights)[0]
    nr_modes = numpy.shape(weights)[1] - 1  # fist column is macro strain
    logger.debug("number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("number of modes detected: {}".format(nr_modes))
    logger.debug("macro strain: ")
    logger.debug(macro_strain)
    logger.debug("mode weights: ")
    logger.debug(weights)

    for t in range(nr_timesteps):
        logger.info("solving fluctuant displacement for timestep {}".format(t))
        filename = "fluctuant_displacement-{:0>3d}".format(t)
        residual = numpy.dot(rhs, weights[t, :])
        uf = numpy.linalg.solve(lhs, rhs)
        numpy.savetxt(filename, uf)

