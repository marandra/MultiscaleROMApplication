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
parser.add_argument('mode_node_file', help="modes_displacement file (BtB matrix")
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

    logger.info("Loading modes_displacement files")
    mode_node_matrix = numpy.loadtxt(args.mode_node_file)

    logger.info("Loading modes_weights files")
    macro_strain = numpy.loadtxt(args.weights_file)[:, 0]
    weights = numpy.loadtxt(args.weights_file)[:, 1:]  # fist column is macro strain
    nr_timesteps = numpy.shape(weights)[0]
    nr_modes = numpy.shape(weights)[1]
    logger.debug("Number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("Number of modes detected: {}".format(nr_modes))
    logger.debug("Macro strain: ")
    logger.debug(macro_strain)
    logger.debug("Mode weights: ")
    logger.debug(weights)

    logger.info("Solving fluctuant displacement")
    for t in range(nr_timesteps):
        logger.info("Timestep {}".format(t))
        filename = "fluctuant_displacement-{:0>3d}".format(t)
        displacement = numpy.dot(mode_node_matrix, weights[t, :])
        displacement_form = numpy.reshape(displacement, (-1, 3))
        nnode = displacement_form.shape[0]
        print(nnode)
        gid_output = numpy.hstack([numpy.arange(1, nnode+1).reshape(-1,1).astype(int), displacement_form])
        numpy.savetxt(filename, gid_output)
