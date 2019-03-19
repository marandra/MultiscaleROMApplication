import postprocess_utilities as util
import configparser
import argparse
import numpy
import logging
import postprocess_utilities as util

#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes matrices for the reconstruction of the stress")
parser.add_argument('rve_data', help="rve post-processing data filename (.json)")
parser.add_argument('integration_weights', help="integration weights filename")
parser.add_argument('energy_modes', help="energy_modes_filename (binary .npy)")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
args = parser.parse_args()

# parse configuration file
#conf = configparser.ConfigParser()
#conf.read(args.config_file)

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S', level=verbosity_level)
logger = logging.getLogger(__name__)

if __name__ == '__main__':

    rve_data = util.read_json(args.rve_data)

    logger.info("Reading reduced set integration points")
    #reduced_ip_set_filename = conf['Parameters']['reduced_ip_set_filename']
    #reduced_ip_set = util.read_numpy_file(reduced_ip_set_filename, 'ascii').astype(int)
    reduced_ip_set = rve_data["ip_global_id"]
    logger.info("Nr ip detected: {}".format(numpy.shape(reduced_ip_set)[0]))
    logger.debug("{}".format(reduced_ip_set))

    logger.info("Reading reduced set integration weights")
    #rve_dataset_filename = conf['Parameters']['rve_dataset_filename']
    #reduced_ip_weights = numpy.array(util.read_json(rve_dataset_filename)['w'])
    reduced_ip_weights = rve_data["ip_weight"]
    logger.info("Nr weights detected: {}".format(numpy.shape(reduced_ip_weights)[0]))
    logger.debug("{}".format(reduced_ip_weights))

    logger.info("Reading energy modes")
    #energy_modes_filename = conf['Parameters']['energy_bases_filename']
    energy_modes_filename = args.energy_modes
    energy_modes_file_format = "binary"
    nr_ips = numpy.shape(reduced_ip_set)[0]
    energy_modes = util.read_numpy_file(energy_modes_filename, energy_modes_file_format)[:, :nr_ips - 1]
    reduced_energy_modes = energy_modes[reduced_ip_set, :]
    logger.info("{} IP detected, reduced to {}".format(numpy.shape(energy_modes)[0], numpy.shape(reduced_energy_modes)[0]))
    logger.info("Nr of modes loaded: {}".format(numpy.shape(reduced_energy_modes)[1]))
    #logger.debug(energy_modes)

    # add a last mode corresponding to a volumetric reduced IP
    logger.info("Reading complete set integration weights")
    #integration_weights_filename = conf['Parameters']['integration_weights_filename']
    integration_weights_filename = args.integration_weights
    integration_weights = util.read_numpy_file(integration_weights_filename, 'ascii')
    volumetric_mode = numpy.sqrt(integration_weights.reshape(-1, 1))
    energy_modes = numpy.hstack([energy_modes, volumetric_mode])
    reduced_integration_weights = util.read_numpy_file(integration_weights_filename, 'ascii')[reduced_ip_set]
    reduced_volumetric_mode = numpy.sqrt(reduced_integration_weights.reshape(-1, 1))
    reduced_energy_modes = numpy.hstack([reduced_energy_modes, reduced_volumetric_mode])
    logger.info("Added volumetric mode")
    logger.info("Nr of modes: {}".format(numpy.shape(energy_modes)[1]))
    logger.info("Nr of modes (reduced): {}".format(numpy.shape(reduced_energy_modes)[1]))
    #logger.debug(reduced_energy_modes)

    logger.info("Computing system")
    logger.debug("-- A = reduced modes T * weights * reduced modes")
    weighted_reduced_energy_modes_transposed = numpy.multiply(reduced_energy_modes.T, reduced_ip_weights.reshape(-1, 1))
    A = numpy.dot(weighted_reduced_energy_modes_transposed, reduced_energy_modes)
    rankA = numpy.linalg.matrix_rank(A)
    logger.debug("A: {}".format(numpy.shape(A)))
    logger.debug("rank A: {}".format(numpy.linalg.matrix_rank(A)))
    if rankA != nr_ips:
        logger.info("Matrix rank not complete. Aborting.")
        exit()
    logger.debug("-- inverse A")
    Ainv = numpy.linalg.inv(A)
    logger.debug("-- modes * invA * reduced modes * weights ")
    A = numpy.dot(Ainv, weighted_reduced_energy_modes_transposed)
    A = numpy.dot(energy_modes, A)

    logger.info("Saving system")
    util.write_numpy_file('reconstruct_stress_ascii.npy', 'ascii', A)
    util.write_numpy_file('reconstruct_stress_binary.npy', 'binary', A)
