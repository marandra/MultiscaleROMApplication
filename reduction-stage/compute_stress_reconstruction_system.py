import postprocess_utilities as util
import configparser
import argparse
import numpy
import logging

#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes matrices for the reconstruction of the stress")
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
args = parser.parse_args()

# parse configuration file
conf = configparser.ConfigParser()
conf.read(args.config_file)

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s', datefmt='%H:%M:%S', level=verbosity_level)
logger = logging.getLogger(__name__)

if __name__ == '__main__':

    logger.info("Reading reduced set integration points")
    reduced_ip_set_filename = conf['Parameters']['reduced_ip_set_filename']
    reduced_ip_set = util.read_numpy_file(reduced_ip_set_filename, 'ascii').astype(int)
    logger.info("Nr ip detected: {}".format(numpy.shape(reduced_ip_set)[0]))
    logger.debug("{}".format(reduced_ip_set))

    logger.info("Reading reduced set integration weights")
    rve_dataset_filename = conf['Parameters']['rve_dataset_filename']
    reduced_ip_weights = numpy.array(util.read_json(rve_dataset_filename)['w'])
    logger.info("Nr weights detected: {}".format(numpy.shape(reduced_ip_weights)[0]))
    logger.debug("{}".format(reduced_ip_weights))

    logger.info("Reading energy modes")
    energy_modes_filename = conf['Parameters']['energy_bases_filename']
    energy_modes_file_format = conf['Parameters']['bases_file_format']
    nr_ips = numpy.shape(reduced_ip_set)[0]
    energy_modes = util.read_numpy_file(energy_modes_filename, energy_modes_file_format)[:, :nr_ips - 1]
    reduced_energy_modes = energy_modes[reduced_ip_set, :]
    logger.info("{} IP detected, reduced to {}".format(numpy.shape(energy_modes)[0], numpy.shape(reduced_energy_modes)[0]))
    logger.info("Nr of modes loaded: {}".format(numpy.shape(reduced_energy_modes)[1]))
    #logger.debug(energy_modes)

    # add a last mode corresponding to a volumetric reduced IP
    logger.info("Reading complete set integration weights")
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    integration_weights = util.read_numpy_file(integration_weights_filename, 'ascii')[reduced_ip_set]
    volumetric_mode = numpy.sqrt(integration_weights.reshape(-1, 1))
    reduced_energy_modes = numpy.hstack([reduced_energy_modes, volumetric_mode])
    logger.info("Added volumetric mode")
    logger.info("Nr of modes: {}".format(numpy.shape(reduced_energy_modes)[1]))
    #logger.debug(reduced_energy_modes)

    logger.info("Computing system")
    logger.debug("- A = reduced modes T * weights * reduced modes")
    weighted_reduced_energy_modes_transposed = numpy.multiply(reduced_energy_modes.T, reduced_ip_weights.reshape(-1, 1))
    A = numpy.dot(weighted_reduced_energy_modes_transposed, reduced_energy_modes)
    logger.debug("A: {}".format(numpy.shape(A)))
    logger.debug("- inverse A")
    Ainv = numpy.linalg.inv(A)
    logger.debug("- modes * invA * reduced modes * weights ")
    A = numpy.dot(Ainv, weighted_reduced_energy_modes_transposed)
    A = numpy.dot(energy_modes, A)

    logger.info("Saving system")
    util.write_numpy_file('reconstruct_stress.npy', 'ascii', A)
