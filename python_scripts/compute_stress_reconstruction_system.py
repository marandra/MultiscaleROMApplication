import argparse
import numpy
import logging
import postprocess_utilities as util

def load_rve_data(rve_data):
    logger.info("Reading reduced set integration points")
    #reduced_ip_set_filename = conf['Parameters']['reduced_ip_set_filename']
    #reduced_ip_set = util.read_numpy_file(reduced_ip_set_filename, 'ascii').astype(int)
    reduced_ip_set = rve_data["ip_global_id"]
    logger.info("Nr ip detected: {}".format(numpy.shape(reduced_ip_set)[0]))
    #logger.debug("{}".format(reduced_ip_set))

    logger.info("Reading reduced set integration weights")
    #rve_dataset_filename = conf['Parameters']['rve_dataset_filename']
    #reduced_ip_weights = numpy.array(util.read_json(rve_dataset_filename)['w'])
    reduced_ip_weights = numpy.array(rve_data["ip_weight"])
    logger.info("Nr weights detected: {}".format(numpy.shape(reduced_ip_weights)[0]))
    #logger.debug("{}".format(reduced_ip_weights))
    return reduced_ip_set, reduced_ip_weights


def load_energy_modes(energy_modes_filename, reduced_ip_set):
    logger.info("Reading energy modes")
    energy_modes_file_format = "binary"
    nr_ips = numpy.shape(reduced_ip_set)[0]
    energy_modes = util.read_numpy_file(energy_modes_filename, energy_modes_file_format)[:, :nr_ips - 1]
    reduced_energy_modes = energy_modes[reduced_ip_set, :]
    logger.info("{} IP detected, reduced to {}".format(numpy.shape(energy_modes)[0], numpy.shape(reduced_energy_modes)[0]))
    logger.info("Nr of modes loaded: {}".format(numpy.shape(reduced_energy_modes)[1]))
    #logger.debug(energy_modes)
    return energy_modes, reduced_energy_modes


def compute_stress_reconstruction_system(reduced_ip_set, reduced_ip_weights, energy_modes, reduced_energy_modes, integration_weights_filename ):
    logger.info("Computing COMPLETE system")
    logger.debug("-- A = reduced modes.T * weights * reduced modes")
    ip_weights = util.read_numpy_file(integration_weights_filename, 'ascii')
    ip_weights_diag = numpy.diag(ip_weights)
    reduced_ip_weights_diag = numpy.diag(reduced_ip_weights)
    weighted_energy_modes_transposed = numpy.dot(energy_modes.T, ip_weights_diag)
    weighted_reduced_energy_modes_transposed = numpy.dot(reduced_energy_modes.T, reduced_ip_weights_diag)
    A = numpy.dot(weighted_energy_modes_transposed, energy_modes)
    logger.debug("-- checking A is not singular")
    rankA = numpy.linalg.matrix_rank(A)
    logger.debug("A: {}".format(numpy.shape(A)))
    logger.debug("rank A: {}".format(numpy.linalg.matrix_rank(A)))
    if rankA != numpy.shape(A)[0]:
        logger.info("Matrix rank not complete (Too many ROQ points?). Aborting.")
        exit()
    logger.debug("-- inverse A")
    Ainv = numpy.linalg.inv(A)
    logger.debug("-- modes * invA * modes.T * weights ")
    aux_1 = numpy.dot(Ainv, weighted_reduced_energy_modes_transposed)
    return aux_1


def compute_system(rve_data_filename, energy_modes_filename, integration_weights_filename):
    rve_data = util.read_json(rve_data_filename)
    reduced_ip_set, reduced_ip_weights = load_rve_data(rve_data)
    energy_modes, reduced_energy_modes = load_energy_modes(energy_modes_filename, reduced_ip_set)
    A = compute_stress_reconstruction_system(reduced_ip_set, reduced_ip_weights, energy_modes, reduced_energy_modes, integration_weights_filename)
    return A


#######################################
# Main
#######################################

# configure logger
#verbosity_level = logging.INFO
verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',     datefmt='%H:%M:%S', level=verbosity_level)
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser(description="Computes     matrices for the reconstruction of the stress")
    parser.add_argument('rve_data', help="rve post-processing     data filename (.json)")
    parser.add_argument('integration_weights',     help="integration weights filename")
    parser.add_argument('energy_modes',     help="energy_modes_filename (binary .npy)")
    parser.add_argument('-v', '--verbose',     action="store_true", help="shows debug information")
    args = parser.parse_args()


    A = compute_system(args.rve_data, args.energy_modes, args.integration_weights)
    logger.info("Saving system")
    util.write_numpy_file('reconstruct_stress_binary.npy', 'binary', A)
