import argparse
import numpy
import logging
import json


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def load_rve_data(rve_data):
    logger.info("Reading reduced set integration points")
    # reduced_ip_set_filename = conf['Parameters']['reduced_ip_set_filename']
    # reduced_ip_set = util.read_numpy_file(reduced_ip_set_filename, 'ascii').astype(int)
    reduced_ip_set = rve_data["ip_global_id"]
    logger.info("Nr ip detected: {}".format(numpy.shape(reduced_ip_set)[0]))
    # logger.debug("{}".format(reduced_ip_set))

    logger.info("Reading reduced set integration weights")
    # rve_dataset_filename = conf['Parameters']['rve_dataset_filename']
    # reduced_ip_weights = numpy.array(util.read_json(rve_dataset_filename)['w'])
    reduced_ip_weights = numpy.array(rve_data["ip_weight"])
    logger.info("Nr weights detected: {}".format(numpy.shape(reduced_ip_weights)[0]))
    # logger.debug("{}".format(reduced_ip_weights))
    return reduced_ip_set, reduced_ip_weights


def load_energy_modes(modes_filename, reduced_ip_set, nr_modes):
    # modes = io_utilities.read_numpy_file(modes_filename, "binary")[:, :nr_modes]
    modes = numpy.load(modes_filename)[:, :nr_modes]
    reduced_modes = modes[reduced_ip_set, :]
    logger.info(
        "Modes matrix {} {} - Reduced modes matrix: {} {}".format(
            numpy.shape(modes)[0],
            numpy.shape(modes)[1],
            numpy.shape(reduced_modes)[0],
            numpy.shape(reduced_modes)[1],
        )
    )
    return modes, reduced_modes


def compute_reconstruction_system(
    reduced_ip_set, reduced_ip_weights, energy_modes, reduced_energy_modes
):
    logger.info("Computing COMPLETE system")
    logger.debug("-- A = reduced modes.T * weights * reduced modes")
    reduced_ip_weights_diag = numpy.diag(reduced_ip_weights)

    weighted_reduced_modes_transposed = numpy.dot(
        reduced_energy_modes.T, reduced_ip_weights_diag
    )
    A = numpy.dot(weighted_reduced_modes_transposed, reduced_energy_modes)

    logger.debug("-- checking A is not singular")
    rankA = numpy.linalg.matrix_rank(A)
    logger.debug("A: {}".format(numpy.shape(A)))
    logger.debug("rank A: {}".format(numpy.linalg.matrix_rank(A)))
    if rankA != numpy.shape(A)[0]:
        logger.info("Matrix rank not complete (Too many ROC points?). Aborting.")
        exit()
    logger.debug("-- inverse A")
    Ainv = numpy.linalg.inv(A)

    logger.debug("-- modes * invA * modes.T * weights ")
    aux_1 = numpy.dot(Ainv, weighted_reduced_modes_transposed)
    aux_2 = numpy.dot(energy_modes, aux_1)
    return aux_2


def compute_system(rve_data_filename, energy_modes_filename, nr_modes):
    rve_data = read_json(rve_data_filename)
    reduced_ip_set, reduced_ip_weights = load_rve_data(rve_data)
    modes, reduced_modes = load_energy_modes(
        energy_modes_filename, reduced_ip_set, nr_modes
    )
    A = compute_reconstruction_system(
        reduced_ip_set, reduced_ip_weights, modes, reduced_modes
    )
    return A


#######################################
# Main
#######################################

# configure logger
# verbosity_level = logging.INFO
verbosity_level = logging.DEBUG
logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=verbosity_level
)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="Computes data necessary for later reconstruction of the damage"
    )
    parser.add_argument("rve_data", help="rve data filename")
    # parser.add_argument('damage_modes',     help="damage_modes_filename (binary .npy)")
    parser.add_argument("r_value_bases", help="r_value bases filename")
    parser.add_argument("nr_modes", help="nr of modes")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="shows debug information"
    )
    args = parser.parse_args()

    # A = compute_system(args.rve_data, args.energy_modes, args.integration_weights)
    A = compute_system(args.rve_data, args.r_value_bases, int(args.nr_modes))
    logger.info("Saving system")
    # io_utilities.write_numpy_file(
    #    "correlation_r_value_{}.npy".format(args.nr_modes), "binary", A
    # )
    numpy.save("correlation_r_value_{}.npy".format(args.nr_modes), A)
