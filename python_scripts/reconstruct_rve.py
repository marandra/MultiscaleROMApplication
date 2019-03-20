import argparse
import logging
import numpy
#import meshio
import postprocess_utilities as util


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="reconstructs following RVE fields:"
   "fluctuant displacement, total displacement and stress")
parser.add_argument('gid_msh_file', help="gid output .msh file")
parser.add_argument('correlation', help="strain-displacemente correlation matric (.npy)")
parser.add_argument('rve_data', help="RVE reconstruction data file (.json)")
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

    logger.info("Loading RVE node info")
    rve_nodes = util.read_gid_msh_nodes(args.gid_msh_file)

    logger.info("Loading strain-displacement correlation data")
    strain_displ_correl = numpy.load(args.correlation)

    logger.info("Loading RVE data")
    data = util.read_json(args.rve_data)
    rve_interpolation_params = numpy.array(data["interpolation_parameters"])
    logger.debug(rve_interpolation_params)
    logger.debug(numpy.shape(rve_interpolation_params))
    rve_macro_strain = numpy.array(data["macro_strain"])
    logger.debug(rve_macro_strain)
    logger.debug(numpy.shape(rve_macro_strain))
    rve_stress = numpy.array(data["stress"])
    logger.debug(rve_stress)
    logger.debug(numpy.shape(rve_stress))

    nr_timesteps = numpy.shape(rve_interpolation_params)[0]
    nr_modes = numpy.shape(rve_interpolation_params)[1]
    logger.debug("Number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("Number of modes detected: {}".format(nr_modes))

    filename = args.gid_msh_file.rsplit(".", 1)[0] + ".res"
    f = util.write_gid_header(filename)
    for t in range(nr_timesteps):
        logger.info("Timestep {}".format(t))

        logger.info("Solving fluctuant displacement")
        displacement = numpy.dot(strain_displ_correl, rve_interpolation_params[t, :])
        displacement = numpy.reshape(displacement, (-1, 3))
        util.write_gid_vector_field(f, "FLUCTUANT_DISPLACEMENT",
                                    numpy.reshape(displacement, (-1, 3)), t)

        logger.info("Solving total displacement")
        strain = rve_macro_strain[t, :]
        s_xx = strain[0]
        s_yy = strain[1]
        s_zz = strain[2]
        s_xy = 0.5 * strain[3]
        s_yz = 0.5 * strain[4]
        s_xz = 0.5 * strain[5]
        strain_tensor = numpy.array([[s_xx, s_xy, s_xz],
                                   [s_xy, s_yy, s_yz],
                                   [s_yz, s_yz, s_zz]])
        comp = numpy.dot(strain_tensor, rve_nodes.T)
        total_displacement = comp.T + displacement
        util.write_gid_vector_field(f, "TOTAL_DISPLACEMENT", total_displacement, t)

        logger.info("Solving stress field")
        r_stress = rve_stress[t, :]

        util.write_gid_vector_field(f, "STRESS_TENSOR", total_displacement, t)

