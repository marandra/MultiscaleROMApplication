import argparse
import logging
import numpy
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities
import meshio


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(
    description="reconstructs fluctuant and total displacements fields"
)
parser.add_argument("mdpa_file", help="the .mdpa model used in training)")
parser.add_argument("correlation_strain", help="strain correlation matrix (.npy)")
parser.add_argument('correlation_energy', help="energy correlation matrix (.npy)")
parser.add_argument("rve_data", help="RVE reconstruction data file (.json)")
parser.add_argument(
    "-v", "--verbose", action="store_true", help="shows debug information"
)
args = parser.parse_args()

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=verbosity_level
)
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logger.info("Loading RVE node info")
    mesh = meshio.read(args.mdpa_file)
    rve_elems = {"hexahedron": mesh.cells["line8"]}
    #rve_elems = {"hexahedron": mesh.cells["line8"], "wedge": mesh.cells["line6"]}
    #rve_elems = [("hexahedron", mesh.cells["line8"]), ("wedge", mesh.cells["line6"]), ("hexahedron", mesh.cells["line8"])]
    rve_nodes = mesh.points

    logger.info("Loading strain-displacement correlation data")
    strain_correl = numpy.load(args.correlation_strain)

    logger.info("Loading energy correlation data")
    energy_correl = numpy.load(args.correlation_energy)

    logger.info("Loading RVE data")
    data = io_utilities.read_json(args.rve_data)
    rve_interpolation_params = numpy.array(data["interpolation_parameters"])
    logger.debug(rve_interpolation_params)
    logger.debug(numpy.shape(rve_interpolation_params))
    rve_macro_strain = numpy.array(data["macro_strain"])
    logger.debug(rve_macro_strain)
    logger.debug(numpy.shape(rve_macro_strain))
    rve_energy = numpy.array(data["strain_energy"])
    logger.debug(numpy.shape(rve_energy))

    nr_timesteps = numpy.shape(rve_interpolation_params)[0]
    nr_modes = numpy.shape(rve_interpolation_params)[1]
    logger.debug("Number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("Number of modes detected: {}".format(nr_modes))

    filename = "rve_reconstructed.xdmf"
    with meshio.XdmfTimeSeriesWriter(filename) as writer:
        writer.write_points_cells(rve_nodes, rve_elems)
        for t in range(nr_timesteps):
            logger.info("Timestep {}".format(t))

            logger.debug("Solving fluctuant displacement")
            displacement = numpy.dot(strain_correl[:, :nr_modes], rve_interpolation_params[t, :])
            displacement = numpy.reshape(displacement, (-1, 3))

            logger.debug("Solving total displacement")
            strain = rve_macro_strain[t, :]
            s_xx = strain[0]
            s_yy = strain[1]
            s_zz = strain[2]
            s_xy = 0.5 * strain[3]
            s_yz = 0.5 * strain[4]
            s_xz = 0.5 * strain[5]
            strain_tensor = numpy.array(
                [[s_xx, s_xy, s_xz], [s_xy, s_yy, s_yz], [s_yz, s_yz, s_zz]]
            )
            comp = numpy.dot(strain_tensor, rve_nodes.T)
            total_displacement = comp.T + displacement

            logger.info("Solving energy field")
            reduced_energy = rve_energy[t, :]
            energy = numpy.dot(energy_correl, reduced_energy)
            print(numpy.shape(energy_correl))
            # TEMP format to XDMF standard
            energy_in_elem = numpy.reshape(energy, (-1, 8))
            energy_list = [list(x) for x in list(energy_in_elem)]

            logger.debug("Writing timestep data")
            writer.write_data(
                t,
                point_data={
                    "FLUCTUANT_DISPLACEMENT": numpy.reshape(displacement, (-1, 3)),
                    "TOTAL_DISPLACEMENT": total_displacement,
                },
                cell_data={
                    #[("triangle", [[0, 1, 2], ...])]
                    #[("hexahedra", [[0, 1, 2, 3, 4, 5, 6, 7], ...]), ("wedge", [[0, 1, 2, 3, 4, 5],[],..])]
                    "STRAIN_ENERGY": {"hexahedron": energy_in_elem}
                    #numpy.reshape(displacement, (-1, 3)),
                },
            )

