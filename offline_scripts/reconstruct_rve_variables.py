import argparse
import logging
import numpy
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.MultiscaleROMApplication
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities
import meshio

def compute_elastic_tensor(E, NU):
    c1 = E / ((1 + NU) * (1 - 2 * NU))
    c2 = c1 * (1 - NU)
    c3 = c1 * NU
    c4 = c1 * 0.5 * (1 - 2 * NU)
    ConstitutiveMatrix = numpy.zeros((6, 6))
    ConstitutiveMatrix[0, 0] = c2
    ConstitutiveMatrix[0, 1] = c3
    ConstitutiveMatrix[0, 2] = c3
    ConstitutiveMatrix[1, 0] = c3
    ConstitutiveMatrix[1, 1] = c2
    ConstitutiveMatrix[1, 2] = c3
    ConstitutiveMatrix[2, 0] = c3
    ConstitutiveMatrix[2, 1] = c3
    ConstitutiveMatrix[2, 2] = c2
    ConstitutiveMatrix[3, 3] = c4
    ConstitutiveMatrix[4, 4] = c4
    ConstitutiveMatrix[5, 5] = c4
    return ConstitutiveMatrix


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(
    description="reconstructs fields"
)
parser.add_argument("mdpa_file", help="the .mdpa model used in training)")
parser.add_argument("correlation_strain", help="strain correlation matrix (.npy)")
parser.add_argument('correlation_energy', help="energy correlation matrix (.npy)")
parser.add_argument("runtime_data", help="multiscale runtime reconstruction data file (.json)")
parser.add_argument("rve_data", help="rve data file (.json)")
parser.add_argument("strain_modes", help="strain modes")
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

    model = KratosMultiphysics.Model()
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()


    logger.info("Loading RVE node info")
    mesh = meshio.read(args.mdpa_file)
    rve_elems = {"hexahedron": mesh.cells["line8"] + 1}
    #rve_elems = {"hexahedron": mesh.cells["line8"], "wedge": mesh.cells["line6"]}
    #rve_elems = [("hexahedron", mesh.cells["line8"]), ("wedge", mesh.cells["line6"]), ("hexahedron", mesh.cells["line8"])]
    rve_nodes = mesh.points

    logger.info("Loading data")
    strain_correl = numpy.load(args.correlation_strain)
    energy_correl = numpy.load(args.correlation_energy)
    data = io_utilities.read_json(args.runtime_data)

    logger.info("Loading rve data")
    rve_data = io_utilities.read_json(args.rve_data)
    props = rve_data["material_parameters"]["properties"]
    material_elastic_tensor = {}
    material_element_list = {}
    for m in props:
        material_name = m["model_part_name"]
        E = m["Material"]["Variables"]["YOUNG_MODULUS"]
        nu = m["Material"]["Variables"]["POISSON_RATIO"]
        material_elastic_tensor[material_name] = compute_elastic_tensor(E, nu)
        material_element_list[material_name] = []
        for elem in model[material_name].Elements:
            material_element_list[material_name].append(elem.Id)
    material_elem_map = {}
    for k, v in	material_element_list.items():
        for idx in v:
            material_elem_map[idx] = k

    modelpart = simulation._GetSolver().GetComputingModelPart()
    nr_comps = 6
    count = 0
    ip_elem_map = {}
    nr_of_ips = {}
    for elem in modelpart.Elements:
        ip_elem_map[elem.Id] = count
        nr_ip = len(elem.GetValuesOnIntegrationPoints(KratosMultiphysics.INTEGRATION_WEIGHT, modelpart.ProcessInfo))
        nr_of_ips[elem.Id] = nr_ip
        count += nr_ip * nr_comps

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
    strain_modes = numpy.load(args.strain_modes)[:, :nr_modes]

    filename = "rve_reconstructed.xdmf"
    with meshio.XdmfTimeSeriesWriter(filename) as writer:
        writer.write_points_cells(rve_nodes, rve_elems)
        for t in range(nr_timesteps):
            logger.info("Timestep {}".format(t))

            logger.debug("Solving fluctuant displacement")
            displacement = numpy.dot(strain_correl[:, :nr_modes], rve_interpolation_params[t, :])
            displacement = numpy.reshape(displacement, (-1, 3))

            logger.debug("Solving total displacement")
            strain_macro = rve_macro_strain[t, :]
            s_xx = strain_macro[0]
            s_yy = strain_macro[1]
            s_zz = strain_macro[2]
            s_xy = 0.5 * strain_macro[3]
            s_yz = 0.5 * strain_macro[4]
            s_xz = 0.5 * strain_macro[5]
            strain_macro_tensor = numpy.array(
                [[s_xx, s_xy, s_xz], [s_xy, s_yy, s_yz], [s_yz, s_yz, s_zz]]
            )
            comp = numpy.dot(strain_macro_tensor, rve_nodes.T)
            total_displacement = comp.T + displacement

            logger.debug("Solving energy field")
            reduced_energy = rve_energy[t, :]
            energy = numpy.dot(energy_correl, reduced_energy) # Vector. cada entrada es un punto de gauss
            # Hasta aqui Ok. Ahora chapuzas para visualizar en nodos
            energy_in_elem = numpy.reshape(energy, (-1, 8))
            mean_energy = numpy.mean(energy_in_elem, axis=1)
            mean_energy = numpy.reshape(mean_energy, (-1,1))

            logger.debug("Solving damage")
            damage_list = []
            elastic_energy_list = []
            strain_global = numpy.dot(strain_modes, rve_interpolation_params[t, :])
            for elem_id, nr_ips in nr_of_ips.items():
                C = material_elastic_tensor[material_elem_map[elem_id]]
                ip_0 = ip_elem_map[elem_id]
                #damage = 0
                mean_elastic_energy = 0
                for i in range(nr_ips):
                    strain = strain_global[ip_0 : ip_0 + nr_comps] + strain_macro
                    aux_1 = numpy.dot(C, strain)
                    elastic_energy = numpy.dot(strain.T, aux_1) / 2
                    #damage +=  (1 - energy_in_elem[elem_id - 1, i] / elastic_energy) / nr_ips
                    mean_elastic_energy +=  elastic_energy / nr_ips
                    ip_0 += nr_comps
                #damage_list.append(damage)
                damage = 1 - mean_energy / mean_elastic_energy
                damage_list.append(damage)
                elastic_energy_list.append(mean_elastic_energy)
            element_damage = numpy.array(damage_list).reshape((-1, 1)) # formatting for meshio
            element_elastic_energy = numpy.array(elastic_energy_list).reshape((-1, 1)) # formatting for meshio

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
                    "dummy_1": {"STRAIN_ENERGY": mean_energy},
                    "dummy_2": {"DAMAGE": element_damage},
                    "dummy_3": {"ELASTIC_ENERGY": element_elastic_energy},
                },
            )

