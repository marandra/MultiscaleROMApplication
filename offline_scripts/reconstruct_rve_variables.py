import argparse
import logging
import numpy
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)
import KratosMultiphysics.MultiscaleROMApplication
import meshio
import math
import json


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def q(r, E, yield_stress, inf_yield_stress, H0, H1):
    r0 = yield_stress / math.sqrt(E)
    q0 = r0  # strain_variable_init
    q1 = inf_yield_stress / math.sqrt(E)  # stress_variable_inf
    r1 = r0 + (q1 - q0) / H0
    if r < r0:
        return q0
    if r >= r0 and r < r1:
        return q0 + H0 * (r - r0)
    # Case r >= r1:
    return q1 + H1 * (r - r1)


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
parser = argparse.ArgumentParser(description="reconstructs fields")
parser.add_argument("mdpa_file", help="the .mdpa model used in training)")
parser.add_argument("correlation_strain", help="strain correlation matrix (.npy)")
parser.add_argument("correlation_r_value", help="r_value correlation matrix (.npy)")
parser.add_argument(
    "runtime_data", help="multiscale runtime reconstruction data file (.json)"
)
parser.add_argument("rve_data", help="rve data file (.json)")
parser.add_argument("strain_modes", help="strain modes")
parser.add_argument("-v", "--verbose", action="store_true", help="sdebug information")
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
    with open("ProjectParameters.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()

    logger.info("Loading RVE node info")
    mesh = meshio.read(args.mdpa_file)
    rve_cells = []
    for cell_block in mesh.cells:
        element_type = cell_block[0]
        #if "hexa" in element_type or "wedge" in element_type:
        if "line8" in element_type:
            rve_cells.append(meshio.CellBlock("hexahedron", cell_block[1]))
        if "line6" in element_type:
            rve_cells.append(meshio.CellBlock("wedge", cell_block[1]))

    logger.info("Loading data")
    strain_correl = numpy.load(args.correlation_strain)
    r_value_correl = numpy.load(args.correlation_r_value)
    data = read_json(args.runtime_data)

    logger.info("Loading rve data")
    rve_data = read_json(args.rve_data)
    props = rve_data["material_parameters"]["properties"]
    material_properties = {}
    material_element_list = {}
    for m in props:
        material_name = m["model_part_name"]
        material_properties[material_name] = {}
        E = m["Material"]["Variables"]["YOUNG_MODULUS"]
        nu = m["Material"]["Variables"]["POISSON_RATIO"]
        yield_stress = m["Material"]["Variables"]["YIELD_STRESS"]
        inf_yield_stress = m["Material"]["Variables"]["INFINITY_YIELD_STRESS"]
        H0 = m["Material"]["Variables"]["HARDENING_MODULI_VECTOR"][0]
        H1 = m["Material"]["Variables"]["HARDENING_MODULI_VECTOR"][1]
        material_properties[material_name]["E"] = E
        material_properties[material_name]["nu"] = nu
        material_properties[material_name]["yield_stress"] = yield_stress
        material_properties[material_name]["inf_yield_stress"] = inf_yield_stress
        material_properties[material_name]["H0"] = H0
        material_properties[material_name]["H1"] = H1
        material_properties[material_name]["C"] = compute_elastic_tensor(E, nu)

        material_element_list[material_name] = []
        for elem in model[material_name].Elements:
            material_element_list[material_name].append(elem.Id)
    material_elem_map = {}
    for k, v in material_element_list.items():
        for idx in v:
            material_elem_map[idx] = k

    modelpart = simulation._GetSolver().GetComputingModelPart()
    nr_comps = 6
    count = 0
    ip_elem_map = {}
    nr_of_ips = {}
    for elem in modelpart.Elements:
        ip_elem_map[elem.Id] = count
        nr_ip = len(
            elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.INTEGRATION_WEIGHT, modelpart.ProcessInfo
            )
        )
        nr_of_ips[elem.Id] = nr_ip
        count += nr_ip * nr_comps

    rve_interpolation_params = numpy.array(data["interpolation_parameters"])
    logger.debug(rve_interpolation_params)
    logger.debug(numpy.shape(rve_interpolation_params))
    rve_macro_strain = numpy.array(data["macro_strain"])
    logger.debug(rve_macro_strain)
    logger.debug(numpy.shape(rve_macro_strain))

    nr_timesteps = numpy.shape(rve_interpolation_params)[0]
    nr_modes = numpy.shape(rve_interpolation_params)[1]
    logger.debug("Number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("Number of modes detected: {}".format(nr_modes))
    strain_modes = numpy.load(args.strain_modes)[:, :nr_modes]

    filename = "rve_reconstructed.xdmf"
    meshio.write_points_cells(filename, mesh.points, rve_cells)
    with meshio.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(mesh.points, rve_cells)
        for t in range(nr_timesteps):
            logger.info("Timestep {}".format(t))

            logger.debug("Solving fluctuant displacement")
            displacement = numpy.dot(
                strain_correl[:, :nr_modes], rve_interpolation_params[t, :]
            )
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
            comp = numpy.dot(strain_macro_tensor, mesh.points.T)
            total_displacement = comp.T + displacement

            logger.debug("Solving damage")
            damage_list = []
            r = numpy.dot(r_value_correl, data["r_value"][t])
            #r_in_elem = numpy.reshape(r, (-1, 8)) #TODO: FIX THIS
            r_in_elem = {}
            for elem_id, nr_ips in nr_of_ips.items():
                r_in_elem[elem_id] = r[:nr_ips]
                r = r[nr_ips:]

            strain_global = numpy.dot(strain_modes, rve_interpolation_params[t, :])
            stress_list = []
            for elem_id, nr_ips in nr_of_ips.items():
                C = material_properties[material_elem_map[elem_id]]["C"]
                E = material_properties[material_elem_map[elem_id]]["E"]
                nu = material_properties[material_elem_map[elem_id]]["nu"]
                yield_stress = material_properties[material_elem_map[elem_id]][
                    "yield_stress"
                ]
                inf_yield_stress = material_properties[material_elem_map[elem_id]][
                    "inf_yield_stress"
                ]
                H0 = material_properties[material_elem_map[elem_id]]["H0"]
                H1 = material_properties[material_elem_map[elem_id]]["H1"]
                r0 = yield_stress / math.sqrt(E)
                ip_0 = ip_elem_map[elem_id]
                damage = 0
                stress = [0, 0, 0, 0, 0, 0]
                #for i in range(nr_ips):
                #    r = r_in_elem[elem_id], #i]
                for r in r_in_elem[elem_id]:
                    if r < r0:
                        r = r0
                    d = 1 - q(r, E, yield_stress, inf_yield_stress, H0, H1) / r
                    # stress
                    strain = strain_global[ip_0 : ip_0 + nr_comps] + strain_macro
                    stress_ip = (1 - d) * numpy.dot(C, strain)
                    stress = stress + stress_ip / nr_ips

                    damage += d / nr_ips
                    ip_0 += nr_comps
                damage_list.append(damage)
                stress_list.append(stress)
            element_damage = numpy.array(damage_list).reshape(
                (-1, 1)
            )  # formatting for meshio

            logger.debug("Writing timestep data")
            writer.write_data(
                t,
                point_data={
                    "FLUCTUANT_DISPLACEMENT": numpy.reshape(displacement, (-1, 3)),
                    "TOTAL_DISPLACEMENT": total_displacement,
                },
                cell_data={
                      "DAMAGE": element_damage,
                      "STRESS": stress_list,
                },
            )
