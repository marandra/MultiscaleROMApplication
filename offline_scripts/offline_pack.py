import json
import numpy
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)
from offline_common import Common

"""
TODO: pending description here.
"""


def skip_calculation(filename, flag_reuse):
    try:
        with open(filename):
            flag_exists = True
    except IOError:
        flag_exists = False
    return flag_exists and flag_reuse


def write_json(filename, data_dict):
    with open(filename, "w") as fo:
        json.dump(data_dict, fo, indent=2)


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def get_properties(rve_modelpart, iw_list):
    # read model materials
    out_prop = []
    for l in iw_list:
        elem_id = int(l[0])
        elem = rve_modelpart.GetElement(elem_id)
        prop_id = elem.Properties.Id
        out_prop.append(prop_id)
    return out_prop


def unpack_ip_data(iw_list):
    out_e = []
    out_ip = []
    out_w = []
    out_gip = []
    for l in iw_list:
        out_e.append(int(l[0]))
        out_ip.append(int(l[1]))
        out_w.append(float(l[2]))
        out_gip.append(int(l[3]))
    return out_e, out_ip, out_w, out_gip


def parse_strain_bases(strain_bases_filename, iw_list, nr_modes):
    strain_bases = numpy.load(strain_bases_filename, mmap_mode="r")
    strain_bases = strain_bases[:, :nr_modes]
    nr_comps = 6
    out_B = []
    for l in iw_list:
        gip = int(l[3])
        index = gip * nr_comps
        B = strain_bases[index : index + nr_comps, :]
        out_B.append(B.tolist())
    return out_B


def create_rve_params_structure(
    strain_bases_filename,
    rve_materials_filename,
    nr_modes,
    reduced_ip_set,
    rve_modelpart,
):
    """ gather and pack IP data for RVE constitutive law """
    rve_params = {}
    out_e, out_lip, out_w, out_gip = unpack_ip_data(reduced_ip_set)
    # required data
    rve_params["ip_global_id"] = out_gip
    rve_params["ip_weight"] = out_w
    rve_params["ip_property_id"] = get_properties(rve_modelpart, reduced_ip_set)
    rve_params["ip_strain_modes"] = parse_strain_bases(
        strain_bases_filename, reduced_ip_set, nr_modes
    )
    rve_params["material_parameters"] = read_json(rve_materials_filename)
    #  metadata
    rve_params["nr_modes"] = nr_modes
    rve_params["nr_reduced_ip"] = len(reduced_ip_set)
    rve_params["ip_element_id"] = out_e  # TODO: check if we use this data
    rve_params["ip_local_id"] = out_lip  # TODO: check if we use this data
    return rve_params


#######################################################################
#######################################################################

if __name__ == "__main__":

    with open("../configuration.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters
    )
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()

    materials_fname = Common().materials_fname
    for p in Common().ip_subsets:
        nr_points = p.GetInt()
        roc_filename = Common().roc_fname(nr_points)
        ip_set = numpy.loadtxt(roc_filename)
        for m in Common().reduced_nr_modes:
            nr_modes = m.GetInt()
            rve_fname = Common().rve_fname(nr_modes, nr_points)
            if skip_calculation(rve_fname, Common().reuse_existing_files):
                print("File {} exists. Skipping calculation".format(rve_fname))
                continue
            print("Generating {}".format(rve_fname))
            rve_params = create_rve_params_structure(
                Common().strain_bases_fname,
                materials_fname,
                nr_modes,
                ip_set,
                modelpart,
            )
            write_json(rve_fname, rve_params)
