import os
import numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.MultiscaleROMApplication import (
    compute_bases,
    compute_ip_weights,
    io_utilities,
    pack_reduced_rve_dataset,
)
import h5py

"""
TODO: pending description here.
"""


def check_consistent_config_values(config):
    # TODO: Use this function
    # Ideas:
    # number of base modes < number of snapshots
    # number of base mode > number of requested modes
    check = True
    return check


def skip_calculation(filename, flag_reuse):
    try:
        with open(filename):
            flag_exists = True
    except IOError:
        flag_exists = False
    return flag_exists and flag_reuse


def create_bases(field_name, nr_elastic_modes, nr_inelastic_modes, trajectory_filename, svd_algorithm):
    bases_fname = "bases_{}_{}m.npy".format(
        field_name, nr_elastic_modes + nr_inelastic_modes
    )
    print("Generating {} bases".format(field_name))
    e_files, i_files = compute_bases.list_of_snapshots(trajectory_filename, field_name)
    compute_bases.generate_bases(
        nr_elastic_modes,
        nr_inelastic_modes,
        e_files,
        i_files,
        bases_fname,
        svd_algorithm=svd_algorithm,
    )
    try:
        os.rename(
            "singular_values_elastic.dat",
            "sv_{}_elastic_{}.dat".format(field_name, nr_elastic_modes),
        )
        os.rename(
            "singular_values_inelastic.dat",
            "sv_{}_inelastic_{}.dat".format(field_name, nr_inelastic_modes),
        )
    except FileNotFoundError:
        pass
    try:
        os.rename(
            "singular_values.dat",
            "sv_{}_{}.dat".format(field_name, nr_elastic_modes + nr_inelastic_modes),
        )
    except FileNotFoundError:
        pass


#######################################################################
#######################################################################

if __name__ == "__main__":

    with open("1_ProjectParameters.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    config = parameters["config_data"]
    config_defaults = KratosMultiphysics.Parameters(
        """{
    "reuse_existing_files": true,
    "svd_algorithm": "standard",
    "rve_materials_filename": "../training/materials.json",
    "trajectory_filename": "../training/trajectory",
    "nr_elastic_modes_energy": 21,
    "nr_inelastic_modes_energy": 200,
    "nr_elastic_modes_strain": 6,
    "nr_inelastic_modes_strain": 100,
    "rve_data_points": [200, -1],
    "rve_data_modes": [20]
    }
    """
    )
    config.ValidateAndAssignDefaults(config_defaults)

    if not check_consistent_config_values(config):
        exit()

    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters
    )
    simulation.Initialize()
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()

    #
    # gather global model part info
    #
    ip_weights = []
    ip_lids = []
    elem_ids = []
    for elem in rve_modelpart.Elements:
        iw_list = elem.GetValuesOnIntegrationPoints(
            Kratos.INTEGRATION_WEIGHT, rve_modelpart.ProcessInfo
        )
        for ip_lid, ip_weight in enumerate(iw_list):
            ip_weights.append(ip_weight[0])
            ip_lids.append(ip_lid)
            elem_ids.append(elem.Id)
    ip_data = [ip_weights, ip_lids, elem_ids]
    nr_ips = len(ip_data[0])

    svd_algorithm = config["svd_algorithm"].GetString()
    trajectory_filename = config["trajectory_filename"].GetString()

    #
    # compute energy bases
    #
    field_name = "ENERGY_FREE"
    nr_elastic_modes = config["nr_elastic_modes_energy"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_energy"].GetInt()
    bases_fname = "bases_{}_{}m.npy".format(
        field_name, nr_elastic_modes + nr_inelastic_modes
    )
    energy_bases_fname = bases_fname
    if not skip_calculation(bases_fname, config["reuse_existing_files"].GetBool()):
        create_bases(field_name, nr_elastic_modes, nr_inelastic_modes, trajectory_filename, svd_algorithm)
    else:
        print("File {} exists. Skipping calculation".format(bases_fname))

    #
    # compute strain bases
    #
    field_name = "STRAIN_FLUCTUANT"
    nr_elastic_modes = config["nr_elastic_modes_strain"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_strain"].GetInt()
    bases_fname = "bases_{}_{}m.npy".format(
        field_name, nr_elastic_modes + nr_inelastic_modes
    )
    strain_bases_fname = bases_fname
    if not skip_calculation(bases_fname, config["reuse_existing_files"].GetBool()):
        create_bases(field_name, nr_elastic_modes, nr_inelastic_modes, trajectory_filename, svd_algorithm)
    else:
        print("File {} exists. Skipping calculation".format(bases_fname))

    #
    # compute DAMAGE bases
    #
    field_name = "DAMAGE"
    nr_elastic_modes = 0
    nr_inelastic_modes = 1400
    bases_fname = "bases_{}_{}m.npy".format(
        field_name, nr_elastic_modes + nr_inelastic_modes
    )
    if not skip_calculation(bases_fname, config["reuse_existing_files"].GetBool()):
        create_bases(field_name, nr_elastic_modes, nr_inelastic_modes, trajectory_filename, svd_algorithm)
    else:
        print("File {} exists. Skipping calculation".format(bases_fname))

    #
    # compute ip set
    #
    for p in config["rve_data_points"]:
        nr_roc_points = p.GetInt()
        if nr_roc_points != -1:  # HPROM case
            set_name = "{}".format(nr_roc_points)
        else:  # ROM case
            set_name = "{}".format("ROM")
        roc_filename = "roc_{}ip".format(set_name)
        if skip_calculation(roc_filename, config["reuse_existing_files"].GetBool()):
            print("File {} exists. Skipping calculation".format(roc_filename))
            continue
        print("Generating {}".format(roc_filename))
        # compute ROC list
        if nr_roc_points != -1:  # HPROM case
            roc_list = compute_ip_weights.compute_hprom_weights(
                ip_data, nr_roc_points, energy_bases_fname
            )
        else:  # ROM case
            roc_list = compute_ip_weights.compute_rom_weights(ip_data)
        with open(roc_filename, "w") as ofile:
            for list in roc_list:
                ofile.write("{} {} {} {}\n".format(list[0], list[1], list[2], list[3]))

    #
    # pack dataset
    #
    rve_materials_filename = config["rve_materials_filename"].GetString()
    for p in config["rve_data_points"]:
        nr_roc_points = p.GetInt()
        if nr_roc_points != -1:  # HPROM case
            set_name = "{}".format(nr_roc_points)
        else:  # ROM case
            set_name = "{}".format("ROM")
        ip_set = numpy.loadtxt("roc_{}ip".format(set_name))
        for m in config["rve_data_modes"]:
            nr_modes = m.GetInt()
            rve_data_filename = "rve_{}m_{}ip.json".format(nr_modes, set_name)
            if skip_calculation(
                rve_data_filename, config["reuse_existing_files"].GetBool()
            ):
                print("File {} exists. Skipping calculation".format(rve_data_filename))
                continue
            print("Generating {}".format(rve_data_filename))
            rve_params = pack_reduced_rve_dataset.create_rve_params_structure(
                strain_bases_fname,
                rve_materials_filename,
                nr_modes,
                ip_set,
                rve_modelpart,
            )
            io_utilities.write_json(rve_data_filename, rve_params)
