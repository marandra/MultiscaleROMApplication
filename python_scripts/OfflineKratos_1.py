# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)
import KratosMultiphysics.MultiscaleROMApplication.compute_bases as bases
import KratosMultiphysics.MultiscaleROMApplication.compute_ip_weights as roc
import KratosMultiphysics.MultiscaleROMApplication.pack_reduced_rve_dataset as pack
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities
import KratosMultiphysics.MultiscaleROMApplication.compute_stress_reconstruction_system as stress_reconstruction

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
    "elastic_snapshots_filename": "elastic_timesteps",
    "snapshot_energy_filename": "snapshot_energy",
    "snapshot_strain_filename": "snapshot_strain",
    "bases_energy_filename": "bases_energy",
    "bases_strain_filename": "bases_strain",
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

    # Read the original model part used in training,
    # and replaces elements and conditions needed by offline process
    # Remove conditions from model part
    #   for condition in rve_modelpart.Conditions:
    #      condition.Set(KratosMultiphysics.TO_ERASE)
    #   rve_modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
    #   print("DEBUG ********")
    #   print(rve_modelpart)

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
    nr_e_snap_filename = config["elastic_snapshots_filename"].GetString()
    #
    # compute energy bases
    #
    nr_strain_components = 1
    nr_elastic_modes = config["nr_elastic_modes_energy"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_energy"].GetInt()
    energy_bases_fname = config[
        "bases_energy_filename"
    ].GetString() + "_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = config["snapshot_energy_filename"].GetString()
    print("Generating ENERGY bases")
    if skip_calculation(energy_bases_fname, config["reuse_existing_files"].GetBool()):
        print("File {} exists. Skipping calculation".format(energy_bases_fname))
    else:
        e_files, i_files = bases.list_of_snapshots(
            trajectory_filename, nr_e_snap_filename, snapshot_filename
        )
        bases.generate_bases(
            nr_ips,
            nr_strain_components,
            nr_elastic_modes,
            nr_inelastic_modes,
            e_files,
            i_files,
            energy_bases_fname,
            svd_algorithm=svd_algorithm,
        )
    #
    # compute strain bases
    #
    nr_strain_components = 6
    nr_elastic_modes = config["nr_elastic_modes_strain"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_strain"].GetInt()
    strain_bases_fname = config[
        "bases_strain_filename"
    ].GetString() + "_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = config["snapshot_strain_filename"].GetString()
    print("Generating STRAIN bases")
    if skip_calculation(strain_bases_fname, config["reuse_existing_files"].GetBool()):
        print("File {} exists. Skipping calculation".format(strain_bases_fname))
    else:
        e_files, i_files = bases.list_of_snapshots(
            trajectory_filename, nr_e_snap_filename, snapshot_filename
        )
        bases.generate_bases(
            nr_ips,
            nr_strain_components,
            nr_elastic_modes,
            nr_inelastic_modes,
            e_files,
            i_files,
            strain_bases_fname,
            svd_algorithm=svd_algorithm,
        )
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
            roc_list = roc.compute_hprom_weights(
                ip_data, nr_roc_points, energy_bases_fname
            )
        else:  # ROM case
            roc_list = roc.compute_rom_weights(ip_data)
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
            rve_params = pack.create_rve_params_structure(
                strain_bases_fname,
                rve_materials_filename,
                nr_modes,
                ip_set,
                rve_modelpart,
            )
            io_utilities.write_json(rve_data_filename, rve_params)

    # generate stress reconstruction system
    # A = stress_reconstruction.compute_system(rve_data_filename, energy_bases_fname, integration_weights_filename)
    # stress_reconstruction.util.write_numpy_file('correlation_stress_{}m_{}ip.npy'.format(nr_modes, nr_roc_points), 'binary', A)

    # from multiscale_rom_analysis import StructuralMechanicsAnalysis
    # settings = Kratos.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    # KratosMultiphysics.ReplaceElementsAndConditionsProcess(rve_modelpart, settings).Execute()
    # print("DEBUG ********")
    # print(rve_modelpart)

    # generate strain-displacement correlation matrix
    # simulation.Run()
