#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from  structural_mechanics_analysis import StructuralMechanicsAnalysis
import compute_bases as bases
import compute_ip_weights as hprom
import pack_reduced_rve_dataset as pack
import compute_stress_reconstruction_system as stress_reconstruction
import numpy

"""
Here description.
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

def check_consistent_config_values(config):
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

    with open("1_ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    config = parameters["config_data"]
    config_defaults = KratosMultiphysics.Parameters('''{
    "reuse_existing_files": true,
    "rve_mdpa_filename": "../training/model.mdpa",
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
    "rve_data_points": [200],
    "rve_data_modes": [20]
    }
    ''')
    config.ValidateAndAssignDefaults(config_defaults)
    print("DEBUG: ")
    print(config)

    if not check_consistent_config_values(config):
        exit()

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model,parameters)
    simulation.Initialize()
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()

    # Read the original model part used in training,
    # and replaces elements and conditions needed by offline process
    # Remove conditions from model part
    #for condition in rve_modelpart.Conditions:
    #    condition.Set(KratosMultiphysics.TO_ERASE)
    #rve_modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
    #print("DEBUG ********")
    #print(rve_modelpart)




    elem = rve_modelpart.GetElement(1)
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, rve_modelpart.ProcessInfo)


    # gather global model part info
    nr_nodes = rve_modelpart.NumberOfNodes()
    nr_elems = rve_modelpart.NumberOfElements()
    nr_ips_per_elem = len(elem.GetIntegrationPoints())
    integration_weights = []
    for elem in rve_modelpart.Elements:
        ip_weights = elem.GetValuesOnIntegrationPoints(Kratos.INTEGRATION_WEIGHT, rve_modelpart.ProcessInfo)
        for ip_weight in ip_weights:
            integration_weights.append(ip_weight[0])
    # optional output
    #integration_weights_filename = "integration_weight"
    #with open(integration_weights_filename, 'w') as ofile:
    #    for ip_weight in integration_weights:
    #        ofile.write("{}\n".format(ip_weight))


    trajectory_filename = config["trajectory_filename"].GetString()
    nr_e_snap_filename = config["elastic_snapshots_filename"].GetString()
    # compute energy bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.STRAIN_ENERGY, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = config["nr_elastic_modes_energy"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_energy"].GetInt()
    energy_bases_fname = config["bases_energy_filename"].GetString() + "_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = config["snapshot_energy_filename"].GetString()
    if skip_calculation(energy_bases_fname, config["reuse_existing_files"].GetBool()):
        print("File {} exists. Skipping calculation".format(energy_bases_fname))
    else:
        e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
        bases.generate_bases(nr_elems, nr_ips_per_elem, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, energy_bases_fname)

    # compute strain bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = config["nr_elastic_modes_strain"].GetInt()
    nr_inelastic_modes = config["nr_inelastic_modes_strain"].GetInt()
    strain_bases_fname = config["bases_strain_filename"].GetString() + "_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = config["snapshot_strain_filename"].GetString()
    if skip_calculation(strain_bases_fname, config["reuse_existing_files"].GetBool()):
        print("File {} exists. Skipping calculation".format(strain_bases_fname))
    else:
        e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
        bases.generate_bases(nr_elems, nr_ips_per_elem, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, strain_bases_fname)

    #
    # computed reduced ip set and pack dataset
    #
    rve_mdpa_filename = config["rve_mdpa_filename"].GetString()
    for p in config["rve_data_points"]:
        nr_roq_points = p.GetInt()
        roq_filename = "roq_{}ip".format(nr_roq_points)
        if skip_calculation(roq_filename, config["reuse_existing_files"].GetBool()):
            print("File {} exists. Skipping calculation".format(roq_filename))
            continue
        # compute ROQ list
        roq_list = hprom.compute_hprom_weights(nr_ips_per_elem, integration_weights, nr_roq_points, energy_bases_fname)
        with open(roq_filename, 'w') as ofile:
            for list in roq_list:
                ofile.write("{} {} {} {}\n".format(list[0], list[1], list[2], list[3]))

    for p in config["rve_data_points"]:
        nr_roq_points = p.GetInt()
        roq_filename = "roq_{}ip".format(nr_roq_points)
        for m in config["rve_data_modes"]:
            nr_modes = m.GetInt()
            reduced_ip_set = numpy.loadtxt(roq_filename)  # TODO: use variable instead or reading file
            rve_params = pack.create_rve_params_structure(strain_bases_fname, rve_mdpa_filename, nr_modes, reduced_ip_set)
            rve_data_filename = "rve_{}m_{}ip.json".format(nr_modes, nr_roq_points)
            pack.util.write_json(rve_data_filename, rve_params)

    # generate stress reconstruction system
    #A = stress_reconstruction.compute_system(rve_data_filename, energy_bases_fname, integration_weights_filename)
    #stress_reconstruction.util.write_numpy_file('correlation_stress_{}m_{}ip.npy'.format(nr_modes, nr_roq_points), 'binary', A)

    #from multiscale_rom_analysis import StructuralMechanicsAnalysis
    #settings = Kratos.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    #KratosMultiphysics.ReplaceElementsAndConditionsProcess(rve_modelpart, settings).Execute()
    #print("DEBUG ********")
    #print(rve_modelpart)

    # generate strain-displacement correlation matrix
    #simulation.Run()
