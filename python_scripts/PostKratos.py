#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from structural_mechanics_analysis import StructuralMechanicsAnalysis
#from multiscale_rom_analysis import StructuralMechanicsAnalysis
import compute_bases as bases
import compute_ip_weights as hprom
import pack_reduced_rve_dataset as pack
import compute_stress_reconstruction_system as stress_reconstruction
import numpy
"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("PostProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model,parameters)
    simulation.Initialize()

    # generate strain-displacement correlation matrix
    #simulation.Run()

    rve_modelpart = simulation._GetSolver().GetComputingModelPart()
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
    integration_weights_filename = "integration_weight"
    with open(integration_weights_filename, 'w') as ofile:
        for ip_weight in integration_weights:
            ofile.write("{}\n".format(ip_weight))


    trajectory_filename = "../training/trajectory"
    nr_e_snap_filename = "elastic_timesteps"
    # compute energy bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.STRAIN_ENERGY, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = 21
    nr_inelastic_modes = 200
    energy_bases_fname = "bases_energy_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = "snapshot_energy"
    e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
#    bases.generate_bases(nr_elems, nr_ips_per_elem, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, energy_bases_fname)

    # compute strain bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = 6
    nr_inelastic_modes = 100
    strain_bases_fname = "bases_strain_{}m.npy".format(nr_elastic_modes + nr_inelastic_modes)
    snapshot_filename = "snapshot_strain"
    e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
#    bases.generate_bases(nr_elems, nr_ips_per_elem, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, strain_bases_fname)

    # computed reduces ip set
    nr_roq_points = 100  # TODO: too many makes rank A != size A
    roq_list = hprom.compute_hprom_weights(nr_ips_per_elem, integration_weights, nr_roq_points, energy_bases_fname)
    # optional output
    roq_filename = "roq_{}ip".format(nr_roq_points)
    with open(roq_filename, 'w') as ofile:
        for list in roq_list:
            ofile.write("{} {} {} {}\n".format(list[0], list[1], list[2], list[3]))

    # pack RVE dataset
    rve_mdpa_filename = "../training/model.mdpa"
    reduced_ip_set = numpy.loadtxt(roq_filename)
    nr_modes = 20
    rve_params = pack.create_rve_params_structure(strain_bases_fname, rve_mdpa_filename, nr_modes, reduced_ip_set)
    rve_data_filename = "rve_{}m_{}ip.json".format(nr_modes, nr_roq_points)
    pack.util.write_json(rve_data_filename, rve_params)

    # generate stress reconstruction system
    A = stress_reconstruction.compute_system(rve_data_filename, energy_bases_fname, integration_weights_filename)
    stress_reconstruction.util.write_numpy_file('correlation_stress.npy', 'binary', A)


