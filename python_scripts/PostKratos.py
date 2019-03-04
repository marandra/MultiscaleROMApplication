#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from structural_mechanics_analysis import StructuralMechanicsAnalysis
import compute_bases as bases

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
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()
    elem = rve_modelpart.GetElement(1)
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, rve_modelpart.ProcessInfo)

    # gather global model part info
    nr_nodes = rve_modelpart.NumberOfNodes()
    nr_elems = rve_modelpart.NumberOfElements()
    nr_ips = len(elem.GetIntegrationPoints())
    with open("integration_weight", 'w') as ofile:
        for elem in rve_modelpart.Elements:
            ip_weights = elem.GetValuesOnIntegrationPoints(Kratos.INTEGRATION_WEIGHT, rve_modelpart.ProcessInfo)
            for ip_weight in ip_weights:
                ofile.write("{}\n".format(ip_weight[0]))

    trajectory_filename = "trajectory"
    nr_e_snap_filename = "elastic_timesteps"

    # compute energy bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.STRAIN_ENERGY, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = 21
    nr_inelastic_modes = 100
    bases_fname = "bases_energy.npy"
    snapshot_filename = "snapshot_energy"
    e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
    print(nr_elems, nr_ips, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)
    bases.generate_bases(nr_elems, nr_ips, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)

    # compute strain bases
    ip_data = elem.GetValuesOnIntegrationPoints(Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, rve_modelpart.ProcessInfo)
    nr_strain_components = len(ip_data[0])
    nr_elastic_modes = 6
    nr_inelastic_modes = 100
    bases_fname = "bases_strain.npy"
    snapshot_filename = "snapshot_strain"
    e_files, i_files = bases.list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
    print(nr_elems, nr_ips, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)
    bases.generate_bases(nr_elems, nr_ips, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)
