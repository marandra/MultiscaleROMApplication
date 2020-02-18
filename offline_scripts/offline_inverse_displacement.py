#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.MultiscaleROMApplication.multiscale_rom_analysis import StructuralMechanicsAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("2_ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    #simulation.Run()

    # we replace .Run() by the code below so we can remove conditions
    # (and in the future replace elements, no we don need to modify model.mdpa)
    simulation.Initialize()
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()

    for condition in rve_modelpart.Conditions:
        condition.Set(KratosMultiphysics.TO_ERASE)
    rve_modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

    #settings = KratosMultiphysics.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    #KratosMultiphysics.ReplaceElementsAndConditionsProcess(rve_modelpart, settings).Execute()

    simulation.RunSolutionLoop()
    simulation.Finalize()

