# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis as periodic_bc_analysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    with open("ProjectParameters.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = periodic_bc_analysis.PBCAnalysis(model, parameters)
    simulation.Run()
