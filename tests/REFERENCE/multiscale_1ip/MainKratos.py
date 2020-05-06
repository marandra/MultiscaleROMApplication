# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import sys
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":


    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        sys.exit("Missing parameters filename.")

    with open(fname, "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
