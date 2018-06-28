# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication
import KratosMultiphysics.StructuralMechanicsApplication

from structural_mechanics_analysis import StructuralMechanicsAnalysis

#class StructuralMechanicsAnalysisMultiscale(StructuralMechanicsAnalysis):
#    """This class prints information abt the computen Eigenvectors
#    It also shows how the "StructuralMechanicsAnalysis" can be customized by deriving from it
#    """
#    pass

if __name__ == "__main__":

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    #simulation = StructuralMechanicsAnalysisMultiscale(model,parameters)
    simulation = StructuralMechanicsAnalysis(model,parameters)
    simulation.Run()
