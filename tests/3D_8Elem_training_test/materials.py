from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.MultiscaleROMApplication import *

def AssignMaterial(Properties):

    #mat = LinearElastic3DLaw()
    mat = SmallDisplacementElastoPlasticJ23DLaw()
    #prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
    Properties[1].SetValue(CONSTITUTIVE_LAW, mat.Clone())
    Properties[2].SetValue(CONSTITUTIVE_LAW, mat.Clone())
