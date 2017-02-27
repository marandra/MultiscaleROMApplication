from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.MultiscaleROMApplication import *


def AssignMaterial(Properties):

    mat = LinearElasticPlasticJ2PlaneStrain2DLaw()
    #mat = LinearIsotropicDamagePlaneStrain2DLaw()
    #mat = LinearElasticPlaneStrain2DLaw()
    #mat = HyperElasticPlaneStrain2DLaw()
    #mat = HyperElasticPlasticJ2PlaneStrain2DLaw()

    Properties[1].SetValue(CONSTITUTIVE_LAW, mat.Clone())
    Properties[2].SetValue(CONSTITUTIVE_LAW, mat.Clone())
        
