# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 
from KratosMultiphysics import *
#from KratosMultiphysics.SolidMechanicsApplication import *
import KratosMultiphysics.MultiscaleROMApplication as msr

def AssignMaterial(Properties):

    #mat = LinearElastic3DLaw()
    mat = msr.LinearIsotropicDamagePlaneStrain2DLaw()
    Properties[1].SetValue(CONSTITUTIVE_LAW, mat.Clone())
    Properties[2].SetValue(CONSTITUTIVE_LAW, mat.Clone())
