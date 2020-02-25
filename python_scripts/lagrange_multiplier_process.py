# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication as kms


def Factory(settings, Model):
    if type(settings) != KratosMultiphysics.Parameters:
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return LagrangeMultiplierProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class LagrangeMultiplierProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """ 
        {
            "model_part_name": "default"
        }
        """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        model_part = Model[settings["model_part_name"].GetString()]
        self.lagrange_utility = kms.LagrangeMultiplierUtility(model_part)

    def ExecuteInitialize(self):
        self.lagrange_utility.Execute()

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def Clear(self):
        pass
