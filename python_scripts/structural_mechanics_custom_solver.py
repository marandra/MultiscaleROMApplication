from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library and applications
import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication as MultiscaleROMApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return StaticMechanicalSolver(model, custom_settings)


class StaticMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(StaticMechanicalSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[Custom MechanicalSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        builder_and_solver = MultiscaleROMApplication.ResidualBasedBlockBuilderAndSolverCustom(linear_solver)
        return builder_and_solver
