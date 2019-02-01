from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from analysis_stage import AnalysisStage

class StructuralMechanicsAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)

    #### Must be defined ####
    def _CreateSolver(self):
        solver_module = __import__("structural_mechanics_custom_solver")
        solver_settings = self.project_parameters["solver_settings"]
        solver = solver_module.CreateSolver(self.model, solver_settings)
        return solver
