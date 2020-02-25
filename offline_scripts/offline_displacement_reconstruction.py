import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import (
    MechanicalSolver,
)


class DisplacementReconstructionSolver(MechanicalSolver):
    def __init__(self, model, custom_settings):
        super(DisplacementReconstructionSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(
            "::[Displacement Reconstruction Solver]:: ", "Construction finished"
        )

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        builder_and_solver = KratosMultiphysics.MultiscaleROMApplication.ResidualBasedBlockBuilderAndSolverCustom(
            linear_solver
        )
        return builder_and_solver


class DisplacementReconstructionAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        super(DisplacementReconstructionAnalysis, self).__init__(
            model, project_parameters
        )

    #### Must be defined ####
    def _CreateSolver(self):
        solver = DisplacementReconstructionSolver(
            self.model, self.project_parameters["solver_settings"]
        )
        return solver


if __name__ == "__main__":

    with open("../configuration_offline_reconstruction.json", "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = DisplacementReconstructionAnalysis(model, parameters)
    # simulation.Run()

    # we replace .Run() by the code below so we can remove conditions
    # (and in the future replace elements, no we don need to modify model.mdpa)
    simulation.Initialize()
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()

    for condition in rve_modelpart.Conditions:
        condition.Set(KratosMultiphysics.TO_ERASE)
    rve_modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

    # settings = KratosMultiphysics.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    # KratosMultiphysics.ReplaceElementsAndConditionsProcess(rve_modelpart, settings).Execute()

    simulation.RunSolutionLoop()
    simulation.Finalize()
