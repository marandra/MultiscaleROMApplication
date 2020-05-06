"""
description here
"""
from  pathlib import Path
import json
import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import (
    MechanicalSolver,
)
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)
from offline_common import Common


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


###############################################################
###############################################################

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        co = Common(root_path=Path(sys.argv[1]))
    else:
        exit("Missing root_path argument.")

    # Read parametres for reconstruction
    with open("../ProjectParameters_correlation.json", "r") as parameter_file:
        parameters_reconstr = KratosMultiphysics.Parameters(parameter_file.read())

    #  Generate auxiliar data structure
    parameters_dict = {
        "problem_data": {
        "problem_name": "High_Fidelity",
        "parallel_type": "OpenMP",
        "start_time": 0.0,
        "end_time": 0.99,
        "echo_level": 1,
    },
    "solver_settings": {
        "model_part_name": "Microstructure",
        "domain_size": 3,
        "echo_level": 1,
        "time_stepping": {},
        "solver_type": "Static",
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": "{}/model".format(co.training_path),
        },
        "material_import_settings": {
            "materials_filename": "{}/materials.json".format(co.training_path)
        },
    },
    }
    parameters_aux = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters_aux
    )
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()
    for elem in modelpart.Elements:
        nr_comp = len(
            elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
            )[0]
        )
        break
    idx_vector = []
    count = 0
    for elem in modelpart.Elements:
        idx_vector.append(count)
        nr_ips = len(
            elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
            )
        )
        count = count + nr_ips * nr_comp
    fname = parameters_reconstr["processes"]["my_processes"][0]["Parameters"][
        "global_index_filename"
    ].GetString()
    with open(fname, "w") as ofile:
        for idx in idx_vector:
            ofile.write("{}\n".format(idx))
    # end of generating auxiliar file

    # Reconstruction
    model = KratosMultiphysics.Model()
    simulation = DisplacementReconstructionAnalysis(model, parameters_reconstr)
    # we replace .Run() by the code below so we can remove conditions
    # (and in the future replace elements, no we don need to modify model.mdpa)
    # simulation.Run()
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()

    for condition in modelpart.Conditions:
        condition.Set(KratosMultiphysics.TO_ERASE)
    modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

    # settings = KratosMultiphysics.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    # KratosMultiphysics.ReplaceElementsAndConditionsProcess(modelpart, settings).Execute()

    simulation.RunSolutionLoop()
    simulation.Finalize()
