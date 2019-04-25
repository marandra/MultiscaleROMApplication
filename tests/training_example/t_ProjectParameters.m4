{
    "problem_data": {
        "problem_name": "High_Fidelity",
        "parallel_type": "OpenMP",
        "start_time": 0.0,
        "end_time": 0.99,
        "echo_level": 0
    },
    "solver_settings": {
        "model_part_name": "Microstructure",
        "domain_size": 3,
        "echo_level": 1,
        "time_stepping": {
            "time_step": 0.025
        },
        "solver_type": "Static",
        "analysis_type": "non_linear",
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": "model"
        },
        "material_import_settings": {
            "materials_filename" : "materials.json"
        },
        "line_search": false,
        "convergence_criterion": "displacement_criterion",
        "displacement_relative_tolerance": 1e-6,
        "displacement_absolute_tolerance": 0.0,
        "residual_relative_tolerance": 1e-4,
        "residual_absolute_tolerance": 0.0,
        "max_iteration": 10,
        "linear_solver_settings": {
            "solver_type": "ExternalSolversApplication.pastix",
            "solution_method": "Iterative",
            "tolerance": 1e-6,
            "gmres_krylov_space_dimension": 100,
            "ilu_level_of_fill": 1,
            "is_symmetric": true,
            "verbosity": 0,
            "scaling": false,
            "block_size": 10000,
            "use_block_matrices_if_possible": true
        },
        "problem_domain_sub_model_part_list": ["RVE"],
        "processes_sub_model_part_list": ["RVE", "PINNED", "FIXED_X", "FIXED_Y", "FIXED_Z", "MATRIX", "INCLUSION", "COHESIVE"],
        "rotation_dofs": false,
        "compute_reactions": false,
        "move_mesh_flag": false,
        "block_builder": true,
        "auxiliary_variables_list": ["LAGRANGE_MULTIPLIER_1", "LAGRANGE_MULTIPLIER_2", "LAGRANGE_MULTIPLIER_3", "LAGRANGE_MULTIPLIER_4", "LAGRANGE_MULTIPLIER_5", "LAGRANGE_MULTIPLIER_6", "LAGRANGE_DISPLACEMENT"]
        },
    "processes": {
        "my_processes": [{
            "python_module": "lagrange_multiplier_process",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "LagrangeMultiplierProcess",
            "Parameters":{
                "model_part_name": "RVE"
                    }
            },{
            "python_module": "calculate_total_displacement_process",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "ComputeTotalDisplacementProcess",
            "Parameters":{
                "model_part_name": "RVE"
                    }
            },{
            "python_module": "write_integration_weights",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteIntegrationWeights",
            "Parameters": {
                "model_part_name": "RVE",
                "filename": "integration_weight"
                }
            },{
            "python_module": "write_flag_timesteps",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteElementOutputScalar",
            "Parameters": {
                "model_part_name": "RVE",
                "filename": "elastic_timesteps",
                "flag_name": "INELASTIC_FLAG",
                "flag_location": "StructuralMechanicsApplication"
                }
            },{
            "python_module": "write_elements_homogenized_output",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteElementsHomogenizedOutput",
            "Parameters": {
                "model_part_name": "RVE",
                "filename": "homogenized_stress.dat",
                "variable_name": "CAUCHY_STRESS_VECTOR"
                }
            },{
            "python_module": "write_snapshot_energy",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteSnapshotEnergy",
            "Parameters": {
                "model_part_name": "RVE",
                "filename": "snapshot_energy",
                "write_mode": "binary"
                }
            },{
            "python_module": "write_snapshot_strain",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteSnapshotStrain",
            "Parameters": {
                "model_part_name": "RVE",
                "filename": "snapshot_strain",
                "write_mode": "binary"
                }
            } 
        ],
        "list_initial_processes": [],
        "list_boundary_processes": [{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "PINNED",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": [0.0, 0.0, 0.0],
                "interval": [0.0, "End"]
            }
        },{
            "python_module": "assign_scalar_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignScalarVariableProcess",
            "Parameters": {
                "model_part_name": "FIXED_X",
                "variable_name": "DISPLACEMENT_X",
                "constrained": true,
                "value": 0.0,
                "interval": [0.0, "End"]
                }
            },{
            "python_module": "assign_scalar_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignScalarVariableProcess",
            "Parameters": {
                "model_part_name": "FIXED_Y",
                "variable_name": "DISPLACEMENT_Y",
                "constrained": true,
                "value": 0.0,
                "interval": [0.0, "End"]
                }
            },{
            "python_module": "assign_scalar_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignScalarVariableProcess",
            "Parameters": {
                "model_part_name": "FIXED_Z",
                "variable_name": "DISPLACEMENT_Z",
                "constrained": true,
                "value": 0.0,
                "interval": [0.0, "End"]
                }
            }],
        "loads_process_list": [{
            "python_module": "impose_initial_strain_process",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "ImposeInitialStrainProcess",
            "Parameters":{
                "mesh_id": 0,
                "model_part_name": "RVE",
                "variable_name": "INITIAL_STRAIN",
                    "initial_strain": M4VAR_INITIALSTRAIN,
                    "lookuptable_time": [0.0, 1.0],
                    "lookuptable_mult": [0.0, 1.0]
                    }
            }]},
    "output_configuration": {
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostAscii",
                "WriteDeformedMeshFlag": "WriteDeformed",
                "WriteConditionsFlag": "WriteConditions",
                "MultiFileFlag": "SingleFile"
                },
            "file_label": "step",
            "output_control_type": "time",
            "output_frequency": 0,
            "body_output": true,
            "node_output": true,
            "skin_output": false,
            "plane_output": [],
            "nodal_results": ["DISPLACEMENT", "LAGRANGE_DISPLACEMENT"],
            "gauss_point_results": ["GREEN_LAGRANGE_STRAIN_TENSOR", "CAUCHY_STRESS_TENSOR", "STRAIN_ENERGY", "PLASTIC_STRAIN"]
            },
        "point_data_configuration": []
        },
    "restart_options": {
        "SaveRestart": false,
        "RestartFrequency": 0,
        "LoadRestart": false,
        "Restart_Step": 0
    },
    "constraints_data": {
        "incremental_load": false,
        "incremental_displacement": false
    }
}
