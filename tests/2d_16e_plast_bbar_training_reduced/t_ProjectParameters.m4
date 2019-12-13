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
        "domain_size": 2,
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
        "convergence_criterion": "residual_criterion",
        "residual_relative_tolerance": 1e-2,
        "residual_absolute_tolerance": 0.0,
        "max_iteration": 10,
        "rotation_dofs": false,
        "compute_reactions": false,
        "move_mesh_flag": false,
        "auxiliary_variables_list": ["LAGRANGE_DISPLACEMENT"]
        },
    "processes": {
        "my_processes": [ ],
        "list_initial_processes": [],
        "list_boundary_processes": [{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "Microstructure.DISPLACEMENT_BC",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": [0.0, 0.0, 0.0],
                "interval": [0.0, "End"]
                }
        }],
        "loads_process_list": [{
            "python_module": "impose_initial_strain_process",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "ImposeInitialStrainProcess",
            "Parameters": {
                "model_part_name": "Microstructure.RVE",
                "variable_name": "INITIAL_STRAIN",
                "initial_strain": M4VAR_INITIALSTRAIN,
                "lookuptable_time": [0.0, 1.0],
                "lookuptable_mult": [0.0, 1.0]
                }
            },{
            "python_module": "calculate_total_displacement_process",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "ComputeTotalDisplacementProcess",
            "Parameters": {
                "model_part_name": "Microstructure.RVE"
	    }
        }]
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
