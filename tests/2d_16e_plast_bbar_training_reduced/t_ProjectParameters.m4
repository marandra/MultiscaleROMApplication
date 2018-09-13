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
        "displacement_relative_tolerance": 1e-4,
        "displacement_absolute_tolerance": 1e-9,
        "residual_relative_tolerance": 1e-4,
        "residual_absolute_tolerance": 1e-9,
        "max_iteration": 10,
        "linear_solver_settings": {
            "solver_type": "SuperLUSolver",
            "scaling": false,
            "verbosity": 0
        },
        "problem_domain_sub_model_part_list": ["RVE"],
        "processes_sub_model_part_list": ["DISPLACEMENT_BC", "RVE", "MATRIX", "INCLUSION"],
        "rotation_dofs": false,
        "compute_reactions": false,
        "move_mesh_flag": false,
        "block_builder": true,
        "auxiliary_variables_list": ["LAGRANGE_DISPLACEMENT"],
        "auxiliary_dofs_list": []
        },
    "constraints_process_list": [{
        "python_module": "assign_vector_variable_process",
        "kratos_module": "KratosMultiphysics",
        "process_name": "AssignVectorVariableProcess",
        "Parameters": {
            "model_part_name": "DISPLACEMENT_BC",
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
            "model_part_name": "RVE",
            "variable_name": "INITIAL_STRAIN",
            "initial_strain": M4VAR_INITIALSTRAIN,
            "lookuptable_time": [0.0, 1.0],
            "lookuptable_mult": [0.0, 1.0]
            }
        },{
        "python_module": "write_elements_output",
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutput",
        "Parameters": {
            "model_part_name": "RVE",
            "filename": "energy",
            "variable_name": "STRAIN_ENERGY",
            "write_mode": "binary"
            }
        },{                                                                     
        "python_module": "write_elements_output",
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutput",                            
        "Parameters": {                                                         
            "model_part_name": "RVE",                                           
            "filename": "strain",                                        
            "variable_name": "GREEN_LAGRANGE_STRAIN_VECTOR",
            "write_mode": "binary",
	    "variable_location": "MultiscaleROMApplication"
            }                                                                   
        },{
        "python_module": "write_elements_output",
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutputScalar",
        "Parameters": {
            "model_part_name": "RVE",
            "filename": "integration_weight",
            "variable_name": "INTEGRATION_WEIGHT",
            "write_frequency": "last_timestep"
            }
        },{
        "python_module": "write_flag_timesteps",
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteGlobalOutputScalarApplication",
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
        "python_module": "calculate_total_displacement_process",
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "ComputeTotalDisplacementProcess",
        "Parameters": {
            "model_part_name": "RVE"
	    }
        }],
    "output_configuration": {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostAscii",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
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
        "point_data_configuration"  : []
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
