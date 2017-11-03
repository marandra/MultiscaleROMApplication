{
    "problem_data": {
        "problem_name": "training",
        "part_name": "DOMAIN",
        "domain_size": 2,
        "nr_time_steps": 250,
        "end_time": 1.0,
        "echo_level": 0
        },
    "solver_settings": {
        "solver_type": "structural_mechanics_static_solver",
        "echo_level": 1,
        "analysis_type": "non_linear",
        "model_import_settings": {
            "input_filename": "model",
            "input_type": "mdpa"
            },
        "material_import_settings": {
            "materials_filename" : "materials.json"
            },
        "line_search": false,
        "convergence_criterion": "residual_criterion",
        "displacement_relative_tolerance": 1e-13,
        "displacement_absolute_tolerance": 1e-13,
        "residual_relative_tolerance": 1e-13,
        "residual_absolute_tolerance": 1e-13,
        "max_iteration": 10,
        "linear_solver_settings": {
            "solver_type": "SuperLUSolver",
            "scaling": false,
            "verbosity": 0
            },
        "problem_domain_sub_model_part_list": ["RVE"],
        "processes_sub_model_part_list": ["DISPLACEMENT_BC", "RVE", "MATRIX", "INCLUSION"],
        "rotation_dofs": false,
        "move_mesh_flag": false
        },
    "constraints_process_list": [{
        "implemented_in_file": "assign_vector_variable_process",
        "implemented_in_module": "KratosMultiphysics",
        "process_name": "AssignVectorVariableProcess",
        "Parameters": {
            "mesh_id": 0,
            "model_part_name": "DISPLACEMENT_BC",
            "variable_name": "DISPLACEMENT",
            "constrained": [true, true, true],
            "value": [0.0, 0.0, 0.0]
            }
    }],
    "loads_process_list": [{
        "implemented_in_file": "impose_initial_strain_process",
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "ImposeInitialStrainProcess",
        "Parameters": {
            "model_part_name": "RVE",
            "variable_name": "INITIAL_STRAIN",
            "initial_strain": M4VAR_INITIALSTRAIN,
            "lookuptable_time": [0.0, 1.0],
            "lookuptable_mult": [0.0, 1.0]
            }
    },{
        "implemented_in_file": "write_elements_output",               
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutput",                            
        "Parameters": {                                                         
            "model_part_name": "RVE",                                           
            "filename": "energy",                                        
            "variable_name": "STRAIN_ENERGY",
	        "write_mode": "binary"
            }                                                                   
        },{                                                                     
        "implemented_in_file": "write_elements_output",
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutput",                            
        "Parameters": {                                                         
            "model_part_name": "RVE",                                           
            "filename": "strain",                                        
            "variable_name": "GREEN_LAGRANGE_STRAIN_VECTOR",
	        "write_mode": "binary"
            }                                                                   
        },{
        "implemented_in_file": "write_elements_output",
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsOutputScalar",
        "Parameters": {
            "model_part_name": "RVE",
	        "filename": "integration_weight",
            "variable_name": "INTEGRATION_WEIGHT",
	        "write_frequency": "last_timestep"
	        }
        },{
        "implemented_in_file": "write_flag_timesteps",
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteGlobalOutputScalarApplication",
        "Parameters": {
            "model_part_name": "RVE",
	        "filename": "elastic_timesteps",
            "flag_name": "INELASTIC_FLAG",
	        "flag_location": "MultiscaleROMApplication"
	        }
        },{
        "implemented_in_file": "write_elements_homogenized_output",
        "implemented_in_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "WriteElementsHomogenizedOutput",
        "Parameters": {
            "model_part_name": "RVE",
	        "filename": "homogenized_stress.dat",
            "variable_name": "CAUCHY_STRESS_VECTOR"
	        }
        }],
    "output_configuration": {},
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
