# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os
import KratosMultiphysics as km
import KratosMultiphysics.SolidMechanicsApplication as som
km.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return MechanicalSolver(main_model_part, custom_settings)


class MechanicalSolver(object):
    def __init__(self, main_model_part, custom_settings): 
        self.main_model_part = main_model_part    
        default_settings = km.Parameters("""
        {
            "solver_type": "solid_mechanics_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "analysis_type": "Non-Linear",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
    def Initialize(self):
        raise Exception("please implement the Custom Initialization of your solver")


    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(km.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(km.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(km.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(km.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(km.INTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(km.EXTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(km.CONTACT_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(km.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(km.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(som.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(som.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(som.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(km.VOLUME_ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(km.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(km.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(km.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(km.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(km.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(som.PRESSURE_REACTION)
        print("[Solver] Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2;

    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            node.AddDof(km.DISPLACEMENT_X, km.REACTION_X);
            node.AddDof(km.DISPLACEMENT_Y, km.REACTION_Y);
            node.AddDof(km.DISPLACEMENT_Z, km.REACTION_Z);
        if self.settings["rotation_dofs"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(km.ROTATION_X, km.TORQUE_X);
                node.AddDof(km.ROTATION_Y, km.TORQUE_Y);
                node.AddDof(km.ROTATION_Z, km.TORQUE_Z);
        if self.settings["pressure_dofs"].GetBool():                
            for node in self.main_model_part.Nodes:
                node.AddDof(km.PRESSURE, som.PRESSURE_REACTION);
        print("[Solver] DOF's ADDED")

    def ImportModelPart(self):
        print("[Solver] Model reading starts.")
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            km.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("[Solver] Import input model part.")
            # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            aux_params = km.Parameters("{}")
            aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            # CheckAndPrepareModelProcess creates the solid_computational_model_part
            import check_and_prepare_model_process_solid
            check_and_prepare_model_process_solid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

            # Constitutive law import
            constitutive_law = __import__(self.settings["model_import_settings"]["materials_filename"].GetString())
            constitutive_law.AssignMaterial(self.main_model_part.Properties);
            print("[Solver] Constitutive law initialized.")

            self.main_model_part.SetBufferSize(self.settings["buffer_size"].GetInt())
            current_buffer_size = self.main_model_part.GetBufferSize()
            if(self.GetMinimumBufferSize() > current_buffer_size):
                current_buffer_size = self.GetMinimumBufferSize()
            self.main_model_part.SetBufferSize( current_buffer_size )
            #fill buffer
            delta_time = self.main_model_part.ProcessInfo[km.DELTA_TIME]
            time = self.main_model_part.ProcessInfo[km.TIME]
            time = time - delta_time * (current_buffer_size)
            self.main_model_part.ProcessInfo.SetValue(km.TIME, time)            
            for size in range(0, current_buffer_size):
                step = size - (current_buffer_size -1)
                self.main_model_part.ProcessInfo.SetValue(km.STEP, step)
                time = time + delta_time
                #delta_time is computed from previous time in process_info
                self.main_model_part.CloneTimeStep(time)
            self.main_model_part.ProcessInfo[km.IS_RESTARTED] = False
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "rest"):
            problem_path = os.getcwd()
            restart_path = os.path.join(problem_path, self.settings["model_import_settings"]["input_filename"].GetString() + "_" + self.settings["model_import_settings"]["input_file_label"].GetString() )
            if(os.path.exists(restart_path+".rest") == False):
                print("    rest file does not exist , check the restart step selected ")
            # set serializer flag
            self.serializer_flag = "SERIALIZER_NO_TRACE"      # binary
            # self.serializer_flag = "SERIALIZER_TRACE_ERROR" # ascii
            # self.serializer_flag = "SERIALIZER_TRACE_ALL"   # ascii
            kratos_serializer_variable = km.KratosGlobals.GetVariable(self.serializer_flag)
            serializer = Serializer(restart_path, kratos_serializer_variable)
            serializer.Load(self.main_model_part.GetModelPartName(), self.main_model_part)
            print("    Load input restart file.")
            self.main_model_part.ProcessInfo[km.IS_RESTARTED] = True
        else:
            raise Exception("Other input options are not yet implemented.")
        print ("[Solver] Model reading finished.")
        
    def GetComputeModelPart(self):
        return self.main_model_part.GetSubModelPart("solid_computational_model_part")
        
    def GetOutputVariables(self):
        pass
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        self.mechanical_solver.Solve()

    def SetEchoLevel(self, level):
        self.mechanical_solver.SetEchoLevel(level)

    def Clear(self):
        self.mechanical_solver.Clear()
        
    def Check(self):
        self.mechanical_solver.Check()
        
    #### Specific internal functions ####
    def _GetSolutionScheme(self, scheme_type, component_wise, compute_contact_forces):
        raise Exception("please implement the Custom Choice of your Scheme (_GetSolutionScheme) in your solver")
    
    def _GetConvergenceCriterion(self):
        # Creation of an auxiliar Kratos parameters object to store the convergence settings
        conv_params = km.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        
        # Construction of the class convergence_criterion
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_params)
        return convergence_criterion.mechanical_convergence_criterion

    def _GetBuilderAndSolver(self, component_wise, block_builder):
        # Creating the builder and solver
        if(component_wise):
            builder_and_solver = som.ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(block_builder):
                # To keep matrix blocks in builder
                builder_and_solver = km.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = km.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)
        
        return builder_and_solver
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search):
        raise Exception("please implement the Custom Choice of your Mechanical Solver (_GetMechanicalSolver) in your solver")
