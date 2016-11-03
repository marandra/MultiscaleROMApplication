# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import operator
# import os
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
import KratosMultiphysics.SolidMechanicsApplication as sol
import KratosMultiphysics.MultiScaleApplication as mss
# import KratosMultiphysics.ExternalSolversApplication
# import KratosMultiphysics.StructuralMechanicsApplication
# import TK_Props
import process_factory
import TK_Rve_V2
import TK_GiD
km.CheckForPreviousImport()


def Factory(settings, Model):
    return SolveRVE(settings["Parameters"], Model)


def analysis(parameters, processes, solver, model_part):
    for process in processes:
        process.ExecuteBeforeSolutionLoop()
    delta_time = parameters["problem_data"]["time_step"].GetDouble()
    time = parameters["problem_data"]["start_time"].GetDouble()
    end_time = parameters["problem_data"]["end_time"].GetDouble()
    while(time <= end_time):
        time = time + delta_time
        model_part.CloneTimeStep(time)
        for process in processes:
            process.ExecuteInitializeSolutionStep()
        solver.Solve()
        for process in processes:
            process.ExecuteFinalizeSolutionStep()
        for process in processes:
            process.ExecuteBeforeOutputStep()
        for process in processes:
            process.ExecuteAfterOutputStep()
    for process in processes:
        process.ExecuteFinalize()


def create_model(parameters):
    domain_size = parameters["problem_data"]["domain_size"].GetInt()
    model_part_name = parameters["problem_data"]["part_name"].GetString()
    model_part = km.ModelPart(model_part_name)
    model_part.ProcessInfo.SetValue(km.DOMAIN_SIZE, domain_size)
    Model = {model_part_name: model_part}
    return Model


def create_solver_complete_model_part(model_part, parameters):
    solver_module = __import__(parameters["solver_settings"]["solver_type"].GetString())
    solver = solver_module.CreateSolver(model_part, parameters["solver_settings"])
    solver.AddVariables()
    solver.ImportModelPart()
    solver.AddDofs()
    constitutive_law_name = parameters["solver_settings"]["model_import_settings"]["constitutive_law"].GetString()
    aux_obj_getter = operator.methodcaller(constitutive_law_name)
    model_part.Properties[1].SetValue(km.CONSTITUTIVE_LAW, aux_obj_getter(sol))
    return solver, model_part


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class SolveRVE(km.Process):

    class RVEStrainSize:
        RVE_PLANE_STRESS = 0
        RVE_PLANE_STRAIN = 1
        RVE_3D = 2
        RVE_THERMAL_PLANE_STRESS = 3
        RVE_THERMAL_3D = 4

    def __init__(self, params, MainModel):
#                 MicroModelPart,
#                 StrainSize,
#                 ResultsIOClass,
#                 ResultsOnNodes = [], 
#                 ResultsOnGaussPoints = [],
#                 RveConstraintHandlerClass = RveConstraintHandler_ZBF_SD,
#                 RveHomogenizerClass = RveHomogenizer,
#                 SchemeClass = RveStaticScheme,
#                 LinearSolverClass = SuperLUSolver,
#                 MaxIterations = 10,
#                 CalculateReactions = False,
#                 ReformDofSetAtEachIteration = False,
#                 MoveMesh = False,
#                 ConvergenceCriteriaClass = ResidualNormCriteria,
#                 ConvergenceRelativeTolerance = 1.0E-6,
#                 ConvergenceAbsoluteTolerance = 1.0E-9,
#                 ConvergenceIsVerbose = False,
#                 TargetElementList = [],
#                 OutputElementList = [],
#                 BoundingPolygonNodesID = None,
#                 # NEW
#                 SecondaryRveModeler = None,
#                 IsSecondary = False):

        main_model_part_name = params['problem_data']['main_part_name'].GetString()
        microscale_model_part_name = params['problem_data']['microscale_part_name'].GetString()
        self.microscale_model_part = MainModel[main_model_part_name].GetSubModelPart(microscale_model_part_name)
        self.params = params
        self.model = create_model(params)
        self.model_part_name = params['problem_data']['part_name'].GetString()
        self.model_part = self.model[self.model_part_name]
        self.solver, self.model_part = create_solver_complete_model_part(self.model_part, params)
        for i in range(params["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = params["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            self.model.update({part_name: self.model_part.GetSubModelPart(part_name)})
        self.processes = process_factory.KratosProcessFactory(self.model)\
            .ConstructListOfProcesses(params["constraints_process_list"])
        self.processes += process_factory.KratosProcessFactory(self.model)\
            .ConstructListOfProcesses(params["loads_process_list"])

        f = operator.attrgetter(params['problem_data']['strain_size'].GetString())
        self.StrainSize = f(self.RVEStrainSize)
        if(self.StrainSize == self.RVEStrainSize.RVE_THERMAL_PLANE_STRESS):
                self.RveAdapterClass = mss.RveThermal2DAdapterV2
                self.RveMaterialClass = mss.RveConstitutiveLawV2Thermal2D
        elif(self.StrainSize == self.RVEStrainSize.RVE_THERMAL_3D):
                self.RveAdapterClass = mss.RveThermal3DAdapterV2
                self.RveMaterialClass = mss.RveConstitutiveLawV2Thermal3D
        elif(self.StrainSize == self.RVEStrainSize.RVE_PLANE_STRESS):
                self.RveAdapterClass = mss.RvePlaneStressAdapterV2
                self.RveMaterialClass = mss.RveConstitutiveLawV2PlaneStress
        elif(self.StrainSize == self.RVEStrainSize.RVE_PLANE_STRAIN):
                raise Exception("Rve Plane Strain Not Yet Implemented")
        else: # RVEStrainSize.RVE_3D
                self.RveAdapterClass = mss.Rve3DAdapterV2
                self.RveMaterialClass = mss.RveConstitutiveLawV23D
        self.BoundingPolygonNodesID = None
        self.RveGeometryDescr = None
        self.ResultsIOClass = TK_GiD.ResultsIO
        self.ResultsOnNodes = [km.DISPLACEMENT, km.REACTION, mss.RVE_FULL_DISPLACEMENT]
        self.ResultsOnGaussPoints = [km.GREEN_LAGRANGE_STRAIN_TENSOR, km.CAUCHY_STRESS_TENSOR]
        self.RveConstraintHandlerClass = mss.RveConstraintHandler_PBF_SD
        self.RveHomogenizerClass       = mss.RveHomogenizer
        self.OutputElementList = [1]
        self.TrackList = {}
        self.Initialized = False
        
 
    def ExecuteInitialize(self):
        self.solver.Initialize()
        print("DEBUG SOLVER RVE")
        print(self.solver.mechanical_solver)
        self.Initialize()
        for process in self.processes:
            process.ExecuteInitialize()
        pass
    
    def ExecuteInitializeSolutionStep(self):
        analysis(self.params, self.processes, self.solver, self.model_part)
        pass

    def ExecuteAfterOutputStep(self):
#       self.__write_output(self.model_part.ProcessInfo[TIME])
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        #t = self.Model.ProcessInfo[km.TIME]
        #if t == self.Model.ProcessInfo[km.END_TIME] or self.__check_write_freq(t):
#        self.write_results() 
        pass
    
    def ExecuteFinalize(self):
#       if(self.Initialized == True):
#           self.__finalize_output()
        pass
                                
   
    def Initialize(self):

        def __generate_rve_constitutive_law():
        # This method generates a new rve constitutive law
        # by cloning the rve model part prototype and creating
        # a new rve constitutive law out of it.
        # This method is meant to be private, do NOT call it explicitly

        
            def port_solver_to_adapter_args():
                scheme = self.solver._GetSolutionScheme(
                    self.solver.settings["analysis_type"].GetString(),
                    self.solver.settings["component_wise"].GetBool(),
                    self.solver.settings["compute_contact_forces"].GetBool())
                print(scheme)
                #scheme = self.solver.mechanical_scheme
                dir(scheme)
                convergence = self.solver.mechanical_convergence_criterion
                linear_solver = self.solver.linear_solver
                #builder_and_solver = self.solver.builder_and_solver
                builder_and_solver = self.solver._GetBuilderAndSolver(
                    False, False)
                return scheme, convergence, linear_solver, builder_and_solver


            modelPartClone = km.ModelPart(self.model_part_name)
            mss.RveCloneModelPart(self.model_part, modelPartClone)
            constraint_handler = self.RveConstraintHandlerClass()
            homogenizer = self.RveHomogenizerClass()
            scheme, convergence, linear_solver, builder_and_solver = port_solver_to_adapter_args()


            strategy_linear = km.ResidualBasedLinearStrategy(
                modelPartClone, scheme, linear_solver,
               #self.solver._GetBuilderAndSolver(False, False),
                False, False, False, False)
            print('DIR STRATEGY')
            dir(strategy_linear)

            strategy_nr = km.ResidualBasedNewtonRaphsonStrategy(
                modelPartClone, scheme, linear_solver,
                convergence,
            #    self.solver._GetBuilderAndSolver(False, False),
                10, False, False, False
           )

            adapter = self.RveAdapterClass() 
            adapter.SetRveData(
                modelPartClone,
                mss.RveMacroscaleData(),
                self.RveGeometryDescr,
                constraint_handler,
                mss.RveLinearSystemOfEquations(linear_solver),
                homogenizer,
                scheme,
                convergence,
                self.solver.mechanical_solver,
                builder_and_solver
                )
        
            rveLaw = self.RveMaterialClass(adapter) # finally generate the constitutive law adapter
            
            for i in range(modelPartClone.GetBufferSize()):
                modelPartClone.CloneTimeStep(0.0)
            return  (rveLaw,modelPartClone) # return a tuple
    

        def __track_rve_constitutive_law(rveLaw, elemID, gpID):
        # This method tracks a rve constitutive law
        # at a given element in a given gauss point.
        # This method is meant to be private, do NOT call it explicitly
                elInfo = TK_Rve_V2.SolidElementInfo(elemID, gpID)
                
                if( next((x for x in self.OutputElementList if x == elemID), None) is not None ):
                        outputFileName = self.model_part.Name + "__" + elInfo.GetStringExtension()
                        rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
                        self.TrackList[elInfo] = (rveLaw, rveLawIO)
                else:
                        self.TrackList[elInfo] = (rveLaw, None)
        
        def __assign_rve_constitutive_law(Element):
        # This method assignes a rve constitutive law
        # at a given element.
        # This method is meant to be private, do NOT call it explicitly
                
                # list of generated rve mdpa clones
                rve_mdpa_clones = []
                
                # get the number of integration points
                elemIntPoints = Element.GetIntegrationPoints()
                num_gp = len(elemIntPoints)
                elem_id = Element.Id
                
                # get a reference to the process into
                pinfo = self.model_part.ProcessInfo
                
                # prepare the list of constitutive laws for the element
                constitutiveLaws = []
                
                # for each element integration point ...
                for gp_id in range(num_gp):
                        
                        # generate a new rve constitutive law
                        rve_law__rve_mdpa__tuple = __generate_rve_constitutive_law()
                        aRveLaw = rve_law__rve_mdpa__tuple[0]
                        constitutiveLaws.append(aRveLaw)
                        
                        # TODO: check what rve law to track...
                        # for the moment let's track them all
                        __track_rve_constitutive_law(aRveLaw, elem_id, gp_id)
                        
                        # store the rve mdpa clone
                        rve_mdpa_clones.append(rve_law__rve_mdpa__tuple[1])
                
                # assign the list of constitutive laws
                Element.SetValuesOnIntegrationPoints(sol.CONSTITUTIVE_LAW_POINTER, constitutiveLaws, pinfo)
                
                return rve_mdpa_clones
        

        def __initialize_output():
        ## Initializes the output for the tracked rves (only if required)
                for key, value in self.TrackList.items():
                        rveIO = value[1]
                        if(rveIO is not None):
                                rveIO.Initialize()
        

        def __write_output(self, currentTime):
        ## Writes the output for the tracked rves (only if required)
                for key, value in self.TrackList.items():
                        rveIO = value[1]
                        if(rveIO is not None):
                                rveIO.Write(currentTime)
        

        def __finalize_output(self):
        ## Finalizes the output for the tracked rves (only if required)
                for key, value in self.TrackList.items():
                        rveIO = value[1]
                        if(rveIO is not None):
                                rveIO.Finalize()
        

        def __print_info(self):
                
                print ("")
                print ("====================================================")
                print ("RveModelerShell - Info:")
                print ("====================================================")
                print ("MODEL PART - PROTOTYPE:")
                print (self.model_part)
                print ("====================================================")
                print ("TRACK LIST:")
                ii = 0
                print ("+--------------------------------------------------------+")
                for key, value in self.TrackList.items():
                        print ("AT[", ii, "]")
                        print ("Info:")
                        print (key)
                        print ("(RveMaterial, IO)")
                        print (value)
                        print ("Micro Model Clone:")
                        micro = value[0].GetModelPart()
                        print (hex(id(micro)))
                        print (micro)
                        print ("+--------------------------------------------------------+")
                        ii+=1

        if(self.Initialized == False):
            # initialize the geometry descriptor
            self.RveGeometryDescr = mss.RveGeometryDescriptor()
            if(self.BoundingPolygonNodesID is not None):
                self.RveGeometryDescr.SetUserCornerNodes(self.BoundingPolygonNodesID)
            self.RveGeometryDescr.Build(self.model_part)
            #print(self.RveGeometryDescr)
            self.stored_rvemdpa_clones=[]

            #for elem_id in self.TargetElementList:
            #    elem = self.microscale_model_part.Elements[elem_id]
            for elem in self.microscale_model_part.Elements:
                elem_rvemdpa_clone_list = __assign_rve_constitutive_law(elem)
                for iclone in elem_rvemdpa_clone_list:
                    self.stored_rvemdpa_clones.append(iclone)
           
            # initialize the output
            __initialize_output()
            
            # set initialization flag
            self.Initialized = True
