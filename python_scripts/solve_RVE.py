# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.MultiScaleApplication as mss # <- check is used
import os
import operator
#from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.MultiScaleApplication import *
#from KratosMultiphysics.SolidMechanicsApplication import *
#from KratosMultiphysics.StructuralMechanicsApplication import *
#import TK_Props
km.CheckForPreviousImport()


def Factory(settings, Model):
#    if(type(settings) != km.Parameters):
#        raise Exception("expected input is Parameters object, encapsulating a json string")
    return SolveRVE(settings["Parameters"], Model)







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
 

    def __init__(self, params, Model):
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

         self.model_part = Model[params['model_part_name'].GetString()]
         f = operator.attrgetter(params['strain_size'].GetString())
         self.StrainSize = f(self.RVEStrainSize)
#         self.BoundingPolygonNodesID = BoundingPolygonNodesID
#         self.IsSecondary = IsSecondary
#         if(self.StrainSize == RVEStrainSize.RVE_THERMAL_PLANE_STRESS):
#         	self.RveAdapterClass = RveThermal2DAdapterV2
#         	self.RveMaterialClass = RveConstitutiveLawV2Thermal2D
#         elif(self.StrainSize == RVEStrainSize.RVE_THERMAL_3D):
#         	self.RveAdapterClass = RveThermal3DAdapterV2
#         	self.RveMaterialClass = RveConstitutiveLawV2Thermal3D
#         elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRESS):
#         	self.RveAdapterClass = RvePlaneStressAdapterV2
#         	self.RveMaterialClass = RveConstitutiveLawV2PlaneStress
#         elif(self.StrainSize == RVEStrainSize.RVE_PLANE_STRAIN):
#         	raise Exception("Rve Plane Strain Not Yet Implemented")
#         else: # RVEStrainSize.RVE_3D):
#         	self.RveAdapterClass = Rve3DAdapterV2
#         	self.RveMaterialClass = RveConstitutiveLawV23D
#         self.RveGeometryDescr = None
#         self.ResultsIOClass = ResultsIOClass
#         self.ResultsOnNodes = ResultsOnNodes
#         self.ResultsOnGaussPoints = ResultsOnGaussPoints
#         self.RveConstraintHandlerClass = RveConstraintHandlerClass
#         self.RveHomogenizerClass       = RveHomogenizerClass
#         self.SchemeClass               = SchemeClass
#         self.LinearSolverClass = LinearSolverClass
#         self.MaxIterations = MaxIterations
#         self.CalculateReactions = CalculateReactions
#         self.ReformDofSetAtEachIteration = ReformDofSetAtEachIteration
#         self.MoveMesh = MoveMesh
#         self.ConvergenceCriteriaClass = ConvergenceCriteriaClass
#         self.ConvergenceRelativeTolerance = ConvergenceRelativeTolerance
#         self.ConvergenceAbsoluteTolerance = ConvergenceAbsoluteTolerance
#         self.ConvergenceIsVerbose = ConvergenceIsVerbose
#         self.TargetElementList = TargetElementList
#         self.OutputElementList = OutputElementList
#         self.TrackList = {}
#         self.Initialized = False
#         self.SecondaryRveModeler = SecondaryRveModeler
#         self.IsSecondary         = IsSecondary
         # if(self.IsSecondary == False):
         	# if(self.SecondaryRveModeler is None):
         		# raise exeption(" -- RVEModelerSolid is the first physic and need SecondaryRveModeler for add RVE_Clone_ModelPart_List -- ")
	
 
    def ExecuteInitialize(self):
#        self.Initialize(self.MicroModelPart)
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
#	self.__write_output(self.MicroModelPart.ProcessInfo[TIME])
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
#	if(self.Initialized == True):
#	    self.__finalize_output()
        pass
				
   
    def Initialize(self):
        if(self.Initialized == False):
            # initialize the geometry descriptor
            self.RveGeometryDescr = RveGeometryDescriptor()
            if(self.BoundingPolygonNodesID is not None):
                self.RveGeometryDescr.SetUserCornerNodes(self.BoundingPolygonNodesID)
            self.RveGeometryDescr.Build(self.model_part)
            #print(self.RveGeometryDescr)
            # generate,assign and track all required rve's
            if(self.IsSecondary == True):
                # il primario ha generato la lista di cloni
                # ho bisogno di sapere in quale id mi trovo della lista
                self.clone_list_counter = 0
                for elem_id in self.TargetElementList:
                    elem = model_part.Elements[elem_id]
                    dummy = self.__assign_rve_constitutive_law(elem)
                self.clone_list_counter = 0 # non necessario ma per sicurazzo lo riazzeriamo
            else:
                # se sono il primario genero una lista di [nelem*ngauss] di rve clones...
                self.stored_rvemdpa_clones=[]
                for elem_id in self.TargetElementList:
                    elem = model_part.Elements[elem_id]
                    elem_rvemdpa_clone_list = self.__assign_rve_constitutive_law(elem)
                    for iclone in elem_rvemdpa_clone_list:
                        self.stored_rvemdpa_clones.append(iclone)
                # ... e la copio nel modeler secondario (che non dovra generarla!!!!!)
                if(self.SecondaryRveModeler is not None):
                    self.SecondaryRveModeler.stored_rvemdpa_clones = self.stored_rvemdpa_clones
           
            # initialize the output
            self.__initialize_output()
            
            # set initialization flag
            self.Initialized = True

	# private methods *******************************************************************************************
	
	## __generate_rve_constitutive_law
	#
	# This method generates a new rve constitutive law
	# by cloning the rve model part prototype and creating
	# a new rve constitutive law out of it.
	# This method is meant to be private, do NOT call it explicitly
	def __generate_rve_constitutive_law(self):
	
		if(self.IsSecondary == True):
			current_rve_primary_clone = self.stored_rvemdpa_clones[ self.clone_list_counter ]
			modelPartClone = ModelPart(self.model_part.Name + "_RVE")
			RveCloneModelPart2Physics(self.model_part, current_rve_primary_clone, modelPartClone) # clone the model part prototype
			
			self.clone_list_counter = self.clone_list_counter + 1
		else:
			modelPartClone = ModelPart(self.model_part.Name + "_RVE")
			RveCloneModelPart(self.model_part, modelPartClone) # clone the model part prototype
		
		msData = RveMacroscaleData() 
		
		linSolver = self.LinearSolverClass() 
		
		timeScheme = self.SchemeClass()
		timeScheme.Check(modelPartClone)
		
		convCriteria = self.ConvergenceCriteriaClass(
			self.ConvergenceRelativeTolerance,
			self.ConvergenceAbsoluteTolerance,
			self.ConvergenceIsVerbose,
			)
			
		constraint_handler = self.RveConstraintHandlerClass()
		
		homogenizer = self.RveHomogenizerClass()
		
		adapter = self.RveAdapterClass() # generate the rve adapter
		
		adapter.SetRveData(
			modelPartClone,
			msData,
			self.RveGeometryDescr,
			constraint_handler,
			RveLinearSystemOfEquations(linSolver),
			homogenizer,
			timeScheme,
			convCriteria
		) # set all data (just for testing...)
		
		rveLaw = self.RveMaterialClass(adapter) # finally generate the constitutive law adapter
		
		if (self.IsSecondary == False):
			for i in range(modelPartClone.GetBufferSize()):
				modelPartClone.CloneTimeStep(0.0)
		return  (rveLaw,modelPartClone) # return a tuple
	
	## __track_rve_constitutive_law
	#
	# This method tracks a rve constitutive law
	# at a given element in a given gauss point.
	# This method is meant to be private, do NOT call it explicitly
	def __track_rve_constitutive_law(self, rveLaw, elemID, gpID):
		elInfo = SolidElementInfo(elemID, gpID)
		
		if( next((x for x in self.OutputElementList if x == elemID), None) is not None ):
			outputFileName = self.model_part.Name + "__" + elInfo.GetStringExtension()
			rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints)
			# if (self.IsSecondary == False):
				# print ("ResultsIOClass Mechanical Mdpa")
				# rveLawIO = self.ResultsIOClass(rveLaw.GetModelPart(), self.MicroModelPartB.Model, outputFileName, self.ResultsOnNodes, self.ResultsOnGaussPoints_ModA, self.ResultsOnGaussPoints_ModB)
			self.TrackList[elInfo] = (rveLaw, rveLawIO)
		else:
			self.TrackList[elInfo] = (rveLaw, None)
	
	## __assign_rve_constitutive_law
	#
	# This method assignes a rve constitutive law
	# at a given element.
	# This method is meant to be private, do NOT call it explicitly
	def __assign_rve_constitutive_law(self, Element):
		
		# list of generated rve mdpa clones
		rve_mdpa_clones = []
		
		# get the number of integration points
		elemIntPoints = Element.GetIntegrationPoints()
		num_gp = len(elemIntPoints)
		elem_id = Element.Id
		
		# get a reference to the process into
		pinfo = self.MicroModelPart.Model.ProcessInfo
		
		# prepare the list of constitutive laws for the element
		constitutiveLaws = []
		
		# for each element integration point ...
		for gp_id in range(num_gp):
			
			# generate a new rve constitutive law
			rve_law__rve_mdpa__tuple = self.__generate_rve_constitutive_law()
			aRveLaw = rve_law__rve_mdpa__tuple[0]
			constitutiveLaws.append(aRveLaw)
			
			# TODO: check what rve law to track...
			# for the moment let's track them all
			self.__track_rve_constitutive_law(aRveLaw, elem_id, gp_id)
			
			# store the rve mdpa clone
			rve_mdpa_clones.append(rve_law__rve_mdpa__tuple[1])
		
		# assign the list of constitutive laws
		Element.SetValuesOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER, constitutiveLaws, pinfo)
		
		return rve_mdpa_clones
	
	## Initializes the output for the tracked rves (only if required)
	def __initialize_output(self):
		for key, value in self.TrackList.items():
			rveIO = value[1]
			if(rveIO is not None):
				rveIO.Initialize()
	
	## Writes the output for the tracked rves (only if required)
	def __write_output(self, currentTime):
		for key, value in self.TrackList.items():
			rveIO = value[1]
			if(rveIO is not None):
				rveIO.Write(currentTime)
	
	## Finalizes the output for the tracked rves (only if required)
	def __finalize_output(self):
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

