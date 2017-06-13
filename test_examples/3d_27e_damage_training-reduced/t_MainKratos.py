# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import time as timer
import operator
import KratosMultiphysics as km
import KratosMultiphysics.SolidMechanicsApplication as sol
import KratosMultiphysics.MultiscaleROMApplication as msr
import process_factory
km.CheckForPreviousImport()


def analysis(parameters, processes, solver, model_part):
    for process in processes:
        process.ExecuteInitialize()
    for process in processes:
        process.ExecuteBeforeSolutionLoop()
    delta_time = parameters["problem_data"]["time_step"].GetDouble()
    time = parameters["problem_data"]["start_time"].GetDouble()
    end_time = parameters["problem_data"]["end_time"].GetDouble()
    tolerance = delta_time / 10.
    while(time <= end_time + tolerance):
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
        time = time + delta_time
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
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_1)
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_2)
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_3)
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_4)
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_5)
    model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_6)
    solver.ImportModelPart()
    solver.AddDofs()
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_1)
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_2)
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_3)
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_4)
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_5)
    model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_6)
    #constitutive_law_name = parameters["solver_settings"]["model_import_settings"]["constitutive_law"].GetString()
    #aux_obj_getter = operator.methodcaller(constitutive_law_name)
    #model_part.Properties[1].SetValue(km.CONSTITUTIVE_LAW, aux_obj_getter(msr))
    #model_part.Properties[1].SetValue(km.CONSTITUTIVE_LAW, aux_obj_getter(sol))
    #model_part.Properties[2].SetValue(km.CONSTITUTIVE_LAW, aux_obj_getter(sol))
    return solver, model_part

parameters = km.Parameters(open("ProjectParameters.json", 'r').read())
Model = create_model(parameters)
model_part = Model[parameters["problem_data"]["part_name"].GetString()]
solver, model_part = create_solver_complete_model_part(model_part, parameters)
#print("Adding Variables for MinimalKineticBC to model part", flush=True)
#model_part.AddNodalSolutionStepVariable(msr.LAGRANGIAN_DOF_1)
#model_part.AddNodalSolutionStepVariable(msr.LAGRANGIAN_DOF_2)
#model_part.AddNodalSolutionStepVariable(msr.LAGRANGIAN_DOF_3)
#model_part.AddNodalSolutionStepVariable(msr.SCALAR_LAGRANGE_MULTIPLIER_1)
#model_part.AddNodalSolutionStepVariable(km.SCALAR_LAGRANGE_MULTIPLIER)
#print("Adding DOFs for MinimalKineticBC to reference node", flush=True)
#model_part.Nodes[2].AddDof(msr.LAGRANGIAN_DOF_1)
#model_part.Nodes[2].AddDof(msr.LAGRANGIAN_DOF_2)
#model_part.Nodes[2].AddDof(msr.LAGRANGIAN_DOF_3)
#model_part.Nodes[2].AddDof(msr.SCALAR_LAGRANGE_MULTIPLIER_1)
#model_part.Nodes[2].AddDof(km.SCALAR_LAGRANGE_MULTIPLIER)
#print(model_part.Nodes, flush=True)
#build sub_model_parts or submeshes (rearrange parts for the application of custom processes)

# initialize GiD  I/O (gid outputs, file_lists)
output_settings = parameters["output_configuration"]
problem_name = parameters["problem_data"]["problem_name"].GetString()

solver.Initialize()


for i in range(parameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = parameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: model_part.GetSubModelPart(part_name)})

processes = process_factory.KratosProcessFactory(Model)\
    .ConstructListOfProcesses(parameters["constraints_process_list"])
processes += process_factory.KratosProcessFactory(Model)\
    .ConstructListOfProcesses(parameters["loads_process_list"])
#processes += process_factory.KratosProcessFactory(Model)\
#    .ConstructListOfProcesses(parameters["loads_rve_process_list"])

t0p = timer.clock()
t0w = timer.time()
analysis(parameters, processes, solver, model_part)
tfp = timer.clock()
tfw = timer.time()
print("Computing Time = {:.2f} s ({:.2f} s wall-time)".format(tfp - t0p, tfw - t0w))
print(timer.ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)
