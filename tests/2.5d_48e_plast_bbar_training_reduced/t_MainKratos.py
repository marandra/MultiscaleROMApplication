# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import time as timer
import operator
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication as msr
import KratosMultiphysics.ExternalSolversApplication
import process_factory
from gid_output_process import GiDOutputProcess

# For Benchmarking purposes
#from run_test_benchmark_results import *

km.CheckForPreviousImport()

def analysis(parameters, processes, gid_output, solver, model_part):
    for process in processes:
        process.ExecuteInitialize()
    gid_output.ExecuteInitialize()
    for process in processes:
        process.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteBeforeSolutionLoop()
    nr_time_steps = parameters["problem_data"]["nr_time_steps"].GetInt()
    end_time = parameters["problem_data"]["end_time"].GetDouble()
    delta_time = end_time / nr_time_steps
    time = delta_time
    tolerance = delta_time / 10.
    while(time <= end_time + tolerance):
        model_part.CloneTimeStep(time)
        for process in processes:
            process.ExecuteInitializeSolutionStep()
        gid_output.ExecuteInitializeSolutionStep()

        solver.Solve()

        for process in processes:
            process.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalizeSolutionStep()
        for process in processes:
            process.ExecuteBeforeOutputStep()
        if gid_output.IsOutputStep():
            gid_output.PrintOutput()
        for process in processes:
            process.ExecuteAfterOutputStep()
        time = time + delta_time
    for process in processes:
        process.ExecuteFinalize()
    gid_output.ExecuteFinalize()

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
    model_part.AddNodalSolutionStepVariable(km.LAGRANGE_DISPLACEMENT)
    solver.ImportModelPart()
    solver.AddDofs()
    return solver, model_part

parameters = km.Parameters(open("ProjectParameters.json", 'r').read())
Model = create_model(parameters)
model_part = Model[parameters["problem_data"]["part_name"].GetString()]
solver, model_part = create_solver_complete_model_part(model_part, parameters)
output_settings = parameters["output_configuration"]
problem_name = parameters["problem_data"]["problem_name"].GetString()
gid_output = GiDOutputProcess(model_part, problem_name, output_settings)

solver.Initialize()


for i in range(parameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = parameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: model_part.GetSubModelPart(part_name)})

processes = process_factory.KratosProcessFactory(Model)\
    .ConstructListOfProcesses(parameters["constraints_process_list"])
processes += process_factory.KratosProcessFactory(Model)\
    .ConstructListOfProcesses(parameters["loads_process_list"])

t0p = timer.clock()
t0w = timer.time()
analysis(parameters, processes, gid_output, solver, model_part)
tfp = timer.clock()
tfw = timer.time()
print("Computing Time = {:.2f} s ({:.2f} s wall-time)".format(tfp - t0p, tfw - t0w))
print(timer.ctime())

