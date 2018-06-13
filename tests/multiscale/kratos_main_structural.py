# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import time as timer
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication
km.CheckForPreviousImport()


def ApplyLoad(initial_value, model_part, time):
    for node in model_part.Nodes:
        factor = 1.
        value = initial_value * factor * time
        node.SetSolutionStepValue(km.POINT_LOAD_Y, 0, value)


## Import define_output
with open("ProjectParameters.json",'r') as parameter_file:
    ProjectParameters = km.Parameters(parameter_file.read())

echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

## Import parallel modules if needed
if (parallel_type == "MPI"):
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.TrilinosApplication import *

## Structure model part definition
main_model_part_name = ProjectParameters["problem_data"]["model_part_name"].GetString()
main_model_part = km.ModelPart(main_model_part_name)
main_model_part.ProcessInfo.SetValue(km.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

## Solver construction
import python_solvers_wrapper_structural
solver = python_solvers_wrapper_structural.CreateSolver(main_model_part, ProjectParameters)
solver.AddVariables()
#main_model_part.AddNodalSolutionStepVariable(km.LAGRANGE_DISPLACEMENT)
solver.ImportModelPart()
solver.AddDofs()

## Initialize GiD I/O
output_post = ProjectParameters.Has("output_configuration")
if (output_post == True):
    if (parallel_type == "OpenMP"):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                      ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                      ProjectParameters["output_configuration"])
    elif (parallel_type == "MPI"):
        from gid_output_process_mpi import GiDOutputProcessMPI
        gid_output = GiDOutputProcessMPI(solver.GetComputingModelPart(),
                                         ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                         ProjectParameters["output_configuration"])

    gid_output.ExecuteInitialize()

## Creation of the Kratos model (build sub_model_parts or submeshes)
StructureModel = km.Model()
StructureModel.AddModelPart(main_model_part)

## Print model_part and properties
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    Logger.PrintInfo("ModelPart", main_model_part)
    for properties in main_model_part.Properties:
        Logger.PrintInfo("Property " + str(properties.Id), properties)

## Processes construction
import process_factory
list_of_processes = process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["loads_process_list"])
if (ProjectParameters.Has("list_other_processes") == True):
    list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
if (ProjectParameters.Has("json_output_process") == True):
    list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["json_output_process"])

if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    count = 0
    for process in list_of_processes:
        count += 1
        Logger.PrintInfo("Process " + str(count), process)

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

## Solver initialization
solver.Initialize()
solver.SetEchoLevel(echo_level)

if (output_post == True):
    gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Writing the full ProjectParameters file before solving
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 0):
    f = open("ProjectParametersOutput.json", 'w')
    f.write(ProjectParameters.PrettyPrintJsonString())
    f.close()

## Stepping and time settings
t0p = timer.clock()
t0w = timer.time()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

if main_model_part.ProcessInfo[km.IS_RESTARTED] == True:
    time = main_model_part.ProcessInfo[km.TIME]
else:
    time = start_time
    main_model_part.ProcessInfo[km.STEP] = 0

if (parallel_type == "OpenMP") or (mpi.rank == 0):
    km.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -START- ")

# Solving the problem (time integration)
while(time <= end_time):

    time = time + delta_time
    main_model_part.ProcessInfo[km.STEP] += 1
    main_model_part.CloneTimeStep(time)

    if (parallel_type == "OpenMP") or (mpi.rank == 0):
        km.Logger.PrintInfo("")
        km.Logger.PrintInfo("STEP: ", main_model_part.ProcessInfo[km.STEP])
        km.Logger.PrintInfo("TIME: ", time)

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    if (output_post == True):
        gid_output.ExecuteInitializeSolutionStep()

    solver.Solve()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    if (output_post == True):
        gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    if (output_post == True) and (gid_output.IsOutputStep()):
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    solver.SaveRestart()

for process in list_of_processes:
    process.ExecuteFinalize()

if (output_post == True):
    gid_output.ExecuteFinalize()

if (parallel_type == "OpenMP") or (mpi.rank == 0):
    tfp = timer.clock()
    tfw = timer.time()
    km.Logger.PrintInfo("")
    km.Logger.PrintInfo("::[KSM Simulation]:: ", "Computing Time = {:.2f} s ({:.2f} s wall-time)".format(tfp - t0p, tfw - t0w))
    km.Logger.PrintInfo("::[KSM Simulation]:: ", "Analysis -END- ")
