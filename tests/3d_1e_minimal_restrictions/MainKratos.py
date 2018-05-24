# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import time as timer
#import operator
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication as msr
import KratosMultiphysics.ExternalSolversApplication
km.CheckForPreviousImport()


def ApplyLoad(initial_value, model_part, time):
    for node in model_part.Nodes:
        factor = 1.
        value = initial_value * factor * time
        node.SetSolutionStepValue(km.POINT_LOAD_Y, 0, value)


########################
# main
########################

ProjectParameters = km.Parameters(open("ProjectParameters.json", 'r').read())
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

## Modelpart definition
main_model_part = km.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(km.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

## Solver construction
import python_solvers_wrapper_structural
solver = python_solvers_wrapper_structural.CreateSolver(main_model_part, ProjectParameters)
solver.AddVariables()
main_model_part.AddNodalSolutionStepVariable(km.LAGRANGE_DISPLACEMENT)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_1)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_2)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_3)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_4)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_5)
main_model_part.AddNodalSolutionStepVariable(msr.LAGRANGE_MULTIPLIER_6)
solver.ImportModelPart()
solver.AddDofs()
lmutility = msr.LagrangeMultiplierUtility(main_model_part)
lmutility.Execute()
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_1)
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_2)
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_3)
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_4)
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_5)
main_model_part.Nodes[1].AddDof(msr.LAGRANGE_MULTIPLIER_6)

## Initialize GiD I/O
output_post = ProjectParameters.Has("output_configuration")
if output_post:
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                  ProjectParameters["problem_data"]["problem_name"].GetString(),
                                  ProjectParameters["output_configuration"])
    gid_output.ExecuteInitialize()

## Creation of the Kratos model (build sub_model_parts or submeshes)
StructureModel = {ProjectParameters["problem_data"]["model_part_name"].GetString(): main_model_part}

## Get the list of the sub_model_parts in where the processes are to be applied
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    StructureModel.update({part_name: main_model_part.GetSubModelPart(part_name)})

## Print model_part and properties
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

## Processes construction
import process_factory
list_of_processes = process_factory.KratosProcessFactory(StructureModel)\
    .ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
list_of_processes += process_factory.KratosProcessFactory(StructureModel)\
    .ConstructListOfProcesses(ProjectParameters["loads_process_list"])
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    for process in list_of_processes:
        print(process)

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

t0p = timer.clock()
t0w = timer.time()
## Temporal loop
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
time = start_time
main_model_part.ProcessInfo[km.TIME_STEPS] = 0

while(time <= end_time - delta_time):
    time += delta_time
    main_model_part.ProcessInfo[km.TIME_STEPS] += 1
    main_model_part.CloneTimeStep(time)
    if (parallel_type == "OpenMP") or (mpi.rank == 0):
        print("")
        print("STEP = ", main_model_part.ProcessInfo[km.TIME_STEPS])
        print("TIME = ", time)

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
    if output_post:
        gid_output.ExecuteInitializeSolutionStep()

    #ApplyLoad(-5250000, main_model_part.GetSubModelPart("LOAD"), time)
    solver.Solve()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    if output_post:
        gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    if output_post and gid_output.IsOutputStep():
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


for process in list_of_processes:
    process.ExecuteFinalize()
if output_post:
    gid_output.ExecuteFinalize()

tfp = timer.clock()
tfw = timer.time()
print("")
print("Computing Time = {:.2f} s ({:.2f} s wall-time)".format(tfp - t0p, tfw - t0w))
print(timer.ctime())

