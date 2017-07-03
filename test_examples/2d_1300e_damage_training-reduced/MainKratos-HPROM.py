# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import pprint as pp
import time as timer
import operator
import KratosMultiphysics as km
import KratosMultiphysics.SolidMechanicsApplication as sol
import KratosMultiphysics.MultiscaleROMApplication as msr
import process_factory
import configparser
import sys
import numpy as np


def analysis(parameters, processes, solver, model_part):
    for process in processes:
        process.ExecuteInitialize()

    conf = configparser.ConfigParser()
    conf.read("reduced_bases.cfg")
    nr_modes = 10
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    roq_weights_filename = conf['Parameters']['roq_weights_filename']
    ngausspoints = int(conf['Parameters']['nr_integration_points'])
    voigtsize = int(conf['Parameters']['nr_strain_components'])

    # TODO this initialization should be done in scheme
    modes_weights = km.Vector(nr_modes)
    for i in range(nr_modes):
        modes_weights[i] = 0.0
    model_part.ProcessInfo[msr.REDUCED_MODES_WEIGHTS] = modes_weights
    model_part.ProcessInfo[msr.NUMBER_REDUCED_MODES] = nr_modes

    # TODO move this to a process
    strain_bases = np.load(strain_bases_filename)[:,:nr_modes]
    nr_dofs = ngausspoints * voigtsize
    index = 0
    for e, elem in enumerate(model_part.Elements):
            BE = km.Matrix(nr_dofs, nr_modes)
            for d in range(nr_dofs):
                for m in range(nr_modes):
                    BE[d, m] = strain_bases[index, m]
                index = index + 1
            elem.SetValue(msr.REDUCED_MODES_MATRIX, BE)
    #for elem in model_part.Elements:
    #        BE = km.Matrix(ngausspoints * voigtsize, nr_modes)
    #        for i in range(ngausspoints * voigtsize):
    #            line = fo.readline().strip().split()[:nr_modes]
    #            for j, value in enumerate(line):
    #                BE[i, j] = float(value)
    #        elem.SetValue(msr.REDUCED_MODES_MATRIX, BE)
    # TODO move this to a process
    with open(roq_weights_filename, "r") as fo:
        for elem in model_part.Elements:
            integration_weights = [float(x) for x in fo.readline().split()]
            elem.SetValue(msr.INTEGRATION_POINT_WEIGHT, integration_weights)
            
    for process in processes:
        process.ExecuteBeforeSolutionLoop()

    print("Finished reading reduced bases")

    with open("output.dat", "w") as fo:
        nr_time_steps = parameters["problem_data"]["nr_time_steps"].GetInt()
        end_time = parameters["problem_data"]["end_time"].GetDouble()
        delta_time = end_time / nr_time_steps
        time = delta_time
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
            # TODO there sould be a process to handle the output of weights
            #print("OUTPUT MODES WEIGHTS:")
            #print(model_part.ProcessInfo[msr.REDUCED_MODES_WEIGHTS])
            #for mode in model_part.ProcessInfo[msr.REDUCED_MODES_WEIGHTS]:
                #print(mode)
                #fo.write("{:17.15f} ".format(mode))
            #fo.write("\n")
            #print("\n")
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
    solver.ImportModelPart()
    solver.AddDofs()
    return solver, model_part

parameters = km.Parameters(open("ProjectParameters.json", 'r').read())
Model = create_model(parameters)
model_part = Model[parameters["problem_data"]["part_name"].GetString()]
solver, model_part = create_solver_complete_model_part(model_part, parameters)
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
