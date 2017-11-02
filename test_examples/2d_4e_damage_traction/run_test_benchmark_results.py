from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
#from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# benchmarking...
import sys
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

def WriteBenchmarkResults(model_part):

    print(benchmarking.InBenchmarkingMode())

    if (benchmarking.InBenchmarkingMode()):
        # Write Reaction Force in Node 4
        ref_Node = 4
        ref_Force = model_part.Nodes[ref_Node].GetSolutionStepValue(REACTION_X)
        #for node in model_part.Nodes:
        #    ir = node.GetSolutionStepValue(ROTATION_Y)
        #    if(ir > r_max):
        #        r_max = ir

        # write
        abs_tol = 1e-9
        rel_tol = 1e-5
        #print(r_max)
        benchmarking.Output(ref_Force, "REACTION_X", abs_tol, rel_tol)
