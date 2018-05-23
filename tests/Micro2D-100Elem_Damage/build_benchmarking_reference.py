from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
#kratos_benchmarking_path = '../../../../benchmarking'
kratos_benchmarking_path = '/home/manuel/Project/Kratos/benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

print("Building reference data for Micro2D 100Elem Damage CL Test")
benchmarking.BuildReferenceData("run_test.py", "benchmark_results.txt")
