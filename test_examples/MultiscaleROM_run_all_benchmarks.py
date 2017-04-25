from __future__ import print_function, absolute_import, division
import os
import sys
sys.path.append('../../../benchmarking')
import benchmarking

def run_case(case_name):
    os.chdir(case_name)
    sys.path.append(os.getcwd())
    successful, Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")
    if successful:
        result = case_name + " OK"
    else:
        result = case_name + " FAILED"
    os.chdir("..")
    return result


case_list = [
    "cube-damage-planestrain-test",
    "cube-plasticity-planestrain-test",
    "cube-plasticity-shear-test",
    "cube-plasticity-traction-test",
    "quad-2D-damage-traction-test",
    "quad-2D-plasticity-traction-test",
]
    
for case in case_list:
    print(run_case(case))

