from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking

# Checking Benchmark Tests - Execution order:
# 0) In case of launch the tests for the first time, we need to get the reference data.
#    launch locally (at each benchmark test folder) the file "build_benchmark_reference.py"
#    this data must be obtained from the reference (or functional) program sources.
#    Otherwise, the test will run with success, even if the results are wrong...
# 1) Launch this script -> python3 MultiscaleROM_run_all_benchmarks.py
#    this routine automatically check if the results are the same comparing with the reference ones.

def Run():
    Msg = ""
    Text = "===== Multiscale ROM Application =====\n"

    # Multiscale ROM Test - Cube Damage Planestrain test
   
    Text += "Cube Damage PlaneStrain test: "
    os.chdir("cube-damage-planestrain-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Cube Damage Planestrain test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Cube Damage Planestrain test test FAILED")

    os.chdir("..")

    # Multiscale ROM Test - Cube Plasticity Planestrain test
   
    Text += "Cube Plasticity PlaneStrain test: "
    os.chdir("cube-plasticity-planestrain-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Cube Plasticity Planestrain test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Cube Plasticity Planestrain test test FAILED")

    os.chdir("..")

    # Multiscale ROM Test - Cube Plasticity Planestrain Shear test
   
    Text += "Cube Plasticity PlaneStrain Shear test: "
    os.chdir("cube-plasticity-shear-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Cube Plasticity Planestrain Shear test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Cube Plasticity Planestrain Shear test test FAILED")

    os.chdir("..")

    # Multiscale ROM Test - Cube Plasticity Planestrain Traction test
   
    Text += "Cube Plasticity PlaneStrain Traction test: "
    os.chdir("cube-plasticity-traction-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Cube Plasticity Planestrain Traction test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Cube Plasticity Planestrain Traction test test FAILED")

    os.chdir("..")

    # Multiscale ROM Test - Quad2D Damage Traction test
   
    Text += "Quad2D Damage Traction test: "
    os.chdir("quad-2D-damage-traction-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Quad 2D Damage Traction test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Quad 2D Damage Traction test test FAILED")

    os.chdir("..")
  
    # Multiscale ROM Test - Qua2D Plasticity Traction test
   
    Text += "Quad2D Plasticity Traction test: "
    os.chdir("quad-2D-plasticity-traction-test")
    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests---")

    print("running the benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("run_test.py", "benchmark_results.txt")

    if(successful==True):
        Text += "OK\n"
        print("Quad 2D Plasticity Traction test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Quad 2D Plasticity Traction test test FAILED")

    os.chdir("..")

    benchmarking.NotifyViaEmail("Status of Kratos examples", Text, ["mcaicedo@cimne.upc.edu"])

    return Text

    
if __name__ == '__main__':
    Run()

