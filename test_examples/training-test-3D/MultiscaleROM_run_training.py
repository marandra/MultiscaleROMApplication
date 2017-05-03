from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys


def Run():
    Msg = ""
    Text = "===== Multiscale ROM Application =====\n"

    os.chdir("RVE_3D_v7_01")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 01---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")




    os.chdir("RVE_3D_v7_02")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 02---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")



    os.chdir("RVE_3D_v7_03")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 03---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")



    os.chdir("RVE_3D_v7_04")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 04---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")



    os.chdir("RVE_3D_v7_05")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 05---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")



    os.chdir("RVE_3D_v7_06")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 06---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")



    os.chdir("RVE_3D_v7_07")

    sys.path.append(os.getcwd())

    print("---start Multiscale ROM application tests 07---")

    os.system("runkratos MainKratos.py > OUTPUT.txt")

    os.chdir("..")


    return 0
    
if __name__ == '__main__':
    Run()

