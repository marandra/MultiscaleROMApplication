# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import read_materials_process
#import KratosMultiphysics.MultiscaleROMApplication as msr
import numpy as np
import json
import time as timer
import copy


def kratos_to_numpy_matrix(K, m, n):
    N = np.empty((m, n))
    for i in range(m):
        for j in range(n):
            N[i,j] = copy.copy(K[i,j])
    return N


def kratos_to_numpy_vector(K, n):
    print("Debug: {}".format(n))
    N = np.empty(n)
    for i in range(n):
        N[i] = copy.copy(K[i])
    return N


def _call_cl(epsilon, cl, geom, process_info, properties):

    N = km.Vector(3)
    #cl.InitializeMaterial(properties, geom, N)
    cl.Check(properties, geom, process_info)
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRAIN, False)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

    F = km.Matrix(3,3)
    F[0,0] = 1.0 + epsilon[0];   F[0,1] = 0.5 * epsilon[3]; F[0,2] = 0.5 * epsilon[5]
    F[1,0] = 0.5 * epsilon[3];   F[1,1] = 1.0 + epsilon[1]; F[1,2] = 0.5 * epsilon[4]
    F[2,0] = 0.5 * epsilon[5];   F[2,1] = 0.5 * epsilon[4]; F[2,2] = 1.0 + epsilon[2]
    FN = kratos_to_numpy_matrix(F, 3, 3)
    detF = np.linalg.det(FN)
    DN_DX = km.Matrix(3,2)
    constitutive_matrix = km.Matrix(cl.GetStrainSize(), cl.GetStrainSize())
    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    for i in range(cl.GetStrainSize()):
        stress_vector[i] = 0.
        strain_vector[i] = epsilon[i]

    #setting the parameters - note that a constitutive law may not need them all!
    cl_params = km.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetDeformationGradientF(F)
    cl_params.SetDeterminantF(detF)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetShapeFunctionsValues(N)
    cl_params.SetShapeFunctionsDerivatives(DN_DX)
    cl_params.SetProcessInfo(process_info)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    cl.CalculateMaterialResponseCauchy(cl_params)
    #print("The Material Response Cauchy")
    #print( "stress = ", cl_params.GetStressVector() )
    #print( "strain = ", cl_params.GetStrainVector() )
    #print( "C      = ", cl_params.GetConstitutiveMatrix() )
    #print("")

    cl.FinalizeSolutionStep(properties, geom, N, process_info)

    CM = cl_params.GetConstitutiveMatrix()
    stress = cl_params.GetStressVector()
    size = cl.GetStrainSize()
    CM_np = np.empty((size, size))
    stress_np = np.empty(size)
    for i in range(size):
        stress_np[i] = stress[i]
        for j in range(size):
            CM_np[i, j] = CM[i, j]

    return stress_np, CM_np
    #return  cl_params.GetStressVector(), cl_params.GetConstitutiveMatrix()


def initialize(iw_list, CL_list, B_list, props_list, model_part, geom):
    print("Initializing CLs")
    for i, cl in enumerate(CL_list):
        cl.InitializeMaterial(model_part.Properties[props_list[i]], geom, km.Vector(3))
    print("")

def calculate_residual(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom):
    nr_points = len(B_list)
    #nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])
    An = np.zeros((nr_modes, nr_modes))
    bn = np.zeros(nr_modes)
    for i in range(nr_points):
        Bn = np.array(B_list[i])
        epsilon = epsilon_h + np.dot(Bn, x)
        sigma_n, C_n = _call_cl(epsilon, CL_list[i], geom, model_part.ProcessInfo, model_part.Properties[props_list[i]])
        #Cn = kratos_to_numpy_matrix(C, nr_comps, nr_comps)
        #print(C_n)
        #print(np.dot(Bn.transpose(), np.dot(Cn, Bn)))
        #print("")
        #sigma_n = kratos_to_numpy_vector(sigma, nr_comps)
        w = iw_list[i]
        An += w * np.dot(Bn.transpose(), np.dot(C_n, Bn))
        bn += w * np.dot(Bn.transpose(), sigma_n)
    return An, bn


def solve(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom):
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])

    it = 0
    tol = 1e-9
    A, res = calculate_residual(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom)
    norm_res = 1
    #print("iteration = {}".format(it))
    #print(A)
    #print(x)
    #print(res)
    print("norm res: {:.3e}".format(norm_res))
    it += 1
    while(norm_res > tol and it < 10):
        #print("iteration = {}".format(it))
        x -= np.linalg.solve(A, res)
        A, res = calculate_residual(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom)
        norm_res = np.linalg.norm(res, ord=2)
        #print(A)
        #print(x)
        #print(res)
        print("norm res: {:.3e}".format(norm_res))
        it += 1
    print("Convergence is achieved (or not)")

##############################################################
if __name__ == "__main__":

    # create model part and assign properties
    model_part = km.ModelPart("Main")
    Model = {"Main" : model_part}
    settings = km.Parameters("""
                {
                    "Parameters": {
                            "materials_filename": "materials.json"
                    }
            }
            """)
    read_materials_process.Factory(settings, Model)
    #print(model_part.Properties[1])
    #print(model_part.Properties[2])

    # get integration point weights
    with open("rve.json", "r") as fi:
        out = json.load(fi)
    iw_list = out['w']
    props_list = out['props_id']
    B_list = out['B']
    
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])
    print("Nr points: {}".format(nr_points))
    print("Nr components: {}".format(nr_comps))
    print("Nr modes: {}".format(nr_modes))
    print("")

    BK_list = []
    CL_list = []
    for i in range(nr_points):
        BK = km.Matrix(nr_comps, nr_modes)
        for c in range(nr_comps):
            for m in range(nr_modes):
                BK[c, m] = B_list[i][c][m]
        BK_list.append(BK)
        CL_list.append(model_part.Properties[props_list[i]][km.CONSTITUTIVE_LAW].Clone())

    initial_strain = km.Vector(nr_comps)
    initial_strain[0] = 0.001
    initial_strain[1] = 0.
    initial_strain[2] = 0.
    initial_strain[3] = 0.
    initial_strain[4] = 0.
    initial_strain[5] = 0.

    #setting of constant values
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])
    node1 = model_part.CreateNewNode(1,0.0,0.0,0.0)
    geom = km.Triangle2D3(node1, node1, node1)

    initialize(iw_list, CL_list, B_list, props_list, model_part, geom)
    x = np.zeros(nr_modes)

    t0 = timer.time()
    nr_time_steps = 100
    end_time = 1.
    delta_time = end_time / nr_time_steps
    time = delta_time
    tolerance = delta_time / 10.
    while(time <= end_time + tolerance):
        print("")
        print("Current time: {:.3f} --- Elapsed time: {:.2f}s".format(time, timer.time() - t0))
        partial_initial_strain = (time / end_time) * initial_strain 
        solve(x, partial_initial_strain, iw_list, CL_list, B_list, props_list, model_part, geom)
        print(x)
        time = time + delta_time

    print("Computing time = {:.2f}s".format(timer.time() - t0))
