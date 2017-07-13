from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import read_materials_process
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


def initialize(CL_list, props_list, model_part, geom):
    print("Initializing CLs")
    for i, cl in enumerate(CL_list):
        cl.InitializeMaterial(model_part.Properties[props_list[i]], geom, km.Vector(3))
    print("")


def finalize_solution_step(CL_list, props_list, model_part, geom):
    print("Finalizing CLs")
    for i, cl in enumerate(CL_list):
        cl.FinalizeSolutionStep(model_part.Properties[props_list[i]], geom,
                km.Vector(3), model_part.ProcessInfo)
    print("")


def calculate_residual(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom):
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])
    An = np.zeros((nr_modes, nr_modes))
    bn = np.zeros(nr_modes)
    homog_stress = np.zeros(nr_comps)
    for i in range(nr_points):
        Bn = np.array(B_list[i])
        epsilon = epsilon_h + np.dot(Bn, x)
        sigma_n, C_n = _call_cl(epsilon, CL_list[i], geom, model_part.ProcessInfo, model_part.Properties[props_list[i]])
        w = iw_list[i]
        An += w * np.dot(Bn.transpose(), np.dot(C_n, Bn))
        bn += w * np.dot(Bn.transpose(), sigma_n)

        homog_stress += w * sigma_n

    return An, bn, homog_stress


def solve(x, epsilon_h, iw_list, CL_list, B_list, props_list, model_part, geom):
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])

    A, res, homog_stress = calculate_residual(x, epsilon_h, iw_list, CL_list,
            B_list, props_list, model_part, geom)
    it = 1
    norm_res = 1
    while(norm_res > 1e-9 and it < 10):
        Dx = -np.linalg.solve(A, res)
        x += Dx
        A, res, homog_stress = calculate_residual(x, epsilon_h, iw_list,
                CL_list, B_list, props_list, model_part, geom)
        norm_res = np.linalg.norm(res, ord=2)
        print("RESIDUAL CRITERION :: norm res: {:.3e}".format(norm_res))
        it += 1
    print("Convergence is achieved (or not)")

    return homog_stress 


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

    initialize(CL_list, props_list, model_part, geom)
    x = np.zeros(nr_modes)

    t0 = timer.time()
    nr_time_steps = 100
    end_time = 1.
    delta_time = end_time / nr_time_steps
    time = delta_time
    tolerance = delta_time / 10.
    homog_stress_list = []
    while(time <= end_time + tolerance):
        print("")
        print("Current time: {:.3f} --- Elapsed time: {:.2f}s".format(time,
            timer.time() - t0))
        partial_initial_strain = (time / end_time) * initial_strain 
        homog_stress_list.append(solve(x, partial_initial_strain, iw_list,
            CL_list, B_list, props_list, model_part, geom))
        finalize_solution_step(CL_list, props_list, model_part, geom)
        #print("OUTPUT MODES WEIGHT")
        #print(x)
        time = time + delta_time
    np.savetxt("homogenized_stress_hpromsolver.dat", homog_stress_list)
    print("Computing time = {:.2f}s".format(timer.time() - t0))
