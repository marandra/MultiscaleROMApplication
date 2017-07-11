# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import read_materials_process
import KratosMultiphysics.MultiscaleROMApplication as msr
import numpy as np
import json
import sys


def kratos_to_numpy_matrix(K, m, n):
    N = np.empty((m, n))
    for i in range(m):
        for j in range(n):
            N[i,j] = K[i,j]
    return N


def kratos_to_numpy_vector(K, n):
    N = np.empty(n)
    for i in range(n):
        N[i] = K[i]
    return N


def _call_cl(epsilon_h, cl, geom, process_info, properties):
    cl.Check(properties, geom, process_info)
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRAIN, False)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

    F = km.Matrix(3,3)
    F[0,0] = 1.0 + epsilon_h[0];   F[0,1] = 0.5 * epsilon_h[3]; F[0,2] = 0.5 * epsilon_h[5]
    F[1,0] = 0.5 * epsilon_h[3];   F[1,1] = 1.0 + epsilon_h[1]; F[1,2] = 0.5 * epsilon_h[4]
    F[2,0] = 0.5 * epsilon_h[5];   F[2,1] = 0.5 * epsilon_h[4]; F[2,2] = 1.0 + epsilon_h[2]
    FN = kratos_to_numpy_matrix(F, 3, 3)
    detF = np.linalg.det(FN)
    N = km.Vector(3)
    DN_DX = km.Matrix(3,2)
    constitutive_matrix = km.Matrix(cl.GetStrainSize(), cl.GetStrainSize())
    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    for i in range(cl.GetStrainSize()):
        stress_vector[i] = 0.
        strain_vector[i] = epsilon_h[i]

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
    #print( "_train = ", cl_params.GetStrainVector() )
    #print( "C      = ", cl_params.GetConstitutiveMatrix() )
    #print("")

    return  cl_params.GetStressVector(), cl_params.GetConstitutiveMatrix()


def calculate_alpha(epsilon_h, iw_list, CL_list, B_list, props_list, model_part):

    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])

    node1 = model_part.CreateNewNode(1,0.0,0.0,0.0)
    geom = km.Triangle2D3(node1, node1, node1)

    N = km.Vector(3)

    for i, cl in enumerate(CL_list):
        cl.InitializeMaterial(model_part.Properties[props_list[i]], geom, N)

    while(True):
        An = np.zeros((nr_modes, nr_modes))
        bn = np.zeros(nr_modes)
        for i in range(nr_points):
            sigma, C = _call_cl(epsilon_h, CL_list[i], geom, model_part.ProcessInfo, model_part.Properties[props_list[i]])
            Bn = np.array(B_list[i])
            Cn = kratos_to_numpy_matrix(C, nr_comps, nr_comps)
            sigma_n = kratos_to_numpy_vector(sigma, nr_comps)
            w = iw_list[i]
            An += w * np.dot(Bn.transpose(), np.dot(Cn, Bn))
            bn += w * np.dot(Bn.transpose(), sigma_n)

        x = np.linalg.solve(An, bn)

        print(An)
        print(bn)
        print(x)



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
    print(model_part.Properties[1])
    print(model_part.Properties[2])
    #model_part.AddProperties(km.Properties(1))
    #properties = model_part.Properties[1]
    #properties.SetValue(km.YOUNG_MODULUS, 21000)
    #properties.SetValue(km.POISSON_RATIO, 0.3)
    #properties.SetValue(km.YIELD_STRESS, 5.5)
    #properties.SetValue(msr.INFINITY_YIELD_STRESS, 100000)
    #properties.SetValue(msr.ISOTROPIC_DAMAGE_MODULUS, 0.01)
    #properties.SetValue(msr.FLOW_RULE_IS_TRACTION_ONLY, False)

    # get integration point weights
    with open("rve.json", "r") as fi:
        out = json.load(fi)
    iw_list = out['w']
    props_list = out['props_id']
    B_list = out['B']
    
    nr_points = len(B_list)
    nr_comps = len(B_list[0])
    nr_modes = len(B_list[0][0])

    print(nr_points)
    print(nr_comps)
    print(nr_modes)

    BK_list = []
    CL_list = []
    for i in range(nr_points):
        BK = km.Matrix(nr_comps, nr_modes)
        for c in range(nr_comps):
            for m in range(nr_modes):
                BK[c, m] = B_list[i][c][m]
        BK_list.append(BK)
        CL_list.append(model_part.Properties[props_list[i]][km.CONSTITUTIVE_LAW].Clone())

    nr_timesteps
    epsilon_h = km.Vector(nr_comps)
    for t in range(nr_timesteps):
        epsilon_h[0] = 0.00001 * t / nr_timesteps
        epsilon_h[1] = 0.00000 * t / nr_timesteps
        epsilon_h[2] = 0.00000 * t / nr_timesteps
        epsilon_h[3] = 0.00000 * t / nr_timesteps
        epsilon_h[4] = 0.00000 * t / nr_timesteps
        epsilon_h[5] = 0.00000 * t / nr_timesteps

        calculate_alpha(epsilon_h, iw_list, CL_list, B_list, props_list, model_part)
