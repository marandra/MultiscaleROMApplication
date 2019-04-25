from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication
import os
import numpy


class Output:
    def __init__(self, filename):
        try:
            os.remove(filename)
        except OSError:
            pass
        self.fo = open(filename, "w")
        self.fo.write("#    {:<78}{:<30}\n".format("Strain (Voigt)", "Stress"))
        self.fo.write("#    {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n#\n"
                      "#    Column\n"
                      "#    {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n#\n"
                      "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "
                      "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
            "XX", "YY", "ZZ", "XY" , "YZ", "XZ", "XX", "YY", "ZZ", "XY" , "YZ", "XZ",
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.))

    def write(self, cl_params):
        strain = cl_params.GetStrainVector()
        stress = cl_params.GetStressVector()
        #line = "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "\
        #       "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
        #    strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
        #    stress[0], stress[1], stress[2], stress[3], stress[4], stress[5])
        line = "sr(00,0)={:<1.4e}; sr(00,1)={:<1.4e}; sr(00,2)={:<1.4e}; sr(00,3)={:<1.4e}; sr(00,4)={:<1.4e}; sr(00,5)={:<1.4e};\n".format(
            strain[0], strain[1], strain[2], strain[3], strain[4], strain[5])
        self.fo.write(line)


def constitutive_law_test(model_part, parameters):
    
    ###############
    # Setting up CL
    ###############

    # Define geometry
    dim = 3
    nnodes = 4
    N = km.Vector(nnodes)
    node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
    node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
    geom = km.Tetrahedra3D4(node1,node2,node3,node4)

    # Material properties
    properties = model_part.Properties[0]
    properties.SetValue(km.YOUNG_MODULUS, 3000)
    properties.SetValue(km.POISSON_RATIO, 0.3)
    properties.SetValue(km.YIELD_STRESS, 0.5)
    properties.SetValue(km.StructuralMechanicsApplication.INFINITY_YIELD_STRESS, 0.7)
    properties.SetValue(km.StructuralMechanicsApplication.HARDENING_MODULI_VECTOR, [0.3, 0.15])

    # Construct a constitutive law
    cl = km.StructuralMechanicsApplication.SmallStrainIsotropicDamageTractionOnly3DLaw()
    cl.Check(properties, geom, model_part.ProcessInfo)

    # Set the parameters to be employed
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    constitutive_matrix = km.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

    # Setting the parameters
    cl_params = km.ConstitutiveLawParameters()
    detF = 1
    F = km.Matrix(3,3)
    for i in range(3):
        for j in range(3):
            F[i, j] = 0
    for i in range(3):
        F[i, i] = 1
    cl_params.SetDeformationGradientF(F)
    cl_params.SetDeterminantF(detF)
    cl_params.SetOptions(cl_options)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetProcessInfo(model_part.ProcessInfo)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    #############
    # Testing
    #############

    cl.InitializeMaterial(properties, geom, N)

    print("Has STRAIN: ", cl.Has(km.STRAIN))
    print("Has INITIAL_STRAIN (=false)", model_part.ProcessInfo.Has(km.INITIAL_STRAIN))
    zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
    model_part.ProcessInfo[km.INITIAL_STRAIN] = zero_vector
    print("Has INITIAL_STRAIN (=true)", model_part.ProcessInfo.Has(km.INITIAL_STRAIN))
    print()

    strain_list = parameters["strain"]
    load_list = parameters["load"]

    strain_str = ""
    stress_str = ""
    cm_str = ""
    strain_energy_str = ""
    damage_variable_str = ""
    zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
    for t, load in enumerate(load_list):
        strain = km.Vector(6);
        for i, s in enumerate(strain_list):
            strain[i] = load * s
        print("Input: {} =  {} * {}".format(strain, load, strain_list))
        model_part.ProcessInfo[km.INITIAL_STRAIN] = strain
        #model_part.ProcessInfo[km.INITIAL_STRAIN] = zero_vector

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        #cl_params.SetStrainVector(strain)
        cl.InitializeMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        #cl_params.SetStrainVector(strain)
        cl.CalculateMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        #cl_params.SetStrainVector(strain)
        cl.FinalizeMaterialResponseCauchy(cl_params)

        v = cl_params.GetStrainVector()
        strain_str += "epr({0},0)={1:e}; epr({0},1)={2:e}; epr({0},2)={3:e}; epr({0},3)={4:e}; epr({0},4)={5:e}; epr({0},5)={6:e};\n".format(t, v[0], v[1], v[2], v[3], v[4], v[5])

        v = cl_params.GetStressVector()
        stress_str += "str({0},0)={1:e}; str({0},1)={2:e}; str({0},2)={3:e}; str({0},3)={4:e}; str({0},4)={5:e}; str({0},5)={6:e};\n".format(t, v[0], v[1], v[2], v[3], v[4], v[5])

        v = cl_params.GetConstitutiveMatrix()
        cm_str += "cmr({0}, 0)={1:e}; cmr({0}, 1)={2:e}; cmr({0}, 2)={3:e}; cmr({0}, 3)={4:e}; cmr({0}, 4)={5:e}; cmr({0}, 5)={6:e};\n".format(t, v[0,0],v[0,1],v[0,2],v[0,3],v[0,4],v[0,5])
        cm_str += "cmr({0}, 6)={1:e}; cmr({0}, 7)={2:e}; cmr({0}, 8)={3:e}; cmr({0}, 9)={4:e}; cmr({0},10)={5:e}; cmr({0},11)={6:e};\n".format(t, v[1,0],v[1,1],v[1,2],v[1,3],v[1,4],v[1,5])
        cm_str += "cmr({0},12)={1:e}; cmr({0},13)={2:e}; cmr({0},14)={3:e}; cmr({0},15)={4:e}; cmr({0},16)={5:e}; cmr({0},17)={6:e};\n".format(t, v[2,0],v[2,1],v[2,2],v[2,3],v[2,4],v[2,5])
        cm_str += "cmr({0},18)={1:e}; cmr({0},19)={2:e}; cmr({0},20)={3:e}; cmr({0},21)={4:e}; cmr({0},22)={5:e}; cmr({0},23)={6:e};\n".format(t, v[3,0],v[3,1],v[3,2],v[3,3],v[3,4],v[3,5])
        cm_str += "cmr({0},24)={1:e}; cmr({0},25)={2:e}; cmr({0},26)={3:e}; cmr({0},27)={4:e}; cmr({0},28)={5:e}; cmr({0},29)={6:e};\n".format(t, v[4,0],v[4,1],v[4,2],v[4,3],v[4,4],v[4,5])
        cm_str += "cmr({0},30)={1:e}; cmr({0},31)={2:e}; cmr({0},32)={3:e}; cmr({0},33)={4:e}; cmr({0},34)={5:e}; cmr({0},35)={6:e};\n".format(t, v[5,0],v[5,1],v[5,2],v[5,3],v[5,4],v[5,5])
        print(v)

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        damage_variable_str += "dvr[{}]={:e}; ".format(t, cl.CalculateValue(cl_params, km.DAMAGE_VARIABLE, float()))

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        strain_energy_str += "ser[{}]={:e}; ".format(t, cl.CalculateValue(cl_params, km.STRAIN_ENERGY, float()))

        zero_vector = km.Vector(6); zero_vector[0] = 0.; zero_vector[1] = 0.; zero_vector[2] = 0.; zero_vector[3] = 0.; zero_vector[4] = 0.; zero_vector[5] = 0.;
        cl_params.SetStrainVector(zero_vector)
        print()

    print(strain_str)
    print(stress_str)
    print(cm_str)
    print(strain_energy_str)
    print(damage_variable_str)

#####################################################################
if __name__ == "__main__":

    model_part = km.Model().CreateModelPart("test")

    print ("Testing traction-only CL: traction load -  unload")
    parameters = {"load": [0.1, 0.3, 1, -1],
                  "strain": [0.001, 0.001, 0.000, 0.001, 0.000, 0.000]
                  }
    constitutive_law_test(model_part, parameters)

    print ("Testing traction-only CL: compression load -  unload")
    parameters = {"load": [-0.1, -0.3, -1, 1],
                  "strain": [0.001, 0.001, 0.000, 0.001, 0.000, 0.000]
                  }
    constitutive_law_test(model_part, parameters)

