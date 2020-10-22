import os
import numpy
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication


def write_header(fn_strain_stress, fn_const_matrix):
    try:
        os.remove(fn_strain_stress)
    except OSError:
        pass
    fs = open(fn_strain_stress, "w")
    fs.write("#    {:<78}{:<30}\n".format("Strain (Voigt)", "Stress"))
    fs.write("#    {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14}\n"
                  "#\n"
                  "#    Column\n"
                  "#    {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14} {:<14}\n"
                  "#\n".format(
        "XX", "YY", "ZZ", "XY" , "YZ", "XZ", "XX", "YY", "ZZ", "XY" , "YZ", "XZ",
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
        ))

    try:
        os.remove(fn_const_matrix)
    except OSError:
        pass
    fc = open(fn_const_matrix, "w")
    col = ""
    index = ""
    for i in range(6):
        for j in range(6):
            index += "{}, {}           ".format(i, j)
            col +=   "{:0>2}             ".format(i*6+j)
    fc.write("#" + index + "\n")
    fc.write("#" + col + "\n")

    return fs, fc


def write(fs, fc, cl_params):
    strain = cl_params.GetStrainVector()
    stress = cl_params.GetStressVector()
    line = "{:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}  "\
           "{:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}  {:<+1.6e}\n".format(
        strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
        stress[0], stress[1], stress[2], stress[3], stress[4], stress[5])
    fs.write(line)

    cm = cl_params.GetConstitutiveMatrix()
    line = ""
    for i in range(6):
        for j in range(6):
            line += "{:<+1.6e}  ".format(cm[i,j])
    fc.write(line + "\n")


def test(model_part, parameters, fs, fc):
    nr_timesteps = parameters["nr_timesteps"]
    load = parameters["load"]
    initial_strain = km.Vector(6)
    initial_strain[0] = parameters["strain"][0]
    initial_strain[1] = parameters["strain"][1]
    initial_strain[2] = parameters["strain"][2]
    initial_strain[3] = parameters["strain"][3]
    initial_strain[4] = parameters["strain"][4]
    initial_strain[5] = parameters["strain"][5]


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
    properties.SetValue(km.YOUNG_MODULUS, 21000)
    properties.SetValue(km.POISSON_RATIO, 0.3)
    properties.SetValue(km.YIELD_STRESS, 5.5)
    properties.SetValue(km.ISOTROPIC_HARDENING_MODULUS, 0.12924)
    properties.SetValue(km.StructuralMechanicsApplication.EXPONENTIAL_SATURATION_YIELD_STRESS, 0.0)
    properties.SetValue(km.HARDENING_EXPONENT, 0.1)

    # Construct a constitutive law
    cl = km.StructuralMechanicsApplication.SmallStrainJ2Plasticity3DLaw()

    cl.Check(properties, geom, model_part.ProcessInfo)
    if(cl.WorkingSpaceDimension() != 3):
        raise Exception("Mismatch between the WorkingSpaceDimension of the "
                        "Constitutive Law and the dimension of the space in "
                        "which the test is performed")


    # Set the parameters to be employed
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    constitutive_matrix = km.Matrix(cl.GetStrainSize(), cl.GetStrainSize())

    # Setting the parameters
    cl_params = km.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetProcessInfo(model_part.ProcessInfo)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)
    F = km.Matrix(3,3)
    F[0,0] = 1; F[0,1] = 0; F[0,2] = 0;
    F[1,0] = 0; F[1,1] = 1; F[1,2] = 0;
    F[2,0] = 0; F[2,1] = 0; F[2,2] = 1;
    cl_params.SetDeformationGradientF(F)
    cl_params.SetDeterminantF(1)

    cl.InitializeMaterial(properties, geom, N)

    print()
    print("Has ACCUMULATED_PLASTIC_STRAIN: ", cl.Has(km.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN))
    print("Has STRAIN_ENERGY: ", cl.Has(km.STRAIN_ENERGY))
    print("Has STRAIN: ", cl.Has(km.STRAIN))
    print("Has INTERNAL_VARIABLES: ", cl.Has(km.INTERNAL_VARIABLES))
    zero_vector = km.Vector(6)
    zero_vector[0] = 0.
    zero_vector[1] = 0.
    zero_vector[2] = 0.
    zero_vector[3] = 0.
    zero_vector[4] = 0.
    zero_vector[5] = 0.
    model_part.ProcessInfo[km.INITIAL_STRAIN] = zero_vector
    print("Has INITIAL_STRAIN ", model_part.ProcessInfo.Has(km.INITIAL_STRAIN))
    print()
    #nr_ts = 21
    ts_array = [numpy.linspace(load[x+0], load[x+1], nr_timesteps) for x in range(len(load)-1)]
    ts_list = []
    for x in ts_array:
        ts_list.extend(x.tolist())
    for strain_mult in ts_list:
        print("DEBUG {}".format(strain_mult))

        cl_params.SetStrainVector(strain_mult * initial_strain)

        # Chauchy
        #model_part.ProcessInfo[km.INITIAL_STRAIN] = cl_params.GetStrainVector()
        model_part.ProcessInfo[km.INITIAL_STRAIN] = strain_mult * initial_strain

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.InitializeMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.CalculateMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.FinalizeMaterialResponseCauchy(cl_params)

        write(fs, fc, cl_params)

        stress = cl_params.GetStressVector()

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        #print("ACCUMULATED_PLASTIC_STRAIN: ", cl.CalculateValue(cl_params, km.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN, float()))

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        #print("STRAIN_ENERGY: ", cl.CalculateValue(cl_params, km.STRAIN_ENERGY, float()))

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        #print("INTERNAL_VARIABLES: ", cl.CalculateValue(cl_params, km.INTERNAL_VARIABLES, km.Vector()))
        #print()


#####################################################################
if __name__ == "__main__":

    model_part = km.Model().CreateModelPart("test")

    fs, fc = write_header("plasticity_traction-strain-stress.dat", "plasticity_traction-const_matrix.dat")
    parameters = {"nr_timesteps": 6,
                  #"load": [0, 1, 0, -1, 0],
                  "load": [0, 1, 0, -1.2, 0, 1.4, 0, -1.6],
                  "strain": [0.001, 0.001, 0.000, 0.001, 0.000, 0.001],
                  }
    #deformation_test = DeformationSmallStrainJ2Plasticity3D(parameters)
    test(model_part, parameters, fs, fc)
    fs.close()
    fc.close()

    fs, fc = write_header("plasticity_compression-strain-stress.dat", "plasticity_compression-const_matrix.dat")
    parameters = {"nr_timesteps": 6,
                  #"load": [0, -1, 0, 1, 0],
                  "load": [0, -1, 0, 1.2, 0, -1.4, 0, 1.6],
                  "strain": [0.001, 0.001, 0.000, 0.001, 0.000, 0.001],
                  }
    #deformation_test = DeformationSmallStrainJ2Plasticity3D(parameters)
    test(model_part, parameters, fs, fc)
    fs.close()
    fc.close()
