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
        line = "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "\
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
            strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
            stress[0], stress[1], stress[2], stress[3], stress[4], stress[5])
        self.fo.write(line)


def generic_constitutive_law_test(model_part, deformation_test, load):

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
    properties.SetValue(km.ISOTROPIC_HARDENING_MODULUS, 0.3)

    # Construct a constitutive law
    cl = km.StructuralMechanicsApplication.SmallStrainIsotropicDamageTractionOnly3DLaw()

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
    constitutive_matrix = km.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

    # Setting the parameters
    cl_params = km.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetProcessInfo(model_part.ProcessInfo)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    cl.InitializeMaterial(properties, geom, N)

    # expected results
    deformation_test.initialize_reference_stress(cl.GetStrainSize())
    output = Output(deformation_test.output_filename)
    print(deformation_test.output_filename)

    print()
    print("Has INELASTIC_FLAG: ", cl.Has(km.StructuralMechanicsApplication.INELASTIC_FLAG))
    print("Has DAMAGE_VARIABLE: ", cl.Has(km.DAMAGE_VARIABLE))
    print("Has STRAIN_ENERGY: ", cl.Has(km.STRAIN_ENERGY))
    print("Has STRAIN: ", cl.Has(km.STRAIN))
    print("Has INITIAL_STRAIN ", model_part.ProcessInfo.Has(km.INITIAL_STRAIN))
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
    nr_ts = 21
    ts_array = [numpy.linspace(load[x+0], load[x+1], nr_ts) for x in range(len(load)-1)]
    ts_list = []
    for x in ts_array:
        ts_list.extend(x.tolist())
    for strain_mult in ts_list:
        deformation_test.set_deformation(cl_params, strain_mult)

        # Chauchy
        model_part.ProcessInfo[km.INITIAL_STRAIN] = cl_params.GetStrainVector()

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

        #output.printout(i, cl_params)
        output.write(cl_params)

        #reference_stress = deformation_test.get_reference_stress(i)
        stress = cl_params.GetStressVector()
        #print("Step ", i)
        #print("Reference: ", reference_stress)
        print("Stress:    ", stress)
        print("INELASTIC_FLAG: ", cl.GetValue(km.StructuralMechanicsApplication.INELASTIC_FLAG, bool()))

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        print("DAMAGE_VARIABLE: ", cl.CalculateValue(cl_params, km.DAMAGE_VARIABLE, float()))

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        print("STRAIN_ENERGY: ", cl.CalculateValue(cl_params, km.STRAIN_ENERGY, float()))

        zero_vector = km.Vector(6)
        zero_vector[0] = 0.
        zero_vector[1] = 0.
        zero_vector[2] = 0.
        zero_vector[3] = 0.
        zero_vector[4] = 0.
        zero_vector[5] = 0.
        cl_params.SetStrainVector(zero_vector)
        #print("STRAIN: ", cl.CalculateValue(cl_params, km.GREEN_LAGRANGE_STRAIN_VECTOR, km.Vector()))
        print()



class DeformationSmallStrainIsotropicDamage3D():
    def __init__(self, parameters):
        self.nr_timesteps = parameters["nr_timesteps"]
        self.load = parameters["load"]
        self.strain_list = parameters["strain"]
        self.output_filename = parameters["output_filename"]

    def initialize_reference_stress(self, strain_size):
        self.reference_stress = km.Vector(strain_size)
        for i in range(strain_size):
            self.reference_stress[i] = 0.0

    def set_deformation(self, cl_params, mult):
        self.strain = mult * self.initial_strain
        detF = 1
        F = km.Matrix(3,3)
        for i in range(3):
            for j in range(3):
                F[i, j] = 0
        for i in range(3):
            F[i, i] = 1
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(detF)
        cl_params.SetStrainVector(self.strain)

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = km.Vector(strain_size)
        self.initial_strain[0] = self.strain_list[0]
        self.initial_strain[1] = self.strain_list[1]
        self.initial_strain[2] = self.strain_list[2]
        self.initial_strain[3] = self.strain_list[3]
        self.initial_strain[4] = self.strain_list[4]
        self.initial_strain[5] = self.strain_list[5]

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(km.Vector(strain_size))
        r_stress[0][0] = 0.57692; r_stress[0][1] = 0.57692; r_stress[0][2] = 0.34615; r_stress[0][3] = 0.11538; r_stress[0][4] = 0.0; r_stress[0][5] = 0.11538
        r_stress[1][0] = 1.15384; r_stress[1][1] = 1.15384; r_stress[1][2] = 0.69231; r_stress[1][3] = 0.23077; r_stress[1][4] = 0.0; r_stress[1][5] = 0.23077
        r_stress[2][0] = 1.73076; r_stress[2][1] = 1.73076; r_stress[2][2] = 1.03850; r_stress[2][3] = 0.34615; r_stress[2][4] = 0.0; r_stress[2][5] = 0.34615
        r_stress[3][0] = 1.94550; r_stress[3][1] = 1.94550; r_stress[3][2] = 1.16730; r_stress[3][3] = 0.38910; r_stress[3][4] = 0.0; r_stress[3][5] = 0.38910
        r_stress[4][0] = 2.11858; r_stress[4][1] = 2.11858; r_stress[4][2] = 1.27120; r_stress[4][3] = 0.42372; r_stress[4][4] = 0.0; r_stress[4][5] = 0.42372
        r_stress[5][0] = 2.29166; r_stress[5][1] = 2.29166; r_stress[5][2] = 1.37500; r_stress[5][3] = 0.45833; r_stress[5][4] = 0.0; r_stress[5][5] = 0.45833
        r_stress[6][0] = 2.46473; r_stress[6][1] = 2.46473; r_stress[6][2] = 1.47880; r_stress[6][3] = 0.49295; r_stress[6][4] = 0.0; r_stress[6][5] = 0.49295
        r_stress[7][0] = 2.63781; r_stress[7][1] = 2.63781; r_stress[7][2] = 1.58270; r_stress[7][3] = 0.52756; r_stress[7][4] = 0.0; r_stress[7][5] = 0.52756
        r_stress[8][0] = 2.68543; r_stress[8][1] = 2.68543; r_stress[8][2] = 1.61130; r_stress[8][3] = 0.53709; r_stress[8][4] = 0.0; r_stress[8][5] = 0.53709
        r_stress[9][0] = 2.68543; r_stress[9][1] = 2.68543; r_stress[9][2] = 1.61130; r_stress[9][3] = 0.53709; r_stress[9][4] = 0.0; r_stress[9][5] = 0.53709
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]


#####################################################################
if __name__ == "__main__":

    model_part = km.Model().CreateModelPart("test")

    parameters = {"nr_timesteps": 20,
                  "load": [0, 1, 0, -1, 0],
                  "strain": [0.001, 0.000, 0.000, 0.000, 0.000, 0.000],
                  "output_filename": "strain-stress-damage-traction-only_traction.dat"
                  }
    load = [0, 1, 0, -1.2, 0, 1.4, 0, -1.6]
    deformation_test = DeformationSmallStrainIsotropicDamage3D(parameters)
    generic_constitutive_law_test(model_part, deformation_test, load)

    parameters = {"nr_timesteps": 20,
                  "load": [0, -1, 0, 1, 0],
                  "strain": [0.001, 0.000, 0.000, 0.000, 0.000, 0.000],
                  "output_filename": "strain-stress-damage-traction-only_compression.dat"
                  }
    load = [0, -1, 0, 1.2, 0, -1.4, 0, 1.6]
    deformation_test = DeformationSmallStrainIsotropicDamage3D(parameters)
    generic_constitutive_law_test(model_part, deformation_test, load)
