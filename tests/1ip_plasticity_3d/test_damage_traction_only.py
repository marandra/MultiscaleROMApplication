from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication
import os

def _create_geometry(model_part, dim):
    # Create new nodes
    node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

    if (dim == 2):
        nnodes = 3

        # Allocate a geometry
        geom = km.Multiphysics.Triangle2D3(node1,node2,node3)
    elif (dim == 3):
        nnodes = 4
        node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

        # Allocate a geometry
        geom = km.Tetrahedra3D4(node1,node2,node3,node4)
    else:
        raise Exception("Error: bad dimension value: ", dim)
    return [geom, nnodes]


def _set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom):
    # Setting the parameters - note that a constitutive law may not need them all!
    cl_params = km.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetDeformationGradientF(F)
    cl_params.SetDeterminantF(detF)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetShapeFunctionsValues(N)
    cl_params.SetShapeFunctionsDerivatives(DN_DX)
    cl_params.SetProcessInfo(model_part.ProcessInfo)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    ## Do all sort of checks
    cl_params.CheckAllParameters() # Can not use this until the geometry is correctly exported to python
    cl_params.CheckMechanicalVariables()
    cl_params.CheckShapeFunctions()
    return cl_params


def _cl_check(cl, properties, geom, model_part, dim):
    cl.Check(properties, geom, model_part.ProcessInfo)

    if(cl.WorkingSpaceDimension() != dim):
        raise Exception("Mismatch between the WorkingSpaceDimension of the "
                        "Constitutive Law and the dimension of the space in "
                        "which the test is performed")


def _set_cl_options(dict_options):
    cl_options = km.Flags()
    if ("USE_ELEMENT_PROVIDED_STRAIN" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN,
                       dict_options["USE_ELEMENT_PROVIDED_STRAIN"])
    if ("COMPUTE_STRESS" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS,
                       dict_options["COMPUTE_STRESS"])
    if ("COMPUTE_CONSTITUTIVE_TENSOR" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR,
                       dict_options["COMPUTE_CONSTITUTIVE_TENSOR"])
    if ("COMPUTE_STRAIN_ENERGY" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRAIN_ENERGY,
                       dict_options["COMPUTE_STRAIN_ENERGY"])
    if ("ISOCHORIC_TENSOR_ONLY" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY,
                       dict_options["ISOCHORIC_TENSOR_ONLY"])
    if ("VOLUMETRIC_TENSOR_ONLY" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY,
                       dict_options["VOLUMETRIC_TENSOR_ONLY"])
    if ("FINALIZE_MATERIAL_RESPONSE" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE,
                       dict_options["FINALIZE_MATERIAL_RESPONSE"])

    # From here below it should be an output, not an input
    if ("FINITE_STRAINS" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.FINITE_STRAINS,
                       dict_options["FINITE_STRAINS"])
    if ("INFINITESIMAL_STRAINS" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.INFINITESIMAL_STRAINS,
                       dict_options["INFINITESIMAL_STRAINS"])
    if ("PLANE_STRAIN_LAW" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.PLANE_STRAIN_LAW,
                       dict_options["PLANE_STRAIN_LAW"])
    if ("PLANE_STRESS_LAW" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.PLANE_STRESS_LAW,
                       dict_options["PLANE_STRESS_LAW"])
    if ("AXISYMMETRIC_LAW" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.AXISYMMETRIC_LAW,
                       dict_options["AXISYMMETRIC_LAW"])
    if ("U_P_LAW" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.U_P_LAW,
                       dict_options["U_P_LAW"])
    if ("ISOTROPIC" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.ISOTROPIC,
                       dict_options["ISOTROPIC"])
    if ("ANISOTROPIC" in dict_options):
        cl_options.Set(km.ConstitutiveLaw.ANISOTROPIC,
                       dict_options["ANISOTROPIC"])
    return cl_options


class Output:
    def __init__(self, filename):
        try:
            os.remove(filename)
        except OSError:
            pass
        self.fo = open(filename, "w")
        self.fo.write("#    {:<78}{:<30}\n".format("Strain", "Stress"))
        self.fo.write("#    "
                      "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                      "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
            "XX", "YY", "ZZ", "XY" , "YZ", "XZ",
            "XX", "YY", "ZZ", "XY" , "YZ", "XZ"))

    def write(self, i, cl_params):
        line = "{:<4} "\
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "\
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
            i,
            cl_params.GetStrainVector()[0],
            cl_params.GetStrainVector()[1],
            cl_params.GetStrainVector()[2],
            cl_params.GetStrainVector()[3],
            cl_params.GetStrainVector()[4],
            cl_params.GetStrainVector()[5],

            cl_params.GetStressVector()[0],
            cl_params.GetStressVector()[1],
            cl_params.GetStressVector()[2],
            cl_params.GetStressVector()[3],
            cl_params.GetStressVector()[4],
            cl_params.GetStressVector()[5]
        )
        self.fo.write(line)

    def printout(self, i, cl_params):
        print("\nStep: {}".format(i))
        print("Strain = ", cl_params.GetStrainVector())
        print("Stress = ", cl_params.GetStressVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())


def generic_constitutive_law_test(model_part, deformation_test):
    # Define geometry
    [geom, nnodes] = _create_geometry(model_part, deformation_test.cl.dim)

    N = km.Vector(nnodes)
    DN_DX = km.Matrix(nnodes, deformation_test.cl.dim)

    # Material properties
    properties = deformation_test.cl.create_properties(model_part)

    # Construct a constitutive law
    cl = deformation_test.cl.create_constitutive_law()
    _cl_check(cl, properties, geom, model_part, deformation_test.cl.dim)

    # Set the parameters to be employed
    dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': True,
                    'COMPUTE_STRESS': True,
                    'COMPUTE_CONSTITUTIVE_TENSOR': True
                    }
    cl_options = _set_cl_options(dict_options)

    # Define deformation gradient
    F = deformation_test.get_init_deformation_gradientF()
    detF = 1.0

    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    constitutive_matrix = km.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

    # Setting the parameters - note that a constitutive law may not need them all!
    cl_params = _set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)

    cl.InitializeMaterial(properties, geom, N)

    # expected results
    deformation_test.initialize_reference_stress(cl.GetStrainSize())

    output = Output("strain-stress.dat")

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
    for i in range(deformation_test.nr_timesteps):
        deformation_test.set_deformation(cl_params, i)

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
        output.write(i, cl_params)

        reference_stress = deformation_test.get_reference_stress(i)
        stress = cl_params.GetStressVector()
        print("Step ", i)
        print("Reference: ", reference_stress)
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
        print("STRAIN: ", cl.CalculateValue(cl_params, km.GREEN_LAGRANGE_STRAIN_VECTOR, km.Vector()))
        print()


class LinearIsotropicDamageTractionOnly3D():
    def __init__(self):
        self.dim = 3
        self.young_modulus = 3000
        self.poisson_ratio = 0.3
        self.yield_stress = 2.0
        self.infinity_yield_stress = 3.0
        self.isotropic_hardening_modulus = 0.3

    def create_properties(self, model_part):
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(km.YOUNG_MODULUS, self.young_modulus)
        properties.SetValue(km.POISSON_RATIO, self.poisson_ratio)
        properties.SetValue(km.YIELD_STRESS, self.yield_stress)
        properties.SetValue(km.StructuralMechanicsApplication.INFINITY_YIELD_STRESS, self.infinity_yield_stress)
        properties.SetValue(km.ISOTROPIC_HARDENING_MODULUS, self.isotropic_hardening_modulus)
        return properties

    def create_constitutive_law(self):
        return km.StructuralMechanicsApplication.LinearIsotropicDamageTractionOnly3DLaw()


class Deformation():
    def __init__(self, parameters):
        self.nr_timesteps = parameters["nr_timesteps"]
        self.strain_list = parameters["strain"]

    def get_init_deformation_gradientF(self):
        self.F = km.Matrix(self.cl.dim,self.cl.dim)
        for i in range(self.cl.dim):
            for j in range(self.cl.dim):
                if(i==j):
                    self.F[i,j] = 1.0
                else:
                    self.F[i,j] = 0.0
        return self.F

    def initialize_reference_stress(self, strain_size):
        self.reference_stress = km.Vector(strain_size)
        for i in range(strain_size):
            self.reference_stress[i] = 0.0

    def set_deformation(self, cl_params, i):
        F = self.get_deformation_gradientF(i)
        detF = self.get_determinantF(i)
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(detF)


class DeformationLinearIsotropicDamageTractionOnly3D(Deformation):
    def __init__(self, parameters):
        Deformation.__init__(self, parameters)
        self.cl = LinearIsotropicDamageTractionOnly3D()

    def set_deformation(self, cl_params, i):
        self.strain = (i+1)/ self.nr_timesteps * self.initial_strain
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

    parameters = {"nr_timesteps": 10,
                 #[XX, YY, ZZ, XY, YZ, XZ]
                  "strain": [-0.001, -0.001, -0.001, -0.001, -0.001, -0.001]
                  }

    model_part = km.Model().CreateModelPart("test")

    # Test plasticity
    deformation_test = DeformationLinearIsotropicDamageTractionOnly3D(parameters)
    generic_constitutive_law_test(model_part, deformation_test)
