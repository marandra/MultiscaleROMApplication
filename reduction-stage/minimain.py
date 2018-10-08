from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication as kmsr
import configparser

import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

def _create_geometry(model_part, dim):
    # Create new nodes
    node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

    if (dim == 2):
        nnodes = 3

        # Allocate a geometry
        geom = KratosMultiphysics.Triangle2D3(node1,node2,node3)
    elif (dim == 3):
        nnodes = 4
        node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

        # Allocate a geometry
        geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3,node4)
    else:
        raise Exception("Error: bad dimension value: ", dim)
    return [geom, nnodes]


def _set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom):
    # Setting the parameters - note that a constitutive law may not need them all!
    cl_params = KratosMultiphysics.ConstitutiveLawParameters()
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
        raise Exception("Mismatch between the WorkingSpaceDimension of the Constitutive Law and the dimension of the space in which the test is performed")

def _set_cl_options(dict_options):
    cl_options = KratosMultiphysics.Flags()
    if ("USE_ELEMENT_PROVIDED_STRAIN" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, dict_options["USE_ELEMENT_PROVIDED_STRAIN"])
    if ("COMPUTE_STRESS" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, dict_options["COMPUTE_STRESS"])
    if ("COMPUTE_CONSTITUTIVE_TENSOR" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, dict_options["COMPUTE_CONSTITUTIVE_TENSOR"])
    if ("COMPUTE_STRAIN_ENERGY" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, dict_options["COMPUTE_STRAIN_ENERGY"])
    if ("ISOCHORIC_TENSOR_ONLY" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, dict_options["ISOCHORIC_TENSOR_ONLY"])
    if ("VOLUMETRIC_TENSOR_ONLY" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, dict_options["VOLUMETRIC_TENSOR_ONLY"])
    if ("FINALIZE_MATERIAL_RESPONSE" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, dict_options["FINALIZE_MATERIAL_RESPONSE"])

    # From here below it should be an otput not an input
    if ("FINITE_STRAINS" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINITE_STRAINS, dict_options["FINITE_STRAINS"])
    if ("INFINITESIMAL_STRAINS" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.INFINITESIMAL_STRAINS, dict_options["INFINITESIMAL_STRAINS"])
    if ("PLANE_STRAIN_LAW" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRAIN_LAW, dict_options["PLANE_STRAIN_LAW"])
    if ("PLANE_STRESS_LAW" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRESS_LAW, dict_options["PLANE_STRESS_LAW"])
    if ("AXISYMMETRIC_LAW" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.AXISYMMETRIC_LAW, dict_options["AXISYMMETRIC_LAW"])
    if ("U_P_LAW" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.U_P_LAW, dict_options["U_P_LAW"])
    if ("ISOTROPIC" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOTROPIC, dict_options["ISOTROPIC"])
    if ("ANISOTROPIC" in dict_options):
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ANISOTROPIC, dict_options["ANISOTROPIC"])
    return cl_options


def _print_cl_output(cl, cl_params, properties, geom, N, model_part):
    print("The Material Response PK2")
    cl.CalculateMaterialResponsePK2(cl_params)
    print("Stress = ", cl_params.GetStressVector())
    print("Strain = ", cl_params.GetStrainVector())
    print("C      = ", cl_params.GetConstitutiveMatrix())

    cl.FinalizeMaterialResponsePK2(cl_params)
    cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

    print("\nThe Material Response Kirchhoff")
    cl.CalculateMaterialResponseKirchhoff(cl_params)
    print("Stress = ", cl_params.GetStressVector())
    print("Strain = ", cl_params.GetStrainVector())
    print("C      = ", cl_params.GetConstitutiveMatrix())

    cl.FinalizeMaterialResponseKirchhoff(cl_params)
    cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

    print("\nThe Material Response Cauchy")
    cl.CalculateMaterialResponseCauchy(cl_params)
    print("Stress = ", cl_params.GetStressVector())
    print("Strain = ", cl_params.GetStrainVector())
    print("C      = ", cl_params.GetConstitutiveMatrix())

    cl.FinalizeMaterialResponseCauchy(cl_params)
    cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

def _generic_constitutive_law_test(model_part, deformation_test):
    # Define geometry
    [geom, nnodes] = _create_geometry(model_part, deformation_test.cl.dim)

    N = KratosMultiphysics.Vector(nnodes)
    DN_DX = KratosMultiphysics.Matrix(nnodes, deformation_test.cl.dim)

    # Material properties
    properties = deformation_test.cl.create_properties(model_part)

    # Construct a constitutive law
    cl = deformation_test.cl.create_constitutive_Law()
    _cl_check(cl, properties, geom, model_part, deformation_test.cl.dim)

    # Set the parameters to be employed
    dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                    'COMPUTE_STRESS': True,
                    'COMPUTE_CONSTITUTIVE_TENSOR': True,
                    'FINITE_STRAINS': True,
                    'ISOTROPIC': True,
                    }
    cl_options = _set_cl_options(dict_options)

    # Define deformation gradient
    F = deformation_test.get_init_deformation_gradientF()
    detF = 1.0

    stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
    strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
    constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

    # Setting the parameters - note that a constitutive law may not need them all!
    cl_params = _set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)
    cl.InitializeMaterial(properties, geom, N)

    # Check the results
    deformation_test.initialize_reference_stress(cl.GetStrainSize())

    for i in range(deformation_test.nr_timesteps):
        deformation_test.set_deformation(cl_params, i)

        # Chauchy
        cl.CalculateMaterialResponseCauchy(cl_params)
        cl.FinalizeMaterialResponseCauchy(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
        reference_stress = deformation_test.get_reference_stress(i)

        stress = cl_params.GetStressVector()
        #DEBUG
        print(stress)

        #tolerance = 1.0e-4
        #for j in range(cl.GetStrainSize()):
        #    if (abs(stress[j]) > tolerance):
        #        assertAlmostEqual((reference_stress[j] - stress[j])/stress[j], 0.0, msg=("Error checking solution " + str(stress[j]) + " different from " + str(reference_stress[j]) + " with tolerance of " + str(tolerance)), delta=tolerance)


class LinearElastic():
   def __init__(self):
       self.young_modulus = 200e9
       self.poisson_ratio = 0.3

   def create_properties(self, model_part):
       prop_id = 0
       properties = model_part.Properties[prop_id]
       properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, self.young_modulus)
       properties.SetValue(KratosMultiphysics.POISSON_RATIO, self.poisson_ratio)
       return properties

class LinearJ2Plasticity(LinearElastic):
    def __init__(self):
        self.young_modulus = 21000
        self.poisson_ratio = 0.3
        self.yield_stress = 5.5
        self.reference_hardening_modulus = 1.0
        self.isotropic_hardening_modulus = 0.12924
        self.infinity_hardening_modulus = 0.0
        self.hardening_exponent = 1.0

    def create_properties(self, model_part):
        properties = LinearElastic.create_properties(self, model_part)
        properties.SetValue(KratosMultiphysics.YIELD_STRESS, self.yield_stress)
        properties.SetValue(KratosMultiphysics.REFERENCE_HARDENING_MODULUS, self.reference_hardening_modulus)
        properties.SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, self.isotropic_hardening_modulus)
        properties.SetValue(KratosMultiphysics.INFINITY_HARDENING_MODULUS, self.infinity_hardening_modulus)
        properties.SetValue(KratosMultiphysics.HARDENING_EXPONENT, self.hardening_exponent)
        return properties

class LinearJ2Plasticity3D(LinearJ2Plasticity):
    def __init__(self):
        LinearJ2Plasticity.__init__(self)
        self.dim = 3

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearJ2Plasticity3DLaw()


class Deformation():
    def __init__(self):
        self.nr_timesteps = 100

    def get_init_deformation_gradientF(self):
        self.F = KratosMultiphysics.Matrix(self.cl.dim,self.cl.dim)
        for i in range(self.cl.dim):
            for j in range(self.cl.dim):
                if(i==j):
                    self.F[i,j] = 1.0
                else:
                    self.F[i,j] = 0.0
        return self.F

    def initialize_reference_stress(self, strain_size):
        self.reference_stress = KratosMultiphysics.Vector(strain_size)
        for i in range(strain_size):
            self.reference_stress[i] = 0.0

    def set_deformation(self, cl_params, i):
        F = self.get_deformation_gradientF(i)
        detF = self.get_determinantF(i)
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(detF)


class DeformationLinearJ2Plasticity(Deformation):
    def __init__(self):
        Deformation.__init__(self)
        self.nr_timesteps = 10

    def set_deformation(self, cl_params, i):
        self.strain = (i+1)/ self.nr_timesteps * self.initial_strain
        cl_params.SetStrainVector(self.strain)

class DeformationLinearJ2Plasticity3D(DeformationLinearJ2Plasticity):
    def __init__(self):
        DeformationLinearJ2Plasticity.__init__(self)
        self.cl = LinearJ2Plasticity3D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.0
        self.initial_strain[3] = 0.001
        self.initial_strain[4] = 0.0
        self.initial_strain[5] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 4.03846; r_stress[0][1] = 4.03846; r_stress[0][2] = 2.42308; r_stress[0][3] = 0.80769; r_stress[0][4] = 0.0; r_stress[0][5] = 0.80769
        r_stress[1][0] = 8.07692; r_stress[1][1] = 8.07692; r_stress[1][2] = 4.84615; r_stress[1][3] = 1.61538; r_stress[1][4] = 0.0; r_stress[1][5] = 1.61538
        r_stress[2][0] = 11.6595; r_stress[2][1] = 11.6595; r_stress[2][2] = 8.18099; r_stress[2][3] = 1.73926; r_stress[2][4] = 0.0; r_stress[2][5] = 1.73926
        r_stress[3][0] = 15.1595; r_stress[3][1] = 15.1595; r_stress[3][2] = 11.681 ; r_stress[3][3] = 1.73926; r_stress[3][4] = 0.0; r_stress[3][5] = 1.73926
        r_stress[4][0] = 18.6595; r_stress[4][1] = 18.6595; r_stress[4][2] = 15.181 ; r_stress[4][3] = 1.73926; r_stress[4][4] = 0.0; r_stress[4][5] = 1.73926
        r_stress[5][0] = 22.1595; r_stress[5][1] = 22.1595; r_stress[5][2] = 18.681 ; r_stress[5][3] = 1.73927; r_stress[5][4] = 0.0; r_stress[5][5] = 1.73927
        r_stress[6][0] = 25.6595; r_stress[6][1] = 25.6595; r_stress[6][2] = 22.181 ; r_stress[6][3] = 1.73927; r_stress[6][4] = 0.0; r_stress[6][5] = 1.73927
        r_stress[7][0] = 29.1595; r_stress[7][1] = 29.1595; r_stress[7][2] = 25.681 ; r_stress[7][3] = 1.73928; r_stress[7][4] = 0.0; r_stress[7][5] = 1.73928
        r_stress[8][0] = 32.6595; r_stress[8][1] = 32.6595; r_stress[8][2] = 29.181 ; r_stress[8][3] = 1.73928; r_stress[8][4] = 0.0; r_stress[8][5] = 1.73928
        r_stress[9][0] = 36.1595; r_stress[9][1] = 36.1595; r_stress[9][2] = 32.681; r_stress[9][3] = 1.73929; r_stress[9][4] = 0.0; r_stress[9][5] = 1.73929
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

#####################################################################

if __name__ == "__main__":
#    KratosUnittest.main()

#    # parse configuration file
#    conf = configparser.ConfigParser()
#    conf.read("reduced_bases.cfg")
#    nr_modes = int(conf['Parameters']['nr_active_modes'])
#
    # Define a model
    model_part = KratosMultiphysics.ModelPart("test")

    deformation_test = DeformationLinearJ2Plasticity3D()

    #self._generic_constitutive_law_test(model_part, deformation_test)
    _generic_constitutive_law_test(model_part, deformation_test)
#
#
#    model_part_rve = km.ModelPart("RVE")
#    node1 = model_part_rve.CreateNewNode(1,0.0,0.0,0.0)
#    geom = km.Triangle2D3(node1, node1, node1) # create point geom
#    Model = {"RVE" : model_part_rve}
#    materials_rve  = km.Parameters("""
#               {
#                   "Parameters": {
#                           "materials_filename": "materials_rve.json"
#                   }
#           }
#           """)
#    read_materials_process.Factory(materials_rve, Model)
#    #read_materials_process.ReadMaterialsProcess(Model, materials_rve)
#
#    # import rve_data json string
#    with open ("rve.json", "r") as myfile:
#        rve_data = km.Parameters(myfile.read())
#
#    cl = kmsr.RVELaw(model_part_rve, rve_data)
#    cl_clone = cl.Clone()
#
#    cl.Check(km.ModelPart("dummy").Properties[1], geom, model_part_rve.ProcessInfo)
#    cl.InitializeMaterial(km.ModelPart("dummy").Properties[1], geom, km.Vector(3))
#
#    nr_comp = cl.GetStrainSize()
#    # creation and init
#    init_strain_macro = km.Vector(nr_comp)
#    homog_stress = km.Vector(nr_comp)
#    homog_constit = km.Matrix(nr_comp, nr_comp)
#    # trajectory 31:
#    init_strain_macro[0] = 0.001
#    init_strain_macro[1] = 0.0
#    init_strain_macro[2] = 0.001
#    init_strain_macro[3] = 0.0
#    init_strain_macro[4] = 0.0
#    init_strain_macro[5] = 0.001
#
#    cl_params = km.ConstitutiveLawParameters()
#    cl_options = km.Flags()
#    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
#    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
#    cl_params.SetOptions(cl_options)
#    cl_params.SetStressVector(homog_stress)
#    cl_params.SetConstitutiveMatrix(homog_constit)
#    cl_params.SetMaterialProperties(model_part_rve.Properties[1])
#
#    nr_timesteps = 250
#    t = dt = 1. / nr_timesteps
#    fo = open("homog_stress.dat",'w')
#    while (t <= 1. + dt / 10.):
#        print("time {:.3f}".format(t))
#        model_part_rve.CloneTimeStep(t)
#        strain_macro = t * init_strain_macro
#        cl_params.SetStrainVector(strain_macro)
#        cl.CalculateMaterialResponseCauchy(cl_params)
#        cl.FinalizeSolutionStep(km.ModelPart("dummy").Properties[1], geom,
#                                km.Vector(3), model_part_rve.ProcessInfo)
#        # Print output
#        #modes_weights = km.Vector(nr_modes)
#        #print(cl.GetValue(kmsr.REDUCED_MODES_WEIGHTS, modes_weights))
#        cl_params.GetStressVector(homog_stress)
#        #homog_stress = cl_params.GetStressVector()
#        #print("{}: {}".format(t, homog_stress))
#        #cl_params.GetConstitutiveMatrix(homog_constit)
#        #print("{}: {}".format(t, homog_constit))
#        t += dt
#        fo.write("{}\n".format(homog_stress[0]))
#    fo.close()
