from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication as MultiscaleROMApplication
import os

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


class Output:
    def __init__(self, filename):
        try:
            os.remove(filename)
        except OSError:
            pass
        self.fo = open(filename, "w")
        self.fo.write("#    {:<78}{:<30}\n".format("Strain", "Stress"))
        self.fo.write("#1   "
                      "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                      "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
            "2", "3", "4", "5" , "6", "7",
            "8", "9", "10", "11" , "12", "13"))
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

    def print(self, i, cl_params):
        print("\nStep: {}".format(i))
        print("Strain = ", cl_params.GetStrainVector())
        print("Stress = ", cl_params.GetStressVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())

def generic_constitutive_law_test(model_part, deformation_test):
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

    output = Output("strain-stress.dat")
    for i in range(deformation_test.nr_timesteps):
        deformation_test.set_deformation(cl_params, i)

        # Chauchy
        cl.CalculateMaterialResponseCauchy(cl_params)
        cl.FinalizeMaterialResponseCauchy(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        #output.print(i, cl_params)
        output.write(i, cl_params)

        reference_stress = deformation_test.get_reference_stress(i)
        stress = cl_params.GetStressVector()


class RVELaw():
   def __init__(self, material_parameters):
       self.dim = 3
       self.material_parameters = material_parameters

   def create_properties(self, model_part):
       prop_id = 0
       properties = model_part.Properties[prop_id]
       return properties

   def create_constitutive_Law(self):
       return MultiscaleROMApplication.RVELaw().Create(self.material_parameters)


class Deformation():
    def __init__(self, parameters):
        self.nr_timesteps = parameters["nr_timesteps"]
        self.strain_list = parameters["strain"]

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


class DeformationRVELaw(Deformation):
    def __init__(self, parameters, material_parameters):
        Deformation.__init__(self, parameters)
        self.cl = RVELaw(material_parameters)

    def set_deformation(self, cl_params, i):
        self.strain = (i+1)/ self.nr_timesteps * self.initial_strain
        cl_params.SetStrainVector(self.strain)

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = self.strain_list[0]
        self.initial_strain[1] = self.strain_list[1]
        self.initial_strain[2] = self.strain_list[2]
        self.initial_strain[3] = self.strain_list[3]
        self.initial_strain[4] = self.strain_list[4]
        self.initial_strain[5] = self.strain_list[5]

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

#####################################################################
if __name__ == "__main__":

    parameters = {"nr_timesteps": 10,
                 #[XX, YY, ZZ, XY, YZ, XZ]
                  #"strain": [0.0, 0.0, 0.0, 0.0005, 0.0, 0.0]
                  "strain": [0.001, 0.001, 0.0, 0.001, 0.0, 0.001]
                  }

    materials_rve  = KratosMultiphysics.Parameters("""
    {
        "name": "RVE Law",
        "Parameters": {
           "rve_materials_filename": "materials_rve.json",
           "rve_data_filename": "rve.json",
           "convergence_criterion": "residual_criterion",
           "residual_relative_tolerance": 1e-2,
           "residual_absolute_tolerance": 1e-9,
           "max_iteration": 10,
           "verbose": 1
        }
    }""")

    model_part = KratosMultiphysics.ModelPart("test")

    # Test RVE
    deformation_test = DeformationRVELaw(parameters, materials_rve)
    generic_constitutive_law_test(model_part, deformation_test)
