from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
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

    def write(self, cl_params, model_part):
        strain = cl_params.GetStrainVector()
        strain = model_part.ProcessInfo[km.INITIAL_STRAIN]
        stress = cl_params.GetStressVector()
        line = "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "\
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
            strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
            stress[0], stress[1], stress[2], stress[3], stress[4], stress[5])
        self.fo.write(line)


def constitutive_law_test(model_part, deformation, idx):

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

    case_params = deformation.cases[idx]

    # Construct a constitutive law
    rve_params = deformation.rve_parameters

    #tmp_param = KratosMultiphysics.Parameters('''{
    #"rve_materials_filename": "ACA 11111",
    #"rve_data_filename": "ACA 22222"
    #}
    #''')

    rve_materials_filename = "{}/{}v{}.json".format(deformation.parameters["rve_materials_path"].GetString(),
                                            deformation.parameters["rve_materials_base_filename"].GetString(),
                                            case_params[3])
    rve_data_filename = "{}/{}{}m_{}ip.json".format(deformation.parameters["rve_data_path"].GetString(),
                                                deformation.parameters["rve_data_base_filename"].GetString(),
                                                case_params[1], case_params[2])
    tmp_param = KratosMultiphysics.Parameters('''{
    "rve_materials_filename": "''' + rve_materials_filename + '''",
    "rve_data_filename": "''' + rve_data_filename + '''"
    }''')

    rve_params["Parameters"].AddMissingParameters(tmp_param)
    cl = km.MultiscaleROMApplication.RVELaw().Create(rve_params)

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
    initial_strain_kratos = km.Vector(6)
    for i in range(6):
        initial_strain_kratos[i] = case_params[4][i]
    #deformation.set_initial_strain(cl.GetStrainSize())
    output = Output(case_params[5])

    print()
    #print("Has DAMAGE_VARIABLE: ", cl.Has(km.DAMAGE_VARIABLE))
    #print("Has STRAIN_ENERGY: ", cl.Has(km.STRAIN_ENERGY))
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
    nr_ts = deformation.nr_timesteps + 1
    load_cycle = deformation.load_cycle
    ts_array = [numpy.linspace(load_cycle[x+0], load_cycle[x+1], nr_ts) for x in range(len(load_cycle)-1)]
    ts_list = []
    for x in ts_array:
        ts_list.extend(x.tolist())
    for strain_mult in ts_list:
        deformation.set_deformation(cl_params, strain_mult, initial_strain_kratos)

        # Chauchy
        model_part.ProcessInfo[km.INITIAL_STRAIN] = cl_params.GetStrainVector()

        # Call sequence
        zero_vector = km.Vector(6)
        for i in range(6):
            zero_vector[i] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.InitializeMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6)
        for i in range(6):
            zero_vector[i] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.CalculateMaterialResponseCauchy(cl_params)

        zero_vector = km.Vector(6)
        for i in range(6):
            zero_vector[i] = 0.
        cl_params.SetStrainVector(zero_vector)
        cl.FinalizeMaterialResponseCauchy(cl_params)

        #output.printout(i, cl_params)
        output.write(cl_params, model_part)

        #debug = model_part.ProcessInfo[km.INITIAL_STRAIN]
        #print(debug)

        #reference_stress = deformation.get_reference_stress(i)
        stress = cl_params.GetStressVector()
        #print("Step ", i)
        #print("Reference: ", reference_stress)


class Deformation():
    def __init__(self, parameters):
        # Constant values
        self.nr_timesteps = parameters["nr_timesteps"].GetInt()
        self.load_cycle = parameters["load_cycle"].GetVector()
        self.rve_parameters = parameters["rve_parameters"]
        self.parameters = parameters

        self.cases = []
        list_modes = [int(x) for x in parameters["rve_data_modes"].GetVector()]
        list_ip = [int(x) for x in parameters["rve_data_points"].GetVector()]
        list_materials = [int(x) for x in parameters["rve_materials"].GetVector()]
        list_traj = [int(x) for x in parameters["trajectories"].GetVector()]
        for t in list_traj:
            for mat in list_materials:
                for i in list_ip:
                    for m in list_modes:
                        if t != -1:
                            filename = "{}/trajectory_{:02}/ProjectParameters.json".format(parameters["trajectories_path"].GetString(), t)
                            with open(filename, 'r') as fp:
                                for line in fp.readlines():
                                    if '''"initial_strain"''' in line:
                                        strain = [float(x) for x in line.split('[')[1].split(']')[0].split(',')]
                                        break
                        else:
                            strain = [x for x in parameters["strain"].GetVector()]
                        filename = "{}_T{:02}_{}m_{}ip_mat{}.dat".format(parameters["output_base_filename"].GetString(), t, m, i, mat)
                        case = [t, m, i, mat, strain, filename]
                        self.cases.append(case)
                        print (case)
        # To be defined inside case loop
        #self.initial_strain_list = strain
        self.nr_cases = len(self.cases)

    #def set_initial_strain(self, strain_size=6):
    #    self.initial_strain_kratos = km.Vector(strain_size)
    #    for i in range(strain_size):
    #        self.initial_strain_kratos[i] = self.initial_strain_list[i]

    def set_deformation(self, cl_params, mult, initial_strain_kratos):
        self.strain = mult * initial_strain_kratos
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


#####################################################################
if __name__ == "__main__":

    model_part = km.Model().CreateModelPart("test")

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    parameters_defaults = KratosMultiphysics.Parameters('''{
    "reuse_existing_files": true,
    "nr_timesteps": 20,
    "load_cycle": [0, 1, 0, -1.2, 0, 1.4, 0, -1.6],
    "trajectories_path": "../training",
    "trajectories": [-1],
    "strain": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "output_base_filename": "strain-stress-rve_traction.dat",
    "rve_materials_path": "rve_materials",
    "rve_materials_base_filename": "rve_material_",
    "rve_materials": [0],
    "rve_data_path": "../offline_data",
    "rve_data_base_filename": "rve_",
    "rve_data_modes": [10],
    "rve_data_points": [100],
    "rve_parameters": {
        "name": "RVE Law",
        "Parameters": {
             "rve_materials_filename": "to be filled by script",
             "rve_data_filename": "to be filled by script",
             "convergence_criterion": "displacement_criterion",
             "residual_relative_tolerance": 1e-6,
             "residual_absolute_tolerance": 0,
             "max_iteration": 10,
             "verbose": 0
        }
    }
}
''')
    parameters.ValidateAndAssignDefaults(parameters_defaults)

    deformation = Deformation(parameters)
    for i in range(deformation.nr_cases):
        constitutive_law_test(model_part, deformation, i)
