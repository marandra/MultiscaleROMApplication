import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import numpy
import matplotlib.pyplot as plt


def VerifyTangent(E, DE, subdivisions=10):
    factors = []
    errors = []
    linear = []
    quadratic = []
    alpha = 2.0
    for i in range(subdivisions):
        alpha /= 2.0
        error_norm = _EstimateError(E, DE, alpha)
        factors.append(alpha)
        errors.append(error_norm)
        linear.append(alpha)
        quadratic.append(alpha * alpha)

    print("DEBUG:")
    print("factors:", factors)
    print("errors:", errors)
    print("linear:", linear)
    print("quadratic:", quadratic)

    plt.plot(factors, errors, "o-")
    plt.plot(factors, linear, "--")
    plt.plot(factors, quadratic, "--")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


def _EstimateError(E, dE, alpha):

    # here compute S(E)
    _, t_S, C = constitutive_law_kratos(E)
    S = t_S[-1]
    dE = numpy.array(dE)
    #calculate stress using linearized approach
    S_linearized = S + alpha * numpy.inner(C, dE)
    E_perturbed = E + alpha * dE
    #here compute S(Eperturbed) = S(E+alpha*DE)
    _, t_S, C = constitutive_law_kratos(E_perturbed)
    #_, t_S, C = constitutive_law_python(E_perturbed)
    S_perturbed = t_S[-1]

    S_error = S_perturbed - S_linearized
    S_error_norm = numpy.linalg.norm(S_error)

    return S_error_norm


def get_step_loads():
    ts_list = []
    load = PARAMETERS["load"]
    for x in range(len(load) - 1):
        for i in numpy.linspace(load[x + 0],
                                load[x + 1],
                                PARAMETERS["nr_timesteps"],
                                endpoint=False):
            ts_list.append(i)
    return ts_list


def constitutive_law_kratos(init_strain):

    node = MODEL_PART.CreateNewNode(1, 0.0, 0.0, 0.0)
    geom = KM.Tetrahedra3D4(node, node, node, node)
    cl_options = KM.Flags()
    cl_options.Set(KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
    cl_options.Set(KM.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
    stress_vector = KM.Vector(CL.GetStrainSize())
    strain_vector = KM.Vector(CL.GetStrainSize())
    constitutive_matrix = KM.Matrix(CL.GetStrainSize(), CL.GetStrainSize())
    cl_params = KM.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetProcessInfo(MODEL_PART.ProcessInfo)
    cl_params.SetMaterialProperties(PROPERTIES)
    #cl_params.SetElementGeometry(geom)
    #F = KM.Matrix(3, 3, 0.0)
    #F[0, 0] = F[1, 1] = F[2, 2] = 1
    #cl_params.SetDeformationGradientF(F)
    #cl_params.SetDeterminantF(1)

    CL.InitializeMaterial(PROPERTIES, geom, KM.Vector(4))

    t_strain = numpy.array([0, 0, 0, 0, 0, 0])
    t_stress = numpy.array([0, 0, 0, 0, 0, 0])
    for mult in get_step_loads():
        s = KM.Vector(list(mult * init_strain))

        cl_params.SetStrainVector(s)
        CL.InitializeMaterialResponseCauchy(cl_params)

        cl_params.SetStrainVector(s)
        CL.CalculateMaterialResponseCauchy(cl_params)

        cl_params.SetStrainVector(s)
        CL.FinalizeMaterialResponseCauchy(cl_params)

        strain = cl_params.GetStrainVector()
        stress = cl_params.GetStressVector()
        tensor = cl_params.GetConstitutiveMatrix()

        t_strain = numpy.vstack((t_strain, strain))
        t_stress = numpy.vstack((t_stress, stress))

    n_tensor = numpy.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            n_tensor[i, j] = tensor[i, j]

    return t_strain, t_stress, n_tensor


def constitutive_law_python(init_strain):
    def hardening_law(r):
        r0 = PROPERTIES[KM.YIELD_STRESS] / numpy.sqrt(
            PROPERTIES[KM.YOUNG_MODULUS])
        H0 = PROPERTIES[SMA.HARDENING_MODULI_VECTOR][0]
        return r0 + H0 * (r - r0)

    def compute_elastic_tensor():
        e = PROPERTIES[KM.YOUNG_MODULUS]
        v = PROPERTIES[KM.POISSON_RATIO]
        c1 = e / ((1.00 + v) * (1 - 2 * v))
        c2 = c1 * (1 - v)
        c3 = c1 * v
        c4 = c1 * 0.5 * (1 - 2 * v)
        cm = numpy.zeros((6, 6))
        cm[0, 0] = c2
        cm[0, 1] = c3
        cm[0, 2] = c3
        cm[1, 0] = c3
        cm[1, 1] = c2
        cm[1, 2] = c3
        cm[2, 0] = c3
        cm[2, 1] = c3
        cm[2, 2] = c2
        cm[3, 3] = c4
        cm[4, 4] = c4
        cm[5, 5] = c4
        return cm

    def calculate_material_response(E):
        H = PROPERTIES[SMA.HARDENING_MODULI_VECTOR][0]
        c = compute_elastic_tensor()
        s_eff = numpy.inner(c, E)
        tao_eps = numpy.sqrt((numpy.inner(E, s_eff)))
        if tao_eps <= R:
            r = R
            q = hardening_law(r)
            d = 1 - (q / r)
            s = (1 - d) * s_eff
            c *= (1 - d)
        else:
            r = tao_eps
            q = hardening_law(r)
            d = 1 - (q / r)
            s = (1 - d) * s_eff
            c *= (1 - d)
            c -= (q - H * r) / (r * r * r) * (s_eff @ s_eff.transpose())
        return r, s, c

    # Main loop
    R = PROPERTIES[KM.YIELD_STRESS] / numpy.sqrt(PROPERTIES[KM.YOUNG_MODULUS])
    t_strain = numpy.array([0, 0, 0, 0, 0, 0])
    t_stress = numpy.array([0, 0, 0, 0, 0, 0])
    for mult in get_step_loads():
        strain = mult * init_strain
        R, stress, tensor = calculate_material_response(strain)
        t_strain = numpy.vstack((t_strain, strain))
        t_stress = numpy.vstack((t_stress, stress))
    return t_strain, t_stress, tensor


#####################################################################
#  main
#####################################################################

if __name__ == "__main__":

    MODEL_PART = KM.Model().CreateModelPart("test")

    PROPERTIES = MODEL_PART.Properties[0]
    PROPERTIES.SetValue(KM.YOUNG_MODULUS, 3000)
    PROPERTIES.SetValue(KM.POISSON_RATIO, 0.3)
    PROPERTIES.SetValue(KM.YIELD_STRESS, 0.5)
    PROPERTIES.SetValue(SMA.INFINITY_YIELD_STRESS, 0.7)
    PROPERTIES.SetValue(SMA.HARDENING_MODULI_VECTOR, [0.05, 0.05])

    #CL = CLA.SmallStrainIsotropicDamageTractionOnly3DLaw()
    CL = CLA.SmallStrainIsotropicDamage3DLaw()

    PARAMETERS = {
        "nr_timesteps": 10,
        #"load": 0.005 * numpy.array([0, 1, 0, -1.2, 0, 1.4, 0, -1.6]),
        "load": 0.000_3 * numpy.array([0, 1]),
        "strain": numpy.array([
            -0.0759,
            0.7483,
            0.1879,
            0.5391,
            0.0063,
            -0.3292,
        ])
    }

    t_strain, t_stress, tensor = constitutive_law_kratos(PARAMETERS["strain"])
    plt.plot(t_strain, t_stress, "go", label="KRATOS")
    t_strain, t_stress, tensor = constitutive_law_python(PARAMETERS["strain"])
    plt.plot(t_strain, t_stress, "c--", label="PYTHON")
    plt.show()

    dstrain = PARAMETERS["load"][-1] / PARAMETERS["nr_timesteps"] * PARAMETERS[
        "strain"]
    VerifyTangent(PARAMETERS["strain"], dstrain, subdivisions=10)
