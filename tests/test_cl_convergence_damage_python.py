import numpy
import matplotlib.pyplot as plt
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
from typing_extensions import final

numpy.set_printoptions(precision=12,
                       threshold=None,
                       edgeitems=None,
                       linewidth=120,
                       suppress=True,
                       nanstr=None,
                       infstr=None,
                       formatter=None,
                       sign=None,
                       floatmode=None,
                       legacy=None)


def get_convergence(final_strain, cl):
    def compute_error(strain, dstrain):
        #_, stresses, tensor = constitutive_law_kratos(strain)
        _, stresses, tensor = cl(strain)
        stress = stresses[-1]
        stress_linearized = stress + numpy.inner(tensor, dstrain)
        _, stresses, _ = constitutive_law_kratos(strain + dstrain)
        stress_perturbated = stresses[-1]
        return numpy.linalg.norm(stress_perturbated - stress_linearized)

    steps = []
    diffs = []
    for step, strain in enumerate(get_strains_cycle(final_strain)):
        alpha = 2
        nliter = 0
        diff = 999
        while diff > 1e-7 and nliter < 10:
            alpha /= 2.0
            nliter += 1
            diff = compute_error(strain, 0.1 * strain * alpha)
            diffs.append(diff)
            steps.append(step)

    return steps, diffs


def test_convergence(final_strain, cl, stop_at=-1):
    def _compute_error(_strain, _dstrain):
        _, stresses, _tensor = cl(_strain)
        _stress = stresses[-1]
        stress_linearized = _stress + numpy.inner(_tensor, _dstrain)
        _, stresses, _ = constitutive_law_kratos(_strain + _dstrain)
        stress_perturbated = stresses[-1]
        return numpy.linalg.norm(stress_perturbated - stress_linearized)

    for iter, strain in enumerate(get_strains_cycle(final_strain)):
        if iter == stop_at:
            break

    alphas = []
    errors = []
    alpha = 2

    for i in range(5):
        alpha /= 2.0
        diff = _compute_error(strain, 0.1 * strain * alpha)
        alphas.append(alpha)
        errors.append(diff)

    return alphas, errors


def get_strains_cycle(strain):
    factors = []
    load = PARAMETERS["load"]
    for x in range(len(load) - 1):
        for i in numpy.linspace(load[x + 0],
                                load[x + 1],
                                PARAMETERS["nr_timesteps"],
                                endpoint=True):
            factors.append(i)
    return [strain * f for f in factors]


def get_load_factors():
    factors = []
    load = PARAMETERS["load"]
    for x in range(len(load) - 1):
        for i in numpy.linspace(load[x + 0],
                                load[x + 1],
                                PARAMETERS["nr_timesteps"],
                                endpoint=True):
            factors.append(i)
    return factors


def constitutive_law_kratos(init_strain):

    N = KM.Vector(4)
    node1 = MODEL_PART.CreateNewNode(1, 0.0, 0.0, 0.0)
    node2 = MODEL_PART.CreateNewNode(2, 1.0, 0.0, 0.0)
    node3 = MODEL_PART.CreateNewNode(3, 0.0, 1.0, 0.0)
    node4 = MODEL_PART.CreateNewNode(4, 0.0, 0.0, 1.0)
    geom = KM.Tetrahedra3D4(node1, node2, node3, node4)
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
    cl_params.SetElementGeometry(geom)
    #F = KM.Matrix(3, 3, 0.0)
    #F[0, 0] = F[1, 1] = F[2, 2] = 1
    #cl_params.SetDeformationGradientF(F)
    #cl_params.SetDeterminantF(1)

    # debugging Plasticity
    #ProcessInfo = MODEL_PART.ProcessInfo
    #ProcessInfo[KM.NL_ITERATION_NUMBER] = 1
    #ProcessInfo[KM.STEP] = 2
    CL.InitializeMaterial(PROPERTIES, geom, N)

    strains = numpy.array([0, 0, 0, 0, 0, 0])
    stresses = numpy.array([0, 0, 0, 0, 0, 0])
    tensors = numpy.zeros((6, 6))
    for factor in get_load_factors():
        s = KM.Vector(list(factor * init_strain))

        cl_params.SetStrainVector(s)
        CL.InitializeMaterialResponseCauchy(cl_params)

        cl_params.SetStrainVector(s)
        CL.CalculateMaterialResponseCauchy(cl_params)

        cl_params.SetStrainVector(s)
        CL.FinalizeMaterialResponseCauchy(cl_params)

        strain_km = cl_params.GetStrainVector()
        stress_km = cl_params.GetStressVector()
        tensor_km = cl_params.GetConstitutiveMatrix()

        strains = numpy.vstack((strains, strain_km))
        stresses = numpy.vstack((stresses, stress_km))
        tensor = numpy.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            tensor[i, j] = tensor_km[i, j]

    return strains, stresses, tensor


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

    # CL main loop
    R = PROPERTIES[KM.YIELD_STRESS] / numpy.sqrt(PROPERTIES[KM.YOUNG_MODULUS])
    strains = numpy.array([0, 0, 0, 0, 0, 0])
    stresses = numpy.array([0, 0, 0, 0, 0, 0])
    for factor in get_load_factors():
        strain = factor * init_strain
        R, stress_km, tensor_km = calculate_material_response(strain)
        strains = numpy.vstack((strains, strain))
        stresses = numpy.vstack((stresses, stress_km))
    tensor = numpy.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            tensor[i, j] = tensor_km[i, j]

    return strains, stresses, tensor


#####################################################################
#  main
#####################################################################

if __name__ == "__main__":

    MODEL_PART = KM.Model().CreateModelPart("test")
    PROPERTIES = MODEL_PART.Properties[0]

    # Settings:

    #     CL = CLA.SmallStrainIsotropicPlasticity3DVonMisesVonMises()
    #     PROPERTIES.SetValue(KM.DENSITY, 7850.0)
    #     PROPERTIES.SetValue(KM.YOUNG_MODULUS, 201000000000.0)
    #     PROPERTIES.SetValue(KM.POISSON_RATIO, 0.263)
    #     PROPERTIES.SetValue(KM.FRACTURE_ENERGY, 10e6)
    #     PROPERTIES.SetValue(KM.YIELD_STRESS, 826958174.5230)
    #     PROPERTIES.SetValue(SMA.HARDENING_CURVE, 3)
    #     PROPERTIES.SetValue(SMA.TANGENCY_REGION2, False)
    #     PROPERTIES.SetValue(SMA.CURVE_FITTING_PARAMETERS, [
    #         826958174.5230, 13615516247.75, -2662617692352.00, 270525894696960.00,
    #         -14210125290012700.00, 370310601181757000.00, -3.78438E+18
    #     ])
    #     PROPERTIES.SetValue(SMA.PLASTIC_STRAIN_INDICATORS, [0.0029, 0.0085642])
    #
    #     # Load steps
    #
    #     PARAMETERS = {
    #         "nr_timesteps":
    #         5,
    #         #"load": 0.005 * numpy.array([0, 1, 0, -1.2, 0, 1.4, 0, -1.6]),
    #         "load":
    #         numpy.array([0, 1]),
    #         "strain":
    #         0.0003 * numpy.array([
    #             -0.0759,
    #             0.7483,
    #             0.1879,
    #             0.5391,
    #             0.0063,
    #             -0.3292,
    #         ])
    #     }
    #

    CL = CLA.SmallStrainIsotropicDamage3DLaw()
    PROPERTIES.SetValue(KM.YOUNG_MODULUS, 3000)
    PROPERTIES.SetValue(KM.POISSON_RATIO, 0.2)
    PROPERTIES.SetValue(KM.YIELD_STRESS, 0.5)
    PROPERTIES.SetValue(SMA.INFINITY_YIELD_STRESS, 0.7)
    PROPERTIES.SetValue(SMA.HARDENING_MODULI_VECTOR, [0.1, 0.1])

    PARAMETERS = {
        "nr_timesteps":
        80,
        #"load": 0.005 * numpy.array([0, 1, 0, -1.2, 0, 1.4, 0, -1.6]),
        "load":
        numpy.array([0, 1]),
        "strain":
        0.0005 * numpy.array([
            -0.0759,
            0.7483,
            0.1879,
            0.5391,
            0.0063,
            -0.3292,
        ])
    }

    cl = constitutive_law_kratos
    cl = constitutive_law_python
    stop_at = 35

    strainsp, stressesp, _ = constitutive_law_kratos(PARAMETERS["strain"])
    strains, stresses, _ = constitutive_law_kratos(PARAMETERS["strain"])
    factors, errors = test_convergence(PARAMETERS["strain"],
                                       cl,
                                       stop_at=stop_at)
    x, y = get_convergence(PARAMETERS["strain"], cl)

    # Create plot

    plt.figure(1)

    plt.subplot(311)
    plt.plot(strains, stresses, "x", label="KRATOS")
    plt.plot(strainsp, stressesp, "-", label="Python")
    plt.xlabel("strain")
    plt.ylabel("stress")

    plt.subplot(312)
    plt.plot(x, [1e-7] * len(x), ":r")
    for i in range(len(x)):
        if x[i] == stop_at:
            plt.scatter(x[i], y[i], marker="x", color="r")
            continue
        plt.scatter(x[i], y[i], marker="x", color="b")
    #plt.plot(factors, factors, "--", factors, numpy.multiply(factors, factors),
    #         "--", factors, errors, "o-")
    #plt.xscale("log")
    plt.yscale("log")
    #plt.ylim((1e-8, 1e-5))
    plt.xlabel("time step")
    plt.ylabel("norm(diff)")

    plt.subplot(313)
    plt.plot(factors, errors[0] * numpy.array(factors), "--", factors,
             errors[0] * numpy.multiply(factors, factors), "--", factors,
             errors, "o-")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("non linear iteration")
    plt.ylabel("norm(diff)")

    plt.show()
