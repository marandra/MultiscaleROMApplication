#include "exponential_isotropic_damage_plane_strain_2D_law.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
// CONSTRUCTOR
ExponentialIsotropicDamagePlaneStrain2DLaw::ExponentialIsotropicDamagePlaneStrain2DLaw()
    : ConstitutiveLaw()
      //, m_initialized(false)
      ,
      m_init_strain()
{
}

// CLONE
ConstitutiveLaw::Pointer ExponentialIsotropicDamagePlaneStrain2DLaw::Clone() const
{
    return ConstitutiveLaw::Pointer(new ExponentialIsotropicDamagePlaneStrain2DLaw());
}

ExponentialIsotropicDamagePlaneStrain2DLaw::SizeType ExponentialIsotropicDamagePlaneStrain2DLaw::WorkingSpaceDimension()
{
    return 2;
}

ExponentialIsotropicDamagePlaneStrain2DLaw::SizeType ExponentialIsotropicDamagePlaneStrain2DLaw::GetStrainSize()
{
    return 3;
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<double>& rThisVariable)
{
    return false;
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    return false;
}

double& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<double>& rThisVariable,
                                                             double& rValue)
{
    return rValue;
}

Vector& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Vector>& rThisVariable,
                                                             Vector& rValue)
{
    if (rThisVariable == INITIAL_STRAIN)
    {
        if (rValue.size() != m_init_strain.size())
            rValue.resize(m_init_strain.size(), false);
        noalias(rValue) = m_init_strain;
    }
    return (rValue);
}

Matrix& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Matrix>& rThisVariable,
                                                             Matrix& rValue)
{
    return rValue;
}

array_1d<double, 3>& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(
    const Variable<array_1d<double, 3>>& rVariable, array_1d<double, 3>& rValue)
{
    return rValue;
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<double>& rVariable,
                                                          const double& rValue,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Vector>& rVariable,
                                                          const Vector& rValue,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == INITIAL_STRAIN)
    {
        if (rValue.size() == m_init_strain.size())
            noalias(m_init_strain) = rValue;
    }
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Matrix>& rVariable,
                                                          const Matrix& rValue,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(
    const Variable<array_1d<double, 3>>& rVariable,
    const array_1d<double, 3>& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::ValidateInput(const Properties& rMaterialProperties)
{
    return true;
}

ExponentialIsotropicDamagePlaneStrain2DLaw::StrainMeasure ExponentialIsotropicDamagePlaneStrain2DLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

ExponentialIsotropicDamagePlaneStrain2DLaw::StressMeasure ExponentialIsotropicDamagePlaneStrain2DLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool ExponentialIsotropicDamagePlaneStrain2DLaw::IsIncremental()
{
    return false;
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::InitializeMaterial(
    const Properties& material_prop, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
{
    r_prev = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    tau_epsilon = 0.;
    m_init_strain = ZeroVector(this->GetStrainSize());
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // update of damage threshold
    r_prev = r;
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& matprops = rValues.GetMaterialProperties();
    Vector& epsilon = rValues.GetStrainVector();
    Vector& sigma_bar = rValues.GetStressVector();
    Vector sigma_bar_pos;
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    double H = matprops[ISOTROPIC_DAMAGE_MODULUS];
    double dpointcoeff;
    double d, q;
    double sigma_xx, sigma_yy, sigma_xy;
    double hyp, sigma_1, sigma_2, angle, cos_a, sin_a;
    bool TRACTION_ONLY = matprops[FLOW_RULE_IS_TRACTION_ONLY];

    noalias(epsilon) -= m_init_strain;
    CalculateConstitutiveMatrix(matprops, constitutive_matrix);
    sigma_bar = prod(constitutive_matrix, epsilon);
    sigma_bar_pos = prod(constitutive_matrix, epsilon);
    //TODO: for use with strain energy computation (see below)
    Matrix constitutive_elastic_matrix = constitutive_matrix;

    // Originally sigma and sigma_positive are the same (symmetrical case).
    // In case of tension-only fluency law, the following block modifies
    // sigma_positive.
    if (TRACTION_ONLY)
    {
        sigma_xx = sigma_bar(0);
        sigma_yy = sigma_bar(1);
        sigma_xy = sigma_bar(2);
        hyp = std::hypot(0.5 * (sigma_xx - sigma_yy), sigma_xy);
        sigma_1 = 0.5 * (sigma_xx + sigma_yy) + hyp;
        sigma_2 = 0.5 * (sigma_xx + sigma_yy) - hyp;
        angle = 0.5 * std::atan2(2.0 * sigma_xy, sigma_xx - sigma_yy);
        cos_a = std::cos(angle);
        sin_a = std::sin(angle);
        sigma_bar_pos(0) = 0.0;
        sigma_bar_pos(1) = 0.0;
        sigma_bar_pos(2) = 0.0;
        if (sigma_1 > 0)
        {
            sigma_bar_pos(0) += sigma_1 * cos_a * cos_a;
            sigma_bar_pos(1) += sigma_1 * sin_a * sin_a;
            sigma_bar_pos(2) += sigma_1 * sin_a * cos_a;
        }
        if (sigma_2 > 0)
        {
            sigma_bar_pos(0) += sigma_2 * sin_a * sin_a;
            sigma_bar_pos(1) += sigma_2 * cos_a * cos_a;
            sigma_bar_pos(2) -= sigma_2 * sin_a * cos_a;
        }
    }

    tau_epsilon = std::sqrt(inner_prod(sigma_bar_pos, epsilon));

    if (tau_epsilon <= r_prev)
    {
        // std::cout << "ELASTIC regime" << std::endl;
        r = r_prev;
        q = CalculateQ(r, matprops);
        d = 1. - q / r;
        constitutive_matrix *= (1 - d);
        sigma_bar *= (1 - d);
    }
    else
    {
        // std::cout << "INELASTIC regime" << std::endl;
        r = tau_epsilon;
        q = CalculateQ(r, matprops);
        d = 1. - q / r;
        dpointcoeff = (q - H * r) / (r * r * r);
        constitutive_matrix *= (1. - d);
        constitutive_matrix -= dpointcoeff * outer_prod(sigma_bar_pos, sigma_bar);
        sigma_bar *= (1. - d);
    }

    //TODO add check of flag here (COMPUTE_STRAIN_ENERGY)
    mStrainEnergy = 0.5 * ((1. - d) * inner_prod(epsilon, prod(constitutive_elastic_matrix, epsilon)));
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                               const GeometryType& rElementGeometry,
                                                               const Vector& rShapeFunctionsValues)
{
}

double ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateQ(double r, const Properties& material_prop)
{
    double H = material_prop[ISOTROPIC_DAMAGE_MODULUS];
    double r0 = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    double q_inf = material_prop[INFINITY_YIELD_STRESS] /
                   std::sqrt(material_prop[YOUNG_MODULUS]);
    double q;

    if (r < r0)
        return r;
    q = q_inf - (q_inf - r0) * std::exp(H * (1 - r / r0));
    return q;
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateConstitutiveMatrix(const Properties& props,
                                                                             Matrix& D)
{
    double E = props[YOUNG_MODULUS];
    double nu = props[POISSON_RATIO];
    double Ebar = E / (1. - nu * nu);
    double nubar = nu / (1. - nu);

    D.clear();

    D(0, 0) = 1;
    D(0, 1) = nubar;
    D(0, 2) = 0;
    D(1, 0) = nubar;
    D(1, 1) = 1;
    D(1, 2) = 0;
    D(2, 0) = 0;
    D(2, 1) = 0;
    D(2, 2) = 0.5 * (1 - nubar);

    D *= Ebar / (1. - nubar * nubar);
}

void ExponentialIsotropicDamagePlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

int ExponentialIsotropicDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties,
                                                      const GeometryType& rElementGeometry,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    if (!rMaterialProperties.Has(YOUNG_MODULUS))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ExponentialIsotropicDamagePlaneStrain2DLaw - "
                           "missing YOUNG_MODULUS",
                           "");
    if (!rMaterialProperties.Has(POISSON_RATIO))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ExponentialIsotropicDamagePlaneStrain2DLaw - "
                           "missing POISSON_RATIO",
                           "");
    if (!rMaterialProperties.Has(INFINITY_YIELD_STRESS))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ExponentialIsotropicDamagePlaneStrain2DLaw - "
                           "missing INFINITY_YIELD_STRESS",
                           "");
    if (rMaterialProperties[INFINITY_YIELD_STRESS] < 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ExponentialIsotropicDamagePlaneStrain2DLaw - "
                           "INFINITY_YIELD_STRESS must be positive",
                           "");
    if (!rMaterialProperties.Has(ISOTROPIC_DAMAGE_MODULUS))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ExponentialIsotropicDamagePlaneStrain2DLaw - "
                           "missing ISOTROPIC_DAMAGE_MODULUS",
                           "");
    return 0;
}

} /* namespace Kratos.*/
