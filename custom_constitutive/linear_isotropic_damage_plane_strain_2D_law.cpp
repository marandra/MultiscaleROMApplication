#include "linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
// CONSTRUCTOR
LinearIsotropicDamagePlaneStrain2DLaw::LinearIsotropicDamagePlaneStrain2DLaw()
    : ConstitutiveLaw()
{
}

// CLONE
ConstitutiveLaw::Pointer LinearIsotropicDamagePlaneStrain2DLaw::Clone() const
{
    return ConstitutiveLaw::Pointer(new LinearIsotropicDamagePlaneStrain2DLaw());
}

LinearIsotropicDamagePlaneStrain2DLaw::SizeType LinearIsotropicDamagePlaneStrain2DLaw::WorkingSpaceDimension()
{
    return 2;
}

LinearIsotropicDamagePlaneStrain2DLaw::SizeType LinearIsotropicDamagePlaneStrain2DLaw::GetStrainSize()
{
    return 3;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    return false;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::Has(const Variable<array_1d<double, 2>>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    return false;
}

double& LinearIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<double>& rThisVariable,
                                                        double& rValue)
{
    if(rThisVariable == STRAIN_ENERGY){
        rValue = mStrainEnergy;
    }
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }
    return rValue;
}

Vector& LinearIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Vector>& rThisVariable,
                                                        Vector& rValue)
{
    return rValue;
}

Matrix& LinearIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Matrix>& rThisVariable,
                                                        Matrix& rValue)
{
    return rValue;
}

array_1d<double, 2>& LinearIsotropicDamagePlaneStrain2DLaw::GetValue(
    const Variable<array_1d<double, 2>>& rVariable, array_1d<double, 2>& rValue)
{
    return rValue;
}

array_1d<double, 3>& LinearIsotropicDamagePlaneStrain2DLaw::GetValue(
    const Variable<array_1d<double, 3>>& rVariable, array_1d<double, 3>& rValue)
{
    return rValue;
}

void LinearIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<double>& rVariable,
                                                     const double& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Vector>& rVariable,
                                                     const Vector& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Matrix>& rVariable,
                                                     const Matrix& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 2>>& rVariable,
                                                     const array_1d<double, 2>& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 3>>& rVariable,
                                                     const array_1d<double, 3>& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

bool LinearIsotropicDamagePlaneStrain2DLaw::ValidateInput(const Properties& rMaterialProperties)
{
    return true;
}

LinearIsotropicDamagePlaneStrain2DLaw::StrainMeasure LinearIsotropicDamagePlaneStrain2DLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

LinearIsotropicDamagePlaneStrain2DLaw::StressMeasure LinearIsotropicDamagePlaneStrain2DLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool LinearIsotropicDamagePlaneStrain2DLaw::IsIncremental()
{
    return false;
}

void LinearIsotropicDamagePlaneStrain2DLaw::InitializeMaterial(const Properties& material_prop,
                                                               const GeometryType& rElementGeometry,
                                                               const Vector& rShapeFunctionsValues)
{
    r_prev = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    tau_epsilon = 0.;
}

void LinearIsotropicDamagePlaneStrain2DLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // update of damage threshold
    r_prev = r;
}

void LinearIsotropicDamagePlaneStrain2DLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& matprops = rValues.GetMaterialProperties();
    Vector& epsilon = rValues.GetStrainVector();
    Vector& sigma_bar = rValues.GetStressVector();
    Vector sigma_bar_pos;
    Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
    double H = matprops[ISOTROPIC_DAMAGE_MODULUS];
    double dpointcoeff;
    double d, q;
    double sigma_xx, sigma_yy, sigma_xy;
    double hyp, sigma_1, sigma_2, angle, cos_a, sin_a;
    bool TRACTION_ONLY = matprops[FLOW_RULE_IS_TRACTION_ONLY];

    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN_VECTOR))
    {
        noalias(epsilon) += rValues.GetProcessInfo()[INITIAL_STRAIN_VECTOR];
    }

    CalculateConstitutiveMatrix(matprops, constitutiveMatrix);
    sigma_bar = prod(constitutiveMatrix, epsilon);
    sigma_bar_pos = prod(constitutiveMatrix, epsilon);

    // for tension-only fluency law:
    // originally sigma and sigma_positive are the same (as it is in the
    // symmetrical case), this block modifies sigma_positive
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
        r = r_prev;
        q = CalculateQ(r, matprops);
        d = 1. - q / r;
        constitutiveMatrix *= (1 - d);
        sigma_bar *= (1 - d);
    }
    else
    {

        mInelasticFlag = 1;
        r = tau_epsilon;
        q = CalculateQ(r, matprops);
        d = 1. - q / r;
        dpointcoeff = (q - H * r) / (r * r * r);
        constitutiveMatrix *= (1. - d);
        constitutiveMatrix -= dpointcoeff * outer_prod(sigma_bar_pos, sigma_bar);
        sigma_bar *= (1. - d);
    }
  //TODO compute mStrainEnergy here
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

void LinearIsotropicDamagePlaneStrain2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues)
{
}

double LinearIsotropicDamagePlaneStrain2DLaw::CalculateQ(double r, const Properties& material_prop)
{
    double H = material_prop[ISOTROPIC_DAMAGE_MODULUS];
    double r0 = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    double q_inf = material_prop[INFINITY_YIELD_STRESS] /
                   std::sqrt(material_prop[YOUNG_MODULUS]);
    double q;

    if (r < r0)
        return r;
    q = r0 + H * (r - r0);
    if ((H > 0 && q > q_inf) || (H < 0 && q < q_inf))
        q = q_inf;
    return q;
}

void LinearIsotropicDamagePlaneStrain2DLaw::CalculateConstitutiveMatrix(const Properties& props,
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

void LinearIsotropicDamagePlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

int LinearIsotropicDamagePlaneStrain2DLaw::Check(const Properties& rMaterialProperties,
                                                 const GeometryType& rElementGeometry,
                                                 const ProcessInfo& rCurrentProcessInfo)
{
    if (!rMaterialProperties.Has(YOUNG_MODULUS))
        KRATOS_THROW_ERROR(
            std::invalid_argument,
            "LinearIsotropicDamagePlaneStrain2DLaw - missing YOUNG_MODULUS", "");
    if (!rMaterialProperties.Has(POISSON_RATIO))
        KRATOS_THROW_ERROR(
            std::invalid_argument,
            "LinearIsotropicDamagePlaneStrain2DLaw - missing POISSON_RATIO", "");
    if (!rMaterialProperties.Has(INFINITY_YIELD_STRESS))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "LinearIsotropicDamagePlaneStrain2DLaw - missing "
                           "INFINITY_YIELD_STRESS",
                           "");
    if (rMaterialProperties[INFINITY_YIELD_STRESS] < 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "LinearIsotropicDamagePlaneStrain2DLaw - "
                           "INFINITY_YIELD_STRESS must be positive",
                           "");
    if (!rMaterialProperties.Has(ISOTROPIC_DAMAGE_MODULUS))
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "LinearIsotropicDamagePlaneStrain2DLaw - missing "
                           "ISOTROPIC_DAMAGE_MODULUS",
                           "");
    // std::cout << "WARNING: ISOTROPIC_DAMAGE_MODULUS check deactivated" <<
    // std::endl;
    // std::cout << "         Fix bug in multiscale application." << std::endl;
    if (rMaterialProperties[ISOTROPIC_DAMAGE_MODULUS] >= 1.)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "LinearIsotropicDamagePlaneStrain2DLaw - "
                           "ISOTROPIC_DAMAGE_MODULUS must be < 1.",
                           "");
    return 0;
}

} /* namespace Kratos.*/
