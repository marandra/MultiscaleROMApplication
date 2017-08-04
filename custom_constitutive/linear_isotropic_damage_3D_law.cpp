#include "linear_isotropic_damage_3D_law.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
// CONSTRUCTOR
LinearIsotropicDamage3DLaw::LinearIsotropicDamage3DLaw()
    : ConstitutiveLaw()
{
}

// COPY CONSTRUCTOR
LinearIsotropicDamage3DLaw::LinearIsotropicDamage3DLaw(const LinearIsotropicDamage3DLaw &rOther)
    : ConstitutiveLaw(rOther)
{
}

// CLONE
ConstitutiveLaw::Pointer LinearIsotropicDamage3DLaw::Clone() const
{
    return ConstitutiveLaw::Pointer(new LinearIsotropicDamage3DLaw());
}

LinearIsotropicDamage3DLaw::SizeType LinearIsotropicDamage3DLaw::WorkingSpaceDimension()
{
    return 3;
}

LinearIsotropicDamage3DLaw::SizeType LinearIsotropicDamage3DLaw::GetStrainSize()
{
    return 6;
}

bool LinearIsotropicDamage3DLaw::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    return false;
}

bool LinearIsotropicDamage3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamage3DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamage3DLaw::Has(const Variable<array_1d<double, 2>>& rThisVariable)
{
    return false;
}

bool LinearIsotropicDamage3DLaw::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    return false;
}

// New method -  for 3D cases
// bool LinearIsotropicDamage3DLaw::Has(const
// Variable<array_1d<double, 6 > >& rThisVariable)
//{
//	return false;
//}

double& LinearIsotropicDamage3DLaw::GetValue(const Variable<double>& rThisVariable,
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

Vector& LinearIsotropicDamage3DLaw::GetValue(const Variable<Vector>& rThisVariable,
                                                        Vector& rValue)
{
    return rValue;
}

Matrix& LinearIsotropicDamage3DLaw::GetValue(const Variable<Matrix>& rThisVariable,
                                                        Matrix& rValue)
{
    return rValue;
}

array_1d<double, 2>& LinearIsotropicDamage3DLaw::GetValue(
    const Variable<array_1d<double, 2>>& rVariable, array_1d<double, 2>& rValue)
{
    return rValue;
}

array_1d<double, 3>& LinearIsotropicDamage3DLaw::GetValue(
    const Variable<array_1d<double, 3>>& rVariable, array_1d<double, 3>& rValue)
{
    return rValue;
}

// New method - for 3D cases (check definition in the kratos core structure)
// array_1d<double, 6 > & LinearIsotropicDamage3DLaw::GetValue(const
// Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
//{
//    return rValue;
//}

void LinearIsotropicDamage3DLaw::SetValue(const Variable<double>& rVariable,
                                                     const double& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::SetValue(const Variable<Vector>& rVariable,
                                                     const Vector& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::SetValue(const Variable<Matrix>& rVariable,
                                                     const Matrix& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::SetValue(const Variable<array_1d<double, 2>>& rVariable,
                                                     const array_1d<double, 2>& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::SetValue(const Variable<array_1d<double, 3>>& rVariable,
                                                     const array_1d<double, 3>& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

// New Method - for 3D cases (check definition in the kratos core structure)
// void LinearIsotropicDamage3DLaw::SetValue(const
// Variable<array_1d<double, 6 > >& rVariable,
//                                          const array_1d<double, 6 > & rValue,
//                                          const ProcessInfo&
//                                          rCurrentProcessInfo)
//{
//}

bool LinearIsotropicDamage3DLaw::ValidateInput(const Properties& rMaterialProperties)
{
    return true;
}

LinearIsotropicDamage3DLaw::StrainMeasure LinearIsotropicDamage3DLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

LinearIsotropicDamage3DLaw::StressMeasure LinearIsotropicDamage3DLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool LinearIsotropicDamage3DLaw::IsIncremental()
{
    return false;
}

void LinearIsotropicDamage3DLaw::InitializeMaterial(const Properties& material_prop,
                                                               const GeometryType& rElementGeometry,
                                                               const Vector& rShapeFunctionsValues)
{
    r_prev = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    tau_epsilon = 0.;
}

void LinearIsotropicDamage3DLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // update of damage threshold
    r_prev = r;
}

void LinearIsotropicDamage3DLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicDamage3DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& matprops = rValues.GetMaterialProperties();
    Vector& epsilon = rValues.GetStrainVector();
    Vector& sigma_bar = rValues.GetStressVector();
    Vector sigma_bar_pos;
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    double H = matprops[ISOTROPIC_DAMAGE_MODULUS];
    double dpointcoeff;
    double d, q;
    // Uncomment in case of solve ONLY_TRACTION surface for 3D cases
    // double sigma_xx, sigma_yy, sigma_zz, sigma_xz, sigma_yz, sigma_xy;
    // double hyp, sigma_1, sigma_2, sigma_3, angle, cos_a, sin_a;
    bool TRACTION_ONLY = matprops[FLOW_RULE_IS_TRACTION_ONLY];

    // TODO(marcelo): to be removed. Use USE_ELEMENT_PROVIDED_STRAIN flag
    //if (rValues.GetProcessInfo().Has(INITIAL_STRAIN_VECTOR))
    //{
    //    noalias(epsilon) += rValues.GetProcessInfo()[INITIAL_STRAIN_VECTOR];
    //}
    KRATOS_WATCH(epsilon)

    CalculateConstitutiveMatrix(matprops, constitutive_matrix);
    sigma_bar = prod(constitutive_matrix, epsilon);
    sigma_bar_pos = prod(constitutive_matrix, epsilon);
    //TODO: for use with strain energy computation (see below)
    Matrix constitutive_elastic_matrix = constitutive_matrix;

    // for tension-only fluency law:
    // originally sigma and sigma_positive are the same (as it is in the
    // symmetrical case), this block modifies sigma_positive
    if (TRACTION_ONLY)
    {
        // Do not use ONLY TRACTION yield surface (on 3D cases) ...
        // Do something to solve the eigenvalue problem  for 3x3 stress tensor
        // cases

        // Compute the invariants of the stress tensor

        /*
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
        if(sigma_1 > 0){
            sigma_bar_pos(0) += sigma_1 * cos_a * cos_a;
            sigma_bar_pos(1) += sigma_1 * sin_a * sin_a;
            sigma_bar_pos(2) += sigma_1 * sin_a * cos_a;
        }
        if(sigma_2 > 0){
            sigma_bar_pos(0) += sigma_2 * sin_a * sin_a;
            sigma_bar_pos(1) += sigma_2 * cos_a * cos_a;
            sigma_bar_pos(2) -= sigma_2 * sin_a * cos_a;
        }
        */

        KRATOS_THROW_ERROR(
            std::invalid_argument,
            "Isotropic 3D Damage - Yield surface not already implemented", "");
    }

    tau_epsilon = std::sqrt(inner_prod(sigma_bar_pos, epsilon));

    if (tau_epsilon <= r_prev)
    {
        r = r_prev;
        q = CalculateQ(r, matprops);
        d = 1. - q / r;
        constitutive_matrix *= (1 - d);
        sigma_bar *= (1 - d);
    }
    else
    {
        mInelasticFlag = 1;
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

void LinearIsotropicDamage3DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicDamage3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

void LinearIsotropicDamage3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues)
{
}

double LinearIsotropicDamage3DLaw::CalculateQ(double r, const Properties& material_prop)
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

void LinearIsotropicDamage3DLaw::CalculateConstitutiveMatrix(const Properties& props,
                                                                        Matrix& rLinearElasticMatrix)
{
    double YoungModulus = props[YOUNG_MODULUS];
    double PoissonCoefficient = props[POISSON_RATIO];
    // double Ebar = E / (1. - nu * nu);
    // double nubar = nu / (1. - nu);

    rLinearElasticMatrix.clear();

    // D(0, 0) = 1;     D(0, 1) = nubar; D(0, 2) = 0;
    // D(1, 0) = nubar; D(1, 1) = 1;     D(1, 2) = 0;
    // D(2, 0) = 0;     D(2, 1) = 0;     D(2, 2) = 0.5 * (1 - nubar);

    // 3D linear elastic constitutive matrix
    rLinearElasticMatrix(0, 0) =
        (YoungModulus * (1.0 - PoissonCoefficient) /
         ((1.0 + PoissonCoefficient) * (1.0 - 2.0 * PoissonCoefficient)));
    rLinearElasticMatrix(1, 1) = rLinearElasticMatrix(0, 0);
    rLinearElasticMatrix(2, 2) = rLinearElasticMatrix(0, 0);

    rLinearElasticMatrix(3, 3) = rLinearElasticMatrix(0, 0) *
                                 (1.0 - 2.0 * PoissonCoefficient) /
                                 (2.0 * (1.0 - PoissonCoefficient));
    rLinearElasticMatrix(4, 4) = rLinearElasticMatrix(3, 3);
    rLinearElasticMatrix(5, 5) = rLinearElasticMatrix(3, 3);

    rLinearElasticMatrix(0, 1) =
        rLinearElasticMatrix(0, 0) * PoissonCoefficient / (1.0 - PoissonCoefficient);
    rLinearElasticMatrix(1, 0) = rLinearElasticMatrix(0, 1);

    rLinearElasticMatrix(0, 2) = rLinearElasticMatrix(0, 1);
    rLinearElasticMatrix(2, 0) = rLinearElasticMatrix(0, 1);

    rLinearElasticMatrix(1, 2) = rLinearElasticMatrix(0, 1);
    rLinearElasticMatrix(2, 1) = rLinearElasticMatrix(0, 1);

    // D *= Ebar / (1. - nubar * nubar);
}

void LinearIsotropicDamage3DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

int LinearIsotropicDamage3DLaw::Check(const Properties& rMaterialProperties,
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
    // Is mandatory to check the value of the hardening modulus in virtue of the
    // Clauss-Duhem inequality (2nd Thermodynamic law)
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
