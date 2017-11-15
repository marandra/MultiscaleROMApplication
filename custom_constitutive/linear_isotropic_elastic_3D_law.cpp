#include "linear_isotropic_elastic_3D_law.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
// CONSTRUCTOR
LinearIsotropicElastic3DLaw::LinearIsotropicElastic3DLaw()
    : ConstitutiveLaw()
{
}

// COPY CONSTRUCTOR
LinearIsotropicElastic3DLaw::LinearIsotropicElastic3DLaw(const LinearIsotropicElastic3DLaw &rOther)
    : ConstitutiveLaw(rOther)
{
}

// CLONE
ConstitutiveLaw::Pointer LinearIsotropicElastic3DLaw::Clone() const
{
    return ConstitutiveLaw::Pointer(new LinearIsotropicElastic3DLaw());
}

LinearIsotropicElastic3DLaw::SizeType LinearIsotropicElastic3DLaw::WorkingSpaceDimension()
{
    return 3;
}

LinearIsotropicElastic3DLaw::SizeType LinearIsotropicElastic3DLaw::GetStrainSize()
{
    return 6;
}

bool LinearIsotropicElastic3DLaw::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    return false;
}

bool LinearIsotropicElastic3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}

bool LinearIsotropicElastic3DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}

bool LinearIsotropicElastic3DLaw::Has(const Variable<array_1d<double, 2>>& rThisVariable)
{
    return false;
}

bool LinearIsotropicElastic3DLaw::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    return false;
}

// New method -  for 3D cases
// bool LinearIsotropicDamage3DLaw::Has(const
// Variable<array_1d<double, 6 > >& rThisVariable)
//{
//	return false;
//}

double& LinearIsotropicElastic3DLaw::GetValue(const Variable<double>& rThisVariable,
                                                        double& rValue)
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }
    return rValue;
}

Vector& LinearIsotropicElastic3DLaw::GetValue(const Variable<Vector>& rThisVariable,
                                                        Vector& rValue)
{
    return rValue;
}

Matrix& LinearIsotropicElastic3DLaw::GetValue(const Variable<Matrix>& rThisVariable,
                                                        Matrix& rValue)
{
    return rValue;
}

array_1d<double, 2>& LinearIsotropicElastic3DLaw::GetValue(
    const Variable<array_1d<double, 2>>& rVariable, array_1d<double, 2>& rValue)
{
    return rValue;
}

array_1d<double, 3>& LinearIsotropicElastic3DLaw::GetValue(
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

void LinearIsotropicElastic3DLaw::SetValue(const Variable<double>& rVariable,
                                                     const double& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::SetValue(const Variable<Vector>& rVariable,
                                                     const Vector& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::SetValue(const Variable<Matrix>& rVariable,
                                                     const Matrix& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::SetValue(const Variable<array_1d<double, 2>>& rVariable,
                                                     const array_1d<double, 2>& rValue,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::SetValue(const Variable<array_1d<double, 3>>& rVariable,
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

bool LinearIsotropicElastic3DLaw::ValidateInput(const Properties& rMaterialProperties)
{
    return true;
}

LinearIsotropicElastic3DLaw::StrainMeasure LinearIsotropicElastic3DLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

LinearIsotropicElastic3DLaw::StressMeasure LinearIsotropicElastic3DLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool LinearIsotropicElastic3DLaw::IsIncremental()
{
    return false;
}

void LinearIsotropicElastic3DLaw::InitializeMaterial(const Properties& material_prop,
                                                               const GeometryType& rElementGeometry,
                                                               const Vector& rShapeFunctionsValues)
{
}

void LinearIsotropicElastic3DLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void LinearIsotropicElastic3DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& matprops = rValues.GetMaterialProperties();
    Vector& epsilon = rValues.GetStrainVector();
    Vector& sigma_bar = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

    // TODO(marcelo): to be removed. Use USE_ELEMENT_PROVIDED_STRAIN flag
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN))
    {
        noalias(epsilon) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    CalculateConstitutiveMatrix(matprops, constitutive_matrix);
    sigma_bar = prod(constitutive_matrix, epsilon);
    //TODO: for use with strain energy computation (see below)
    Matrix constitutive_elastic_matrix = constitutive_matrix;

    //TODO add check of flag here (COMPUTE_STRAIN_ENERGY)
    mStrainEnergy = 0.5 * (inner_prod(epsilon, prod(constitutive_elastic_matrix, epsilon)));
}

double& LinearIsotropicElastic3DLaw::CalculateValue(Parameters& rParameterValues,
                                                              const Variable<double>& rThisVariable, double& rValue)
{
    //const Properties& MaterialProperties  = rParameterValues.GetMaterialProperties();
    //Vector& StrainVector                  = rParameterValues.GetStrainVector();
    //Vector& StressVector                  = rParameterValues.GetStressVector();
    //const double& E          = MaterialProperties[YOUNG_MODULUS];
    //const double& NU    = MaterialProperties[POISSON_RATIO];

    //if (rThisVariable == STRAIN_ENERGY)
    //{
    //    CalculateCauchyGreenStrain(rParameterValues, StrainVector);
    //    CalculatePK2Stress( StrainVector, StressVector, E, NU );

    //    rValue = 0.5 * inner_prod(StrainVector,StressVector); // Strain energy = 0.5*E:C:E
    //}

    if(rThisVariable == STRAIN_ENERGY){
        rValue = mStrainEnergy;
    }
    return( rValue );
}

void LinearIsotropicElastic3DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearIsotropicElastic3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

void LinearIsotropicElastic3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues)
{
}

double LinearIsotropicElastic3DLaw::CalculateQ(double r, const Properties& material_prop)
{
}

void LinearIsotropicElastic3DLaw::CalculateConstitutiveMatrix(const Properties& props,
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

void LinearIsotropicElastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

int LinearIsotropicElastic3DLaw::Check(const Properties& rMaterialProperties,
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
    return 0;
}
} /* namespace Kratos.*/
