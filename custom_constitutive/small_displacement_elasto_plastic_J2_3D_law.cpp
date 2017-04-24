#include "small_displacement_elasto_plastic_J2_3D_law.hpp"
#include "multiscale_rom_application_variables.h"
#include <solid_mechanics_application_variables.h>

namespace Kratos {
// CONSTRUCTOR
// ElastoplasticJ2PlanStrain2DLaw
SmallDisplacementElastoPlasticJ23DLaw::SmallDisplacementElastoPlasticJ23DLaw()
    : ConstitutiveLaw() {}

// CLONE
ConstitutiveLaw::Pointer SmallDisplacementElastoPlasticJ23DLaw::Clone() const {
  return ConstitutiveLaw::Pointer(new SmallDisplacementElastoPlasticJ23DLaw());
}

SmallDisplacementElastoPlasticJ23DLaw::SizeType
SmallDisplacementElastoPlasticJ23DLaw::WorkingSpaceDimension() {
  return 2;
}

SmallDisplacementElastoPlasticJ23DLaw::SizeType
SmallDisplacementElastoPlasticJ23DLaw::GetStrainSize() {
  return 6;
}

bool SmallDisplacementElastoPlasticJ23DLaw::Has(
    const Variable<double> &rThisVariable) {
  return false;
}

bool SmallDisplacementElastoPlasticJ23DLaw::Has(
    const Variable<Vector> &rThisVariable) {
  return false;
}

bool SmallDisplacementElastoPlasticJ23DLaw::Has(
    const Variable<Matrix> &rThisVariable) {
  return false;
}

bool SmallDisplacementElastoPlasticJ23DLaw::Has(
    const Variable<array_1d<double, 2>> &rThisVariable) {
  return false;
}

bool SmallDisplacementElastoPlasticJ23DLaw::Has(
    const Variable<array_1d<double, 3>> &rThisVariable) {
  return false;
}

double &SmallDisplacementElastoPlasticJ23DLaw::GetValue(
    const Variable<double> &rThisVariable, double &rValue) {
  return rValue;
}

Vector &SmallDisplacementElastoPlasticJ23DLaw::GetValue(
    const Variable<Vector> &rThisVariable, Vector &rValue) {
  return rValue;
}

Matrix &SmallDisplacementElastoPlasticJ23DLaw::GetValue(
    const Variable<Matrix> &rThisVariable, Matrix &rValue) {
  return rValue;
}

array_1d<double, 2> &SmallDisplacementElastoPlasticJ23DLaw::GetValue(
    const Variable<array_1d<double, 2>> &rVariable,
    array_1d<double, 2> &rValue) {
  return rValue;
}

array_1d<double, 3> &SmallDisplacementElastoPlasticJ23DLaw::GetValue(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rValue) {
  return rValue;
}

void SmallDisplacementElastoPlasticJ23DLaw::SetValue(
    const Variable<double> &rVariable, const double &rValue,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::SetValue(
    const Variable<Vector> &rVariable, const Vector &rValue,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::SetValue(
    const Variable<Matrix> &rVariable, const Matrix &rValue,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::SetValue(
    const Variable<array_1d<double, 2>> &rVariable,
    const array_1d<double, 2> &rValue, const ProcessInfo &rCurrentProcessInfo) {
}

void SmallDisplacementElastoPlasticJ23DLaw::SetValue(
    const Variable<array_1d<double, 3>> &rVariable,
    const array_1d<double, 3> &rValue, const ProcessInfo &rCurrentProcessInfo) {
}

bool SmallDisplacementElastoPlasticJ23DLaw::ValidateInput(
    const Properties &rMaterialProperties) {
  return true;
}

SmallDisplacementElastoPlasticJ23DLaw::StrainMeasure
SmallDisplacementElastoPlasticJ23DLaw::GetStrainMeasure() {
  return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

SmallDisplacementElastoPlasticJ23DLaw::StressMeasure
SmallDisplacementElastoPlasticJ23DLaw::GetStressMeasure() {
  return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool SmallDisplacementElastoPlasticJ23DLaw::IsIncremental() { return false; }

void SmallDisplacementElastoPlasticJ23DLaw::InitializeMaterial(
    const Properties &material_prop, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues) {
  // mPlasticStrainOld.resize(4);
  mPlasticStrainOld = ZeroVector(this->GetStrainSize());
  mAccumulatedPlasticStrainOld = 0.0;
}

void SmallDisplacementElastoPlasticJ23DLaw::InitializeSolutionStep(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeSolutionStep(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo) {
  mPlasticStrainOld = mPlasticStrain;
  mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
}

void SmallDisplacementElastoPlasticJ23DLaw::InitializeNonLinearIteration(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeNonLinearIteration(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo) {}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateMaterialResponsePK1(
    Parameters &rValues) {
  CalculateMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateMaterialResponsePK2(
    Parameters &rValues) {
  CalculateMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateMaterialResponseKirchhoff(
    Parameters &rValues) {
  CalculateMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateMaterialResponseCauchy(
    Parameters &rValues) {
  const Properties &matprops = rValues.GetMaterialProperties();
  Vector &epsilon = rValues.GetStrainVector();
  Vector &StressVector = rValues.GetStressVector();
  Matrix ElasticityTensor;
  Matrix &TangentTensor =
      rValues.GetConstitutiveMatrix(); // TODO find proper getter
  double hardening_modulus = matprops[ISOTROPIC_HARDENING_MODULUS];
  double delta_k = matprops[INFINITY_HARDENING_MODULUS];
  double hardening_exponent = matprops[HARDENING_EXPONENT];
  double trial_yield_function;

  double E = matprops[YOUNG_MODULUS];
  double poisson_ratio = matprops[POISSON_RATIO];
  double mu = E / (2. + 2. * poisson_ratio);
  double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

  if (rValues.GetProcessInfo().Has(INITIAL_STRAIN_VECTOR)) {
    noalias(epsilon) += rValues.GetProcessInfo()[INITIAL_STRAIN_VECTOR];
  }

  mPlasticStrain = mPlasticStrainOld;
  mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

  ElasticityTensor.resize(6, 6);
  CalculateElasticityTensor(matprops, ElasticityTensor);
  Vector sigma_trial;
  sigma_trial.resize(6);
  sigma_trial = prod(ElasticityTensor, epsilon - mPlasticStrainOld);

  // StressTrialDev = sigma - 1/3 tr(sigma) * I
  Vector StressTrialDev;
  StressTrialDev.resize(6);
  StressTrialDev = sigma_trial;

  double trace =
      0.33333333333333333 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
  StressTrialDev(0) -= trace;
  StressTrialDev(1) -= trace;
  StressTrialDev(2) -= trace;
  double norm_dev_stress =
      std::sqrt(StressTrialDev(0) * StressTrialDev(0) +
                StressTrialDev(1) * StressTrialDev(1) +
                StressTrialDev(2) * StressTrialDev(2) +
                2. * StressTrialDev(3) * StressTrialDev(3) +
                2. * StressTrialDev(4) * StressTrialDev(4) +
                2. * StressTrialDev(5) * StressTrialDev(5));
  trial_yield_function = this->yieldFunction(norm_dev_stress, matprops);
  double dgamma = 0;

  if (trial_yield_function <= 0.) {
    StressVector = sigma_trial;
    TangentTensor = ElasticityTensor;
  } else {
    Vector YieldFunctionNormalVector = StressTrialDev / norm_dev_stress;
    if (delta_k != 0.0 && hardening_exponent != 0.0) {
      // Exponential softening
      dgamma = GetDeltaGamma(norm_dev_stress, matprops);
    } else {
      // Linear softening
      dgamma = trial_yield_function /
               (2. * mu * (1. + (hardening_modulus / (3. * mu))));
    }

    StressVector(0) =
        volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
        StressTrialDev(0) - 2. * mu * dgamma * YieldFunctionNormalVector(0);
    StressVector(1) =
        volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
        StressTrialDev(1) - 2. * mu * dgamma * YieldFunctionNormalVector(1);
    StressVector(2) =
        volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
        StressTrialDev(2) - 2. * mu * dgamma * YieldFunctionNormalVector(2);
    StressVector(3) =
        StressTrialDev(3) - 2. * mu * dgamma * YieldFunctionNormalVector(3);
    StressVector(4) =
        StressTrialDev(4) - 2. * mu * dgamma * YieldFunctionNormalVector(4);
    StressVector(5) =
        StressTrialDev(5) - 2. * mu * dgamma * YieldFunctionNormalVector(5);

    mPlasticStrain(0) =
        mPlasticStrainOld(0) + dgamma * YieldFunctionNormalVector(0);
    mPlasticStrain(1) =
        mPlasticStrainOld(1) + dgamma * YieldFunctionNormalVector(1);
    mPlasticStrain(2) =
        mPlasticStrainOld(2) + dgamma * YieldFunctionNormalVector(2);
    mPlasticStrain(3) =
        mPlasticStrainOld(3) + dgamma * YieldFunctionNormalVector(3) * 2;
    mPlasticStrain(4) =
        mPlasticStrainOld(4) + dgamma * YieldFunctionNormalVector(4) * 2;
    mPlasticStrain(5) =
        mPlasticStrainOld(5) + dgamma * YieldFunctionNormalVector(5) * 2;

    mAccumulatedPlasticStrain =
        mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;

    // Actualizar derivada del modulo de hardening-softening

    // Computar tensor tangente
    CalculateTangentTensor(dgamma, norm_dev_stress, YieldFunctionNormalVector,
                           matprops, TangentTensor);
  }
}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeMaterialResponsePK1(
    Parameters &rValues) {
  FinalizeMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeMaterialResponsePK2(
    Parameters &rValues) {
  FinalizeMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeMaterialResponseKirchhoff(
    Parameters &rValues) {
  FinalizeMaterialResponseCauchy(rValues);
}

void SmallDisplacementElastoPlasticJ23DLaw::FinalizeMaterialResponseCauchy(
    Parameters &rValues) {}

void SmallDisplacementElastoPlasticJ23DLaw::ResetMaterial(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const Vector &rShapeFunctionsValues) {}

double SmallDisplacementElastoPlasticJ23DLaw::GetDeltaGamma(
    double norm_s_trial, const Properties &rMaterialProperties) {
  double E = rMaterialProperties[YOUNG_MODULUS];
  double poisson_ratio = rMaterialProperties[POISSON_RATIO];
  double yield_stress = rMaterialProperties[YIELD_STRESS];
  double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
  double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
  double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
  double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
  double tolerance = 1e-6 * yield_stress;
  double mu = E / (2. * (1. + poisson_ratio));
  double dgamma = 0.0;
  double norm_yieldfunction = 1.0;
  mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;
  while (norm_yieldfunction > tolerance) {
    double k_new =
        yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrain) +
        delta_k *
            (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    double kp_new =
        theta * hardening_modulus +
        delta_k * (hardening_exponent *
                   std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    double yieldfunction =
        -0.8164965809277260 * k_new + norm_s_trial - 2. * mu * dgamma;
    double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
    dgamma = dgamma - yieldfunction / derivative_yieldfunction;
    mAccumulatedPlasticStrain =
        mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;
    norm_yieldfunction = std::abs(yieldfunction);
  }
  // TODO handle the case when no convergence is achieved.
  return dgamma;
}

double SmallDisplacementElastoPlasticJ23DLaw::yieldFunction(
    const double norm_dev_stress, const Properties &rMaterialProperties) {
  double yield_stress = rMaterialProperties[YIELD_STRESS];
  double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
  double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
  double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
  double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
  double k_old =
      yield_stress +
      (theta * hardening_modulus * mAccumulatedPlasticStrainOld) +
      (delta_k) *
          (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

  return norm_dev_stress - k_old * 0.8164965809277260; // sqrt(2/3)
}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateElasticityTensor(
    const Properties &props, Matrix &D) {
  double E = props[YOUNG_MODULUS];
  double poisson_ratio = props[POISSON_RATIO];
  double lambda =
      E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
  double mu = E / (2. + 2. * poisson_ratio);

  D.clear();
  D.resize(6, 6, false);
  D = ZeroMatrix(6, 6);

  D(0, 0) = lambda + 2. * mu;
  D(0, 1) = lambda;
  D(0, 2) = lambda;
  D(1, 0) = lambda;
  D(1, 1) = lambda + 2. * mu;
  D(1, 2) = lambda;
  D(2, 0) = lambda;
  D(2, 1) = lambda;
  D(2, 2) = lambda + 2. * mu;
  D(3, 3) = mu;
  D(4, 4) = mu;
  D(5, 5) = mu;
}

void SmallDisplacementElastoPlasticJ23DLaw::CalculateTangentTensor(
    double dgamma, double norm_s_trial, const Vector &N_new,
    const Properties &props, Matrix &D) {
  double hardening_modulus = props[ISOTROPIC_HARDENING_MODULUS];
  double theta = props[REFERENCE_HARDENING_MODULUS];
  double delta_k = props[INFINITY_HARDENING_MODULUS];
  double hardening_exponent = props[HARDENING_EXPONENT];
  double E = props[YOUNG_MODULUS];
  double poisson_ratio = props[POISSON_RATIO];
  double mu = E / (2. + 2. * poisson_ratio);
  double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

  double kp_new =
      (theta * hardening_modulus) +
      delta_k * (hardening_exponent *
                 std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

  double theta_new = 1 - (2. * mu * dgamma) / norm_s_trial;
  double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

  D(0, 0) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
            (2 * mu * theta_new_b * (N_new(0) * N_new(0)));
  D(0, 1) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(0) * N_new(1)));
  D(0, 2) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(0) * N_new(2)));
  D(0, 3) = -(2 * mu * theta_new_b * (N_new(0) * N_new(3)));
  D(0, 4) = -(2 * mu * theta_new_b * (N_new(0) * N_new(4)));
  D(0, 5) = -(2 * mu * theta_new_b * (N_new(0) * N_new(5)));

  D(1, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(1) * N_new(0)));
  D(1, 1) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
            (2 * mu * theta_new_b * (N_new(1) * N_new(1)));
  D(1, 2) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(1) * N_new(2)));
  D(1, 3) = -(2 * mu * theta_new_b * (N_new(1) * N_new(3)));
  D(1, 4) = -(2 * mu * theta_new_b * (N_new(1) * N_new(4)));
  D(1, 5) = -(2 * mu * theta_new_b * (N_new(1) * N_new(5)));

  D(2, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(2) * N_new(0)));
  D(2, 1) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
            (2 * mu * theta_new_b * (N_new(2) * N_new(1)));
  D(2, 2) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
            (2 * mu * theta_new_b * (N_new(2) * N_new(2)));
  D(2, 3) = -(2 * mu * theta_new_b * (N_new(2) * N_new(3)));
  D(2, 4) = -(2 * mu * theta_new_b * (N_new(2) * N_new(4)));
  D(2, 5) = -(2 * mu * theta_new_b * (N_new(2) * N_new(5)));

  D(3, 0) = -(2 * mu * theta_new_b * (N_new(3) * N_new(0)));
  D(3, 1) = -(2 * mu * theta_new_b * (N_new(3) * N_new(1)));
  D(3, 2) = -(2 * mu * theta_new_b * (N_new(3) * N_new(2)));
  D(3, 3) = mu * theta_new - (2 * mu * theta_new_b * (N_new(3) * N_new(3)));
  D(3, 4) = -(2 * mu * theta_new_b * (N_new(3) * N_new(4)));
  D(3, 5) = -(2 * mu * theta_new_b * (N_new(3) * N_new(5)));

  D(4, 0) = -(2 * mu * theta_new_b * (N_new(4) * N_new(0)));
  D(4, 1) = -(2 * mu * theta_new_b * (N_new(4) * N_new(1)));
  D(4, 2) = -(2 * mu * theta_new_b * (N_new(4) * N_new(2)));
  D(4, 3) = -(2 * mu * theta_new_b * (N_new(4) * N_new(3)));
  D(4, 4) = mu * theta_new - (2 * mu * theta_new_b * (N_new(4) * N_new(4)));
  D(4, 5) = -(2 * mu * theta_new_b * (N_new(4) * N_new(5)));

  D(5, 0) = -(2 * mu * theta_new_b * (N_new(5) * N_new(0)));
  D(5, 1) = -(2 * mu * theta_new_b * (N_new(5) * N_new(1)));
  D(5, 2) = -(2 * mu * theta_new_b * (N_new(5) * N_new(2)));
  D(5, 3) = -(2 * mu * theta_new_b * (N_new(5) * N_new(3)));
  D(5, 4) = -(2 * mu * theta_new_b * (N_new(5) * N_new(4)));
  D(5, 5) = mu * theta_new - (2 * mu * theta_new_b * (N_new(5) * N_new(5)));
}

void SmallDisplacementElastoPlasticJ23DLaw::GetLawFeatures(
    Features &rFeatures) {
  rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
  rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
  rFeatures.mOptions.Set(ISOTROPIC);
  rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
  rFeatures.mStrainSize = GetStrainSize();
  rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

int SmallDisplacementElastoPlasticJ23DLaw::Check(
    const Properties &rMaterialProperties, const GeometryType &rElementGeometry,
    const ProcessInfo &rCurrentProcessInfo) {
  if (!rMaterialProperties.Has(YOUNG_MODULUS))
    KRATOS_THROW_ERROR(
        std::invalid_argument,
        "SmallDisplacementsElastoPlasticJ23DLaw - missing YOUNG_MODULUS", "");
  if (!rMaterialProperties.Has(POISSON_RATIO))
    KRATOS_THROW_ERROR(
        std::invalid_argument,
        "SmallDisplacementsElastoPlasticJ23DLaw - missing POISSON_RATIO", "");
  return 0;
}
} /* namespace Kratos.*/
