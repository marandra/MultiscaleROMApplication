#include <solid_mechanics_application_variables.h>
#include "linear_elastic_plastic_J2_plane_strain_2D_law.hpp"
#include "multiscale_rom_application_variables.h"
#include "../../../kratos/includes/variables.h"

namespace Kratos
{
    // CONSTRUCTOR
    //ElastoplasticJ2PlanStrain2DLaw
    LinearElasticPlasticJ2PlaneStrain2DLaw::LinearElasticPlasticJ2PlaneStrain2DLaw() 
    	: ConstitutiveLaw()
    	, mInitStrain()
    {
    }
    
    // CLONE
    ConstitutiveLaw::Pointer LinearElasticPlasticJ2PlaneStrain2DLaw::Clone() const
    {
    	return ConstitutiveLaw::Pointer(new LinearElasticPlasticJ2PlaneStrain2DLaw());
    }
    
    LinearElasticPlasticJ2PlaneStrain2DLaw::SizeType LinearElasticPlasticJ2PlaneStrain2DLaw::WorkingSpaceDimension()
    {
    	return 2;
    }
    
    LinearElasticPlasticJ2PlaneStrain2DLaw::SizeType LinearElasticPlasticJ2PlaneStrain2DLaw::GetStrainSize()
    {
    	return 4;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::Has(const Variable<double>& rThisVariable)
    {
    	return false;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::Has(const Variable<Vector>& rThisVariable)
    {
    	return false;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
    	return false;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
    {
    	return false;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
    	return false;
    }
    
    double& LinearElasticPlasticJ2PlaneStrain2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
    	return rValue;
    }
    
    Vector& LinearElasticPlasticJ2PlaneStrain2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
    	if (rThisVariable == INITIAL_STRAIN) {
    		if (rValue.size() != mInitStrain.size())
    		    rValue.resize(mInitStrain.size());
    		noalias(rValue) = mInitStrain;
    	}
    	return (rValue);
    }
    
    Matrix& LinearElasticPlasticJ2PlaneStrain2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
    	return rValue;
    }
    
    array_1d<double, 2 > & LinearElasticPlasticJ2PlaneStrain2DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
    {
    	return rValue;
    }
    
    array_1d<double, 3 > & LinearElasticPlasticJ2PlaneStrain2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
    	return rValue;
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::SetValue(const Variable<double>& rVariable,
    	const double& rValue,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::SetValue(const Variable<Vector >& rVariable,
    	const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    	if (rVariable == INITIAL_STRAIN) {
    		if (rValue.size() == mInitStrain.size())
    			noalias(mInitStrain) = rValue;
    	}
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::SetValue(const Variable<Matrix >& rVariable,
    	const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
    	const array_1d<double, 2 > & rValue,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
    	const array_1d<double, 3 > & rValue,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
	return true;
    }
    
    LinearElasticPlasticJ2PlaneStrain2DLaw::StrainMeasure LinearElasticPlasticJ2PlaneStrain2DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }
    
    LinearElasticPlasticJ2PlaneStrain2DLaw::StressMeasure LinearElasticPlasticJ2PlaneStrain2DLaw::GetStressMeasure()
    {
    	return ConstitutiveLaw::StressMeasure_Cauchy;
    }
    
    bool LinearElasticPlasticJ2PlaneStrain2DLaw::IsIncremental()
    {
    	return false;
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::InitializeMaterial(
    	const Properties& material_prop,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues)
    {
        //mPlasticStrainOld.resize(4);
        mPlasticStrainOld = ZeroVector(this->GetStrainSize());
        mAccumulatedPlasticStrainOld = 0.;
        mInitStrain = ZeroVector(this->GetStrainSize());
    }
	    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
        mPlasticStrainOld = mPlasticStrain;
        mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
    {
    	CalculateMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
    {
    	CalculateMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
    {
    	CalculateMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
    {
    	const Properties& matprops = rValues.GetMaterialProperties();
    	Vector& epsilon = rValues.GetStrainVector();
    	Vector& StressVector = rValues.GetStressVector();
        Matrix ElasticityTensor;
        Matrix& TangentTensor = rValues.GetConstitutiveMatrix(); //TODO find proper getter
        mInitStrain = matprops[INITIAL_STRAIN];
        KRATOS_WATCH(mInitStrain)
        //noalias(epsilon) += mInitStrain;

        double hardening_modulus =  matprops[ISOTROPIC_HARDENING_MODULUS];
        double delta_k =  matprops[INFINITY_HARDENING_MODULUS];
        double hardening_exponent =  matprops[HARDENING_EXPONENT];
        double trial_yield_function;

        double E = matprops[YOUNG_MODULUS];
        double poisson_ratio = matprops[POISSON_RATIO];
        double mu = E / (2. + 2. * poisson_ratio);
        double volumetric_modulus = E / (3 * (1 - 2 * poisson_ratio));

        mPlasticStrain = mPlasticStrainOld;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

        ElasticityTensor.resize(4, 4);
        CalculateElasticityTensor(matprops, ElasticityTensor);
        Vector sigma_trial;
        sigma_trial.resize(4);
        sigma_trial = prod(ElasticityTensor, epsilon - mPlasticStrainOld);

        // StressTrialDev = sigma - 1/3 tr(sigma) * I
        Vector StressTrialDev;
        StressTrialDev.resize(4);
        StressTrialDev = sigma_trial;

        double trace = 0.33333333333333333 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
        StressTrialDev(0) -= trace;
        StressTrialDev(1) -= trace;
        StressTrialDev(2) -= trace;
        double norm_dev_stress = std::sqrt(StressTrialDev(0) * StressTrialDev(0) + StressTrialDev(1) * StressTrialDev(1)
                                           + StressTrialDev(2) * StressTrialDev(2) + 2. * StressTrialDev(3) * StressTrialDev(3));
        trial_yield_function = this->yieldFunction(norm_dev_stress, matprops);
        double dgamma = 0;

        if (trial_yield_function <= 0.){
            StressVector = sigma_trial;
            TangentTensor = ElasticityTensor;
        }
        else{
            Vector YieldFunctionNormalVector = StressTrialDev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0){
                // Exponential softening
                dgamma = GetDeltaGamma (norm_dev_stress, matprops);
            }
            else{
                // Linear softening
                dgamma = trial_yield_function / (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            StressVector(0) = volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) + StressTrialDev(0)
                              - 2. * mu * dgamma * YieldFunctionNormalVector(0);
            StressVector(1) = volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) + StressTrialDev(1)
                              - 2. * mu * dgamma * YieldFunctionNormalVector(1);
            StressVector(2) = volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) + StressTrialDev(2)
                              - 2. * mu * dgamma * YieldFunctionNormalVector(2);
            StressVector(3) = StressTrialDev(3) - 2. * mu * dgamma * YieldFunctionNormalVector(3);

            mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * YieldFunctionNormalVector(0);
            mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * YieldFunctionNormalVector(1);
            mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * YieldFunctionNormalVector(2);
            mPlasticStrain(3) = mPlasticStrainOld(3) + dgamma * YieldFunctionNormalVector(3) * 2;

            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;


            // Actualizar derivada del modulo de hardening-softening

            // Computar tensor tangente
            CalculateTangentTensor(dgamma, norm_dev_stress, YieldFunctionNormalVector, matprops, TangentTensor);

             }

    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
    {
    	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
    {
    	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
    {
    	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
    {
    }
    
    void LinearElasticPlasticJ2PlaneStrain2DLaw::ResetMaterial(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues)
    {
    }

    double LinearElasticPlasticJ2PlaneStrain2DLaw::GetDeltaGamma(double norm_s_trial, const Properties& rMaterialProperties){
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
        while (norm_yieldfunction > tolerance){
            double k_new = yield_stress + (theta*hardening_modulus * mAccumulatedPlasticStrain) + delta_k * (1 - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
            double kp_new = theta * hardening_modulus + delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
            double yieldfunction = -0.8164965809277260 * k_new + norm_s_trial - 2. * mu * dgamma;
            double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
            dgamma = dgamma - yieldfunction / derivative_yieldfunction;
            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;
            norm_yieldfunction = std::abs(yieldfunction);
        }
        //TODO handle the case when no convergence is achieved.
        return dgamma;
    }

    double LinearElasticPlasticJ2PlaneStrain2DLaw::yieldFunction(const double norm_dev_stress, const Properties& rMaterialProperties){
        double yield_stress = rMaterialProperties[YIELD_STRESS];
        double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
        double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
        double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
        double hardening_exponent =  rMaterialProperties[HARDENING_EXPONENT];
        double k_old = yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrainOld)
                       + (delta_k) * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

        return norm_dev_stress - k_old * 0.8164965809277260; //sqrt(2/3)
    }


    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateElasticityTensor(
            const Properties& props, Matrix& D){
        double E = props[YOUNG_MODULUS];
        double poisson_ratio = props[POISSON_RATIO];
        double lambda = E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
        double mu = E / (2. + 2. * poisson_ratio);

        D.clear();

        D(0, 0) = lambda + 2. * mu; D(0, 1) = lambda;           D(0, 2) = lambda;           D(0, 3) = 0.;
        D(1, 0) = lambda;           D(1, 1) = lambda + 2. * mu; D(1, 2) = lambda;           D(1, 3) = 0.;
        D(2, 0) = lambda;           D(2, 1) = lambda;           D(2, 2) = lambda + 2. * mu; D(2, 3) = 0.;
        D(3, 0) = 0.;               D(3, 1) = 0.;               D(3, 2) = 0.;               D(3, 3) = mu;
    }

    void LinearElasticPlasticJ2PlaneStrain2DLaw::CalculateTangentTensor(double dgamma, double norm_s_trial,
                                                                        const Vector& N_new, const Properties& props,
                                                                        Matrix& D){
        double hardening_modulus = props[ISOTROPIC_HARDENING_MODULUS];
        double theta = props[REFERENCE_HARDENING_MODULUS];
        double delta_k = props[INFINITY_HARDENING_MODULUS];
        double hardening_exponent = props[HARDENING_EXPONENT];
    	double E = props[YOUNG_MODULUS];
    	double poisson_ratio = props[POISSON_RATIO];
	    double mu = E / (2. + 2. * poisson_ratio);
        double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

        double kp_new = (theta * hardening_modulus) + delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

        double theta_new   = 1 - (2. * mu * dgamma) / norm_s_trial;
        double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

       D(0, 0) = volumetric_modulus + (2 * mu * theta_new * 2./3.)    - (2 * mu * theta_new_b * (N_new(0) * N_new(0)));
       D(0, 1) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(0) * N_new(1)));
       D(0, 2) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(0) * N_new(2)));
       D(0, 3) =                                                      - (2 * mu * theta_new_b * (N_new(0) * N_new(3)));

       D(1, 0) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(1) * N_new(0)));
       D(1, 1) = volumetric_modulus + (2 * mu * theta_new * 2./3.)    - (2 * mu * theta_new_b * (N_new(1) * N_new(1)));
       D(1, 2) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(1) * N_new(2)));
       D(1, 3) =                                                      - (2 * mu * theta_new_b * (N_new(1) * N_new(3)));

       D(2, 0) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(2) * N_new(0)));
       D(2, 1) = volumetric_modulus + (2 * mu * theta_new * (-1./3.)) - (2 * mu * theta_new_b * (N_new(2) * N_new(1)));
       D(2, 2) = volumetric_modulus + (2 * mu * theta_new * 2./3.)    - (2 * mu * theta_new_b * (N_new(2) * N_new(2)));
       D(2, 3) =                                                      - (2 * mu * theta_new_b * (N_new(2) * N_new(3)));

       D(3, 0) =                                                      - (2 * mu * theta_new_b * (N_new(3) * N_new(0)));
       D(3, 1) =                                                      - (2 * mu * theta_new_b * (N_new(3) * N_new(1)));
       D(3, 2) =                                                      - (2 * mu * theta_new_b * (N_new(3) * N_new(2)));
       D(3, 3) =                                       mu * theta_new - (2 * mu * theta_new_b * (N_new(3) * N_new(3)));
    }

    void LinearElasticPlasticJ2PlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
    {
    	rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    	rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    	rFeatures.mOptions.Set(ISOTROPIC);
    	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    	rFeatures.mStrainSize = GetStrainSize();
    	rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }
    
    int LinearElasticPlasticJ2PlaneStrain2DLaw::Check(
    	const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    		if(!rMaterialProperties.Has(YOUNG_MODULUS)) 
    		    KRATOS_THROW_ERROR(std::invalid_argument, "LinearElasticPlasticJ2PlaneStrain2DLaw - missing YOUNG_MODULUS", "");
    		if(!rMaterialProperties.Has(POISSON_RATIO)) 
    		    KRATOS_THROW_ERROR(std::invalid_argument, "LinearElasticPlasticJ2PlaneStrain2DLaw - missing POISSON_RATIO", "");
    		return 0;
    }



} /* namespace Kratos.*/
