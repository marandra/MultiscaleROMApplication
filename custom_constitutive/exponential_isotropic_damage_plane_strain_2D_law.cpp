#include "exponential_isotropic_damage_plane_strain_2D_law.hpp"
#include "multiscale_rom_application_variables.h"
//#include "../SolidMechanicsApplication/solid_mechanics_application_variables.h"

namespace Kratos
{
    // CONSTRUCTOR
    ExponentialIsotropicDamagePlaneStrain2DLaw::ExponentialIsotropicDamagePlaneStrain2DLaw() 
    	: ConstitutiveLaw()
    	//, m_initialized(false)
    	//, m_init_gradT()
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
    
    bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
    {
    	return false;
    }
    
    bool ExponentialIsotropicDamagePlaneStrain2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
    	return false;
    }
    
    double& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
    	return rValue;
    }
    
    Vector& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
    	//std::stringstream ss;
    	//ss << "LinearIsotropicDamagePlaneStrain2DLaw::GetValue" << std::endl;
    	//if (rThisVariable == INITIAL_TEMP_GRAD) {
    	//	if (rValue.size() != m_init_gradT.size())
    	//		rValue.resize(m_init_gradT.size());
    	//	noalias(rValue) = m_init_gradT;
    	//}
    	//if (rThisVariable == FLUX_RVE || rThisVariable == HEAT_FLUX_RVE) {
    	//	if (rValue.size() != mStressVector.size())
    	//		rValue.resize(mStressVector.size());
    	//	noalias(rValue) = mStressVector;
    	//}
    	//if (rThisVariable == HEAT_FLUX_RVE) { //For Output
    	//	if (rValue.size() != 6)
    	//		rValue.resize(6);
    	//	rValue(0) = mStressVector(0); // / 1.0e6; //[W/mm^2]
    	//	rValue(1) = mStressVector(1); // / 1.0e6;
    	//	rValue(2) = mStressVector(2); // / 1.0e6;
    	//	rValue(3) = 0.0;
    	//	rValue(4) = 0.0;
    	//	rValue(5) = 0.0;
    
    	//	//ss << "HEAT_FLUX_RVE = " << rValue << ", " << std::endl;
    	//	//std::cout << ss.str();
    	//}
    	return rValue;
    }
    
    Matrix& ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
    	return rValue;
    }
    
    array_1d<double, 2 > & ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
    {
    	return rValue;
    }
    
    array_1d<double, 3 > & ExponentialIsotropicDamagePlaneStrain2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
    	return rValue;
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<double>& rVariable,
    	const double& rValue,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Vector >& rVariable,
    	const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    	//if (rVariable == INITIAL_TEMP_GRAD) {
    	//	if (rValue.size() == m_init_gradT.size())
    	//		noalias(m_init_gradT) = rValue;
    	//}
    	//if (rVariable == FLUX_RVE)
    	//{
    	//	if (rValue.size() == mStressVector.size())
    	//		noalias(mStressVector) = rValue;
    	//}
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<Matrix >& rVariable,
    	const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
    	const array_1d<double, 2 > & rValue,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
    	const array_1d<double, 3 > & rValue,
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
    	const Properties& material_prop,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues)
    {
	double nu = material_prop[POISSON_RATIO];
	r_prev = std::sqrt(1 - nu * nu) * material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
	tau_e = 0.;
    }
	    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
        // update of damage threshold
	r_prev = r;
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
    {
    //	CalculateMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
    {
    //	CalculateMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
    {
    //	CalculateMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
    {
    	// get some references
    	const Properties& matprops = rValues.GetMaterialProperties();
    	Vector& strain = rValues.GetStrainVector();
    	Vector& stress = rValues.GetStressVector();
    	Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
	double H = matprops[ISOTROPIC_HARDENING_MODULUS];
	double dpointcoeff;
	double d, q;

            // sigma_bar = C : epsilon
    	    CalculateConstitutiveMatrix(matprops, constitutiveMatrix);
            stress = prod(constitutiveMatrix, strain);
            // tau_epsilon = sqrt(epsilon : sigma_bar)
            tau_e = std::sqrt(inner_prod(strain, stress));

	    // r = r_prev
	    // d = 1 - q(r) / r
	    // sigma = (1 - d) * sigma_bar
	    // C_tan = (1 - d) * C
            if (tau_e <= r_prev) {
                r = r_prev;
                q = CalculateQ(r, matprops);
                d = 1. - q / r;
                stress *= (1 - d); 
                constitutiveMatrix *= (1 - d); 
	    }
	    // r = tau_e
	    // d = 1 - q(r) / r
	    // sigma = (1 - d) * sigma_bar
	    // C_tan = (1 - d) * C - q(r)-H/r3 * sigma x sigma
	    else {
                r = tau_e;
                q = CalculateQ(r, matprops);
                d = 1. - q / r;
                stress *= (1. - d); 
                dpointcoeff = (q - H * r)/(r * r * r);
                constitutiveMatrix *= (1. - d); 
                constitutiveMatrix -= dpointcoeff * outer_prod(stress, stress); 
	    }

	    //std::cout << "DEBUG strain " << strain<< std::endl;
	    //std::cout << "DEBUG stress " << stress<< std::endl;
	    //std::cout << "DEBUG tau            " << tau_e << std::endl;
	    //std::cout << "DEBUG r " << r_prev << std::endl;
	    //std::cout << "DEBUG C_sec " << constitutiveMatrix << std::endl;
	    //std::cout << "DEBUG curve " << tau_e << " " << r << " " << q << " " << std::endl;

    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
    {
    //	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
    {
    //	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
    {
    //	FinalizeMaterialResponseCauchy(rValues);
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
    {
    }
    
    void ExponentialIsotropicDamagePlaneStrain2DLaw::ResetMaterial(const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const Vector& rShapeFunctionsValues)
    {
    }
    
    double ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateQ(double r,
    	const Properties& material_prop) {

	double H = material_prop[ISOTROPIC_HARDENING_MODULUS];
	double nu = material_prop[POISSON_RATIO];
	double r0 = std::sqrt(1 - nu * nu) * material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    	double q_inf = std::sqrt(1 - nu * nu) * material_prop[INFINITY_YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
        double q;

	if (r < r0)
	    return r;
	q = q_inf - (q_inf - r0) * std::exp(H * (1 - r / r0));
	return q;
    }

    void ExponentialIsotropicDamagePlaneStrain2DLaw::CalculateConstitutiveMatrix(
        const Properties& props, Matrix& D)
    {
    	double E = props[YOUNG_MODULUS];
    	double nu = props[POISSON_RATIO];
	double Ebar = E / (1. - nu * nu);
	double nubar = nu / (1. - nu);

    	D.clear();
    
        D(0, 0) = 1;     D(0, 1) = nubar; D(0, 2) = 0; 
        D(1, 0) = nubar; D(1, 1) = 1;     D(1, 2) = 0; 
        D(2, 0) = 0;     D(2, 1) = 0;     D(2, 2) = 0.5 * (1 - nubar); 

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
    
    int ExponentialIsotropicDamagePlaneStrain2DLaw::Check(
    	const Properties& rMaterialProperties,
    	const GeometryType& rElementGeometry,
    	const ProcessInfo& rCurrentProcessInfo)
    {
    		if(!rMaterialProperties.Has(YOUNG_MODULUS)) 
    		    KRATOS_THROW_ERROR(std::invalid_argument, "ExponentialIsotropicDamagePlaneStrain2DLaw - missing YOUNG_MODULUS", "");
    		if(!rMaterialProperties.Has(POISSON_RATIO)) 
    		    KRATOS_THROW_ERROR(std::invalid_argument, "ExponentialIsotropicDamagePlaneStrain2DLaw - missing POISSON_RATIO", "");
    		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS)) 
    		KRATOS_THROW_ERROR(std::invalid_argument, "ExponentialIsotropicDamagePlaneStrain2DLaw - missing ISOTROPIC_HARDENING_MODULUS", "");
    		if(!rMaterialProperties.Has(INFINITY_YIELD_STRESS)) 
    		KRATOS_THROW_ERROR(std::invalid_argument, "ExponentialIsotropicDamagePlaneStrain2DLaw - missing INFINITY_YIELD_STRESS", "");
    		if(rMaterialProperties[INFINITY_YIELD_STRESS] < 0) 
    		    KRATOS_THROW_ERROR(std::invalid_argument, "ExponentialIsotropicDamagePlaneStrain2DLaw - INFINITY_YIELD_STRESS must be positive", "");
    		return 0;
    }
    
} /* namespace Kratos.*/
