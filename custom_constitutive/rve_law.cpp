#include "rve_law.h"

namespace Kratos
{
// CLONE
// ConstitutiveLaw::Pointer RVELaw::Clone() const
//{
//	return ConstitutiveLaw::Pointer(new RVELaw());
//}
RVELaw::RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param)
    : mpRVEModelPart(mpModelPart)
{
    auto w_list = param["w"];
    auto B_list = param["B"];
    auto prop_id = param["props_id"];
    unsigned int nr_points = B_list.size();
    unsigned int nr_comps = B_list[0].size();
    unsigned int nr_modes = B_list[0][0].size();

    for (unsigned int i = 0; i < nr_points; i++)
    {
        Matrix BK(nr_comps, nr_modes);
        for (unsigned int c = 0; c < nr_comps; c++)
        {
            for (unsigned int m = 0; m < nr_modes; m++)
            {
                mB_list.push_back(BK);
                BK(c, m) = B_list[i][c][m].GetDouble();
            }
        }
        KRATOS_WATCH(BK)
        mB_list.push_back(BK);

        mIW_list.push_back(w_list[i].GetDouble());

        // NOTE: here the rMaterialProperties come from the COARSE SCALE
        // we are assuming however it also contains the materials to be used in
        // the small scale
        auto prop = mpRVEModelPart->pGetProperties(prop_id[i].GetInt());
        mprop_list.push_back(prop);

        ConstitutiveLaw::Pointer porigin_cl = prop->GetValue(CONSTITUTIVE_LAW);
        ConstitutiveLaw::Pointer pcl = porigin_cl->Clone();
        KRATOS_WATCH(*pcl)
        mCL_list.push_back(pcl);
    }
}

void RVELaw::InitializeMaterial(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const Vector& rShapeFunctionsValues)
{
    KRATOS_WATCH("inside initialize material")
}


/*
void LinearIsotropicDamage3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
  const Properties& matprops = rValues.GetMaterialProperties();
  Vector& epsilon = rValues.GetStrainVector();
  Vector& sigma_bar = rValues.GetStressVector();
  Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();




}
*/
// CL functions
// RVELaw::SizeType HomogenizedRVEResponse2D::WorkingSpaceDimension()
// {
// 	return 2;
// }
//
// RVELaw::SizeType HomogenizedRVEResponse2D::GetStrainSize()
// {
// 	return 3;
// }
//
// bool RVELaw::Has(const Variable<double>& rThisVariable)
// {
// 	return false;
// }
//
// bool RVELaw::Has(const Variable<Vector>& rThisVariable)
// {
// 	return false;
// }
//
// bool RVELaw::Has(const Variable<Matrix>& rThisVariable)
// {
// 	return false;
// }
//
// bool RVELaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
// {
// 	return false;
// }
//
// bool RVELaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
// {
// 	return false;
// }
//
// double& RVELaw::GetValue(const Variable<double>& rThisVariable, double&
// rValue)
// {
// 	return rValue;
// }
//
// Vector& RVELaw::GetValue(const Variable<Vector>& rThisVariable, Vector&
// rValue)
// {
// 	return rValue;
// }
//
// Matrix& RVELaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix&
// rValue)
// {
// 	return rValue;
// }
//
// array_1d<double, 2 > & RVELaw::GetValue(const Variable<array_1d<double, 2 >
// >& rVariable, array_1d<double, 2 > & rValue)
// {
// 	return rValue;
// }
//
// array_1d<double, 3 > & RVELaw::GetValue(const Variable<array_1d<double, 3 >
// >& rVariable, array_1d<double, 3 > & rValue)
// {
// 	return rValue;
// }
//
// void RVELaw::SetValue(const Variable<double>& rVariable,
// 	const double& rValue,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::SetValue(const Variable<Vector >& rVariable,
// 	const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::SetValue(const Variable<Matrix >& rVariable,
// 	const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
// 	const array_1d<double, 2 > & rValue,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
// 	const array_1d<double, 3 > & rValue,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// bool RVELaw::ValidateInput(const Properties& rMaterialProperties)
// {
//	return true;
// }
//
// RVELaw::StrainMeasure HomogenizedRVEResponse2D::GetStrainMeasure()
// {
// 	return ConstitutiveLaw::StrainMeasure_Infinitesimal;
// }
//
// RVELaw::StressMeasure HomogenizedRVEResponse2D::GetStressMeasure()
// {
// 	return ConstitutiveLaw::StressMeasure_Cauchy;
// }
//
// bool RVELaw::IsIncremental()
// {
// 	return false;
// }
//
// void RVELaw::InitializeMaterial(
// 	const Properties& material_prop,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues)
// {
//	double nu = material_prop[POISSON_RATIO];
//	r_prev = std::sqrt(1 - nu * nu) * material_prop[YIELD_STRESS] /
// std::sqrt(material_prop[YOUNG_MODULUS]);
//	tau_e = 0.;
// }
//
// void RVELaw::InitializeSolutionStep(const Properties& rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
//     // update of damage threshold
//	r_prev = r;
// }
//
// void RVELaw::InitializeNonLinearIteration(const Properties&
// rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::FinalizeNonLinearIteration(const Properties&
// rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// }
//
// void RVELaw::CalculateMaterialResponsePK1(Parameters& rValues)
// {
// //	CalculateMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::CalculateMaterialResponsePK2(Parameters& rValues)
// {
// //	CalculateMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
// {
// //	CalculateMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::CalculateMaterialResponseCauchy(Parameters& rValues)
// {
// 	// get some references
// 	const Properties& matprops = rValues.GetMaterialProperties();
// 	Vector& strain = rValues.GetStrainVector();
// 	Vector& stress = rValues.GetStressVector();
// 	Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
//	double H = matprops[ISOTROPIC_DAMAGE_MODULUS];
//	double dpointcoeff;
//	double d, q;

//         // sigma_bar = C : epsilon
// 	    CalculateConstitutiveMatrix(matprops, constitutiveMatrix);
//         stress = prod(constitutiveMatrix, strain);
//         // tau_epsilon = sqrt(epsilon : sigma_bar)
//         tau_e = std::sqrt(inner_prod(strain, stress));

//	    // r = r_prev
//	    // d = 1 - q(r) / r
//	    // sigma = (1 - d) * sigma_bar
//	    // C_tan = (1 - d) * C
//         if (tau_e <= r_prev) {
//             r = r_prev;
//             q = CalculateQ(r, matprops);
//             d = 1. - q / r;
//             stress *= (1 - d);
//             constitutiveMatrix *= (1 - d);
//	    }
//	    // r = tau_e
//	    // d = 1 - q(r) / r
//	    // sigma = (1 - d) * sigma_bar
//	    // C_tan = (1 - d) * C - q(r)-H/r3 * sigma x sigma
//	    else {
//             r = tau_e;
//             q = CalculateQ(r, matprops);
//             d = 1. - q / r;
//             stress *= (1. - d);
//             dpointcoeff = (q - H * r)/(r * r * r);
//             constitutiveMatrix *= (1. - d);
//             constitutiveMatrix -= dpointcoeff * outer_prod(stress, stress);
//	    }

//	    //std::cout << "DEBUG strain " << strain<< std::endl;
//	    //std::cout << "DEBUG stress " << stress<< std::endl;
//	    //std::cout << "DEBUG tau            " << tau_e << std::endl;
//	    //std::cout << "DEBUG r " << r_prev << std::endl;
//	    //std::cout << "DEBUG C_sec " << constitutiveMatrix << std::endl;
//	    //std::cout << "DEBUG curve " << tau_e << " " << r << " " << q << " " <<
// std::endl;

// }
//
// void RVELaw::FinalizeMaterialResponsePK1(Parameters& rValues)
// {
// //	FinalizeMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::FinalizeMaterialResponsePK2(Parameters& rValues)
// {
// //	FinalizeMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
// {
// //	FinalizeMaterialResponseCauchy(rValues);
// }
//
// void RVELaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
// {
// }
//
// void RVELaw::ResetMaterial(const Properties& rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const Vector& rShapeFunctionsValues)
// {
// }
//
// double RVELaw::CalculateQ(double r,
// 	const Properties& material_prop) {

//	double H = material_prop[ISOTROPIC_DAMAGE_MODULUS];
//	double nu = material_prop[POISSON_RATIO];
//	double r0 = std::sqrt(1 - nu * nu) * material_prop[YIELD_STRESS] /
// std::sqrt(material_prop[YOUNG_MODULUS]);
// 	double q_inf = std::sqrt(1 - nu * nu) * material_prop[INFINITY_YIELD_STRESS]
// / std::sqrt(material_prop[YOUNG_MODULUS]);
//     double q;

//	if (r < r0)
//	    return r;
//	q = q_inf - (q_inf - r0) * std::exp(H * (1 - r / r0));
//	return q;
// }

// void RVELaw::CalculateConstitutiveMatrix(
//     const Properties& props, Matrix& D)
// {
// 	double E = props[YOUNG_MODULUS];
// 	double nu = props[POISSON_RATIO];
//	double Ebar = E / (1. - nu * nu);
//	double nubar = nu / (1. - nu);

// 	D.clear();
//
//     D(0, 0) = 1;     D(0, 1) = nubar; D(0, 2) = 0;
//     D(1, 0) = nubar; D(1, 1) = 1;     D(1, 2) = 0;
//     D(2, 0) = 0;     D(2, 1) = 0;     D(2, 2) = 0.5 * (1 - nubar);

//	D *= Ebar / (1. - nubar * nubar);
// }
//
// void RVELaw::GetLawFeatures(Features& rFeatures)
// {
// 	rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
// 	rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
// 	rFeatures.mOptions.Set(ISOTROPIC);
// 	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
// 	rFeatures.mStrainSize = GetStrainSize();
// 	rFeatures.mSpaceDimension = WorkingSpaceDimension();
// }
//
// int RVELaw::Check(
// 	const Properties& rMaterialProperties,
// 	const GeometryType& rElementGeometry,
// 	const ProcessInfo& rCurrentProcessInfo)
// {
// 		if(!rMaterialProperties.Has(YOUNG_MODULUS))
// 		    KRATOS_THROW_ERROR(std::invalid_argument, "RVELaw - missing
// YOUNG_MODULUS", "");
// 		if(!rMaterialProperties.Has(POISSON_RATIO))
// 		    KRATOS_THROW_ERROR(std::invalid_argument, "RVELaw - missing
// POISSON_RATIO", "");
// 		if(!rMaterialProperties.Has(ISOTROPIC_DAMAGE_MODULUS))
// 		KRATOS_THROW_ERROR(std::invalid_argument, "RVELaw - missing
// ISOTROPIC_DAMAGE_MODULUS", "");
// 		if(!rMaterialProperties.Has(INFINITY_YIELD_STRESS))
// 		KRATOS_THROW_ERROR(std::invalid_argument, "RVELaw - missing
// INFINITY_YIELD_STRESS", "");
// 		if(rMaterialProperties[INFINITY_YIELD_STRESS] < 0)
// 		    KRATOS_THROW_ERROR(std::invalid_argument, "RVELaw -
// INFINITY_YIELD_STRESS must be positive", "");
// 		return 0;
// }

} /* namespace Kratos.*/
