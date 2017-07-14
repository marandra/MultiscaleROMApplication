#include "rve_law.h"

namespace Kratos
{
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

        // NOTE: here the rMaterialProperties come from the MACROSCALE.
        // We are assuming, however, that it also contains the materials
        // to be used in the microscale
        auto prop = mpRVEModelPart->pGetProperties(prop_id[i].GetInt());
        mprop_list.push_back(prop);
        ConstitutiveLaw::Pointer pcl = prop->GetValue(CONSTITUTIVE_LAW)->Clone();
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


void RVELaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
  const Properties& matprops = rValues.GetMaterialProperties();
  Vector& epsilon_h = rValues.GetStrainVector();
  Vector& sigma_bar = rValues.GetStressVector();
  Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

  KRATOS_WATCH("inside calculate material response")

    unsigned int nr_points = mB_list.size();
    unsigned int nr_comps = mB_list[0].size();
    unsigned int nr_modes = mB_list[0][0].size();

    A, res, homog_stress = CalculateResidual(x, epsilon_h, sigma_bar, constitutive_matrix, props_list, model_part, geom)
    it = 1
    norm_res = 1
    while(norm_res > 1e-9 and it < 10):
        Dx = -np.linalg.solve(A, res)
        x += Dx
        A, res, homog_stress = calculate_residual(x, epsilon_h, iw_list,
                CL_list, B_list, props_list, model_part, geom)
        norm_res = np.linalg.norm(res, ord=2)
        print("RESIDUAL CRITERION :: norm res: {:.3e}".format(norm_res))
        it += 1
    print("Convergence is achieved (or not)")

}

void CalculateResidual(x, epsilon_h, sigma_bar, constitutive_matrix, props_list, model_part, geom)
{
  //mIW_list[i];
  //mCL_list[i];
  unsigned int nr_points = mB_list.size();
  unsigned int nr_comps = mB_list[0].size();
  unsigned int nr_modes = mB_list[0][0].size();

  Matrix A = ZeroMatrix(nr_modes, nr_comps);
  Vector b = ZeroVector(nr_modes);

  for (unsigned int i = 0; i < nr_points; i++)
  {
    Matrix B = mB_list[i];
    Vector epsilon = epsilon_h + prod(B, x);

    Vector N = ZeroVector(3);
    Matrix F(3,3) = ZeroMatrix(3, 3);
    F(0,0) = 1.0 + epsilon(0);   F(0,1) = 0.5 * epsilon(3); F(0,2) = 0.5 * epsilon(5w);
    F(1,0) = 0.5 * epsilon(3);   F(1,1) = 1.0 + epsilon(1); F(1,2) = 0.5 * epsilon(4w);
    F(2,0) = 0.5 * epsilon(5);   F(2,1) = 0.5 * epsilon(4); F(2,2) = 1.0 + epsilon(2w);
    //TODO compute det(F)
    detF = 1.;

    Matrix DN_DX(3,2);
    constitutive_matrix = km.Matrix(cl.GetStrainSize(), cl.GetStrainSize())
    stress_vector = km.Vector(cl.GetStrainSize())
    strain_vector = km.Vector(cl.GetStrainSize())
    for i in range(cl.GetStrainSize()):
    stress_vector[i] = 0.
    strain_vector[i] = epsilon[i]

#setting the parameters - note that a constitutive law may not need them all!
    cl_params = km.ConstitutiveLawParameters()
    cl_params.SetOptions(cl_options)
    cl_params.SetDeformationGradientF(F)
    cl_params.SetDeterminantF(detF)
    cl_params.SetStrainVector(strain_vector)
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetShapeFunctionsValues(N)
    cl_params.SetShapeFunctionsDerivatives(DN_DX)
    cl_params.SetProcessInfo(process_info)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    cl.CalculateMaterialResponseCauchy(cl_params)

    CM = cl_params.GetConstitutiveMatrix()
    stress = cl_params.GetStressVector()
    size = cl.GetStrainSize()
    CM_np = np.empty((size, size))
    stress_np = np.empty(size)
    for i in range(size):
    stress_np[i] = stress[i]
    for j in range(size):
    CM_np[i, j] = CM[i, j]

    return stress_np, CM_np


  }
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
