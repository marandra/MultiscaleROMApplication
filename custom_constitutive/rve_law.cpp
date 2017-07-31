#include "rve_law.h"
#include "custom_utilities/qr_utility.h"

namespace Kratos
{
RVELaw::RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param)
    : mpRVEModelPart(mpModelPart)
{
    auto w_list = param["w"];
    auto B_list = param["B"];
    auto prop_id_list = param["props_id"];
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
                BK(c, m) = B_list[i][c][m].GetDouble();
            }
        }
        Properties::Pointer prop = mpRVEModelPart->pGetProperties(prop_id_list[i].GetInt());
        ConstitutiveLaw::Pointer pcl = prop->GetValue(CONSTITUTIVE_LAW)->Clone();
        mB_vec.push_back(BK);
        mIW_vec.push_back(w_list[i].GetDouble());
        mCL_vec.push_back(pcl);
        mPropId_vec.push_back(prop_id_list[i].GetInt());
    }
}

// constructor used by Clone()
RVELaw::RVELaw(ModelPart::Pointer mpModelPart,
               std::vector<Matrix> B_list,
               std::vector<double> IW_list,
               std::vector<ConstitutiveLaw::Pointer> CL_list,
               std::vector<int> prop_id_list)
    : mpRVEModelPart(mpModelPart), mB_vec(B_list), mIW_vec(IW_list),
      mCL_vec(CL_list), mPropId_vec(prop_id_list)
{
    KRATOS_WATCH("Inside Clone Constructor");
    KRATOS_WATCH(mB_vec[0]);
    KRATOS_WATCH(mIW_vec[0]);
    KRATOS_WATCH(mCL_vec[0]);
    KRATOS_WATCH(mPropId_vec[0]);
}

// Destructor
RVELaw::~RVELaw() {}

ConstitutiveLaw::Pointer RVELaw::Clone() const
{
    RVELaw::Pointer pnewCL = boost::make_shared<RVELaw>(mpRVEModelPart,
                                                        mB_vec, mIW_vec,
                                                        mCL_vec, mPropId_vec);
    return pnewCL;
}

void RVELaw::InitializeMaterial(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const Vector& rShapeFunctionsValues)
{
    KRATOS_WATCH("inside initialize material")
    typedef std::vector<ConstitutiveLaw::Pointer>::size_type t_vsize;
    for (t_vsize i = 0; i < mCL_vec.size(); i++)
    {
        KRATOS_WATCH("initialized CL")
        mCL_vec[i]->InitializeMaterial(mPropId_vec[i], rElementGeometry,
                                       rShapeFunctionsValues);
    }
}

void RVELaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_WATCH("inside calculate material response")

    const auto nr_points = mB_vec.size();
    const auto nr_modes = mB_vec[0].size2();
    const auto nr_comps = GetStrainSize(); //TODO check == mB_vec[0].size1()
    KRATOS_WATCH(nr_points);
    KRATOS_WATCH(nr_modes);
    KRATOS_WATCH(nr_comps);

    const Properties& mat_props = rValues.GetMaterialProperties();
    Vector & strain = rValues.GetStrainVector();
    Vector & stress_bar = rValues.GetStressVector();
    Matrix & constit_matrix = rValues.GetConstitutiveMatrix();
    KRATOS_WATCH(strain);
    KRATOS_WATCH(stress_bar);
    KRATOS_WATCH(constit_matrix);
    Vector x = ZeroVector(nr_modes);
    Vector res = ZeroVector(nr_modes);
    Matrix A = ZeroMatrix(nr_modes, nr_modes);
    KRATOS_WATCH(res);
    KRATOS_WATCH(A);

    // row_major, col_mayor:order of the input matrix.
    // Should be col_major for the best performance.
    //enum storage_order {
    //    row_major,
    //    col_major
    //};
    QR<double, storage_order::row_major> QR_decomposition; // QR decomposition object
    KRATOS_WATCH(storage_order::row_major);
    KRATOS_WATCH(storage_order::col_major);
    CalculateResidual(A, res, x, rValues);
    // constit_matrix, props_list, model_part, geomParameters& rValues)
    // A, res, homog_stress = CalculateResidual(x, strain, stress_bar,
    // constit_matrix, props_list, model_part, geomParameters& rValues)
    auto tmp_nr_modes = 2;
    A(0,0) = 1.; A(0,1) = 2.;
    A(1,0) = 3.; A(1,1) = 4.;
    res(0) = 1.; res(1) = 2.;
    double aux_qr_A[tmp_nr_modes][tmp_nr_modes];
    double aux_qr_res[tmp_nr_modes];
    double aux_qr_Dx[tmp_nr_modes];
    for (auto ii = 0; ii < tmp_nr_modes; ii++){
        for (auto jj = 0; jj < tmp_nr_modes; jj++){
            aux_qr_A[ii][jj] = A(ii, jj);
        }
        aux_qr_res[ii] = res(ii);
    }
    KRATOS_WATCH(A);
    KRATOS_WATCH(aux_qr_A[0][0]);
    KRATOS_WATCH(aux_qr_A[0][1]);
    KRATOS_WATCH(aux_qr_A[1][0]);
    KRATOS_WATCH(aux_qr_A[1][1]);
    KRATOS_WATCH(aux_qr_res[0]);
    KRATOS_WATCH(aux_qr_res[1]);
    QR_decomposition.compute(tmp_nr_modes, tmp_nr_modes, &(aux_qr_A[0][0]));
    KRATOS_WATCH(aux_qr_A[0][0]);
    KRATOS_WATCH(aux_qr_A[0][1]);
    KRATOS_WATCH(aux_qr_A[1][0]);
    KRATOS_WATCH(aux_qr_A[1][1]);
    QR_decomposition.solve(&(aux_qr_res[0]), &(aux_qr_Dx[0]));
    KRATOS_WATCH(aux_qr_Dx[0]);
    KRATOS_WATCH(aux_qr_Dx[1]);

   /*
    int it = 1;
    double norm_res = 1.;
    while (norm_res > 1e-9 and it < 10)
    {
        Vector& Dx();
        QR_decomposition.compute(nr_modes, nr_modes, &(*A)(0, 0));
        QR_decomposition.solve(&(*res)(0), &(*Dx)(0));
    x -= Dx;
    A, res,
        homog_stress = calculate_residual(x, strain, iw_list, CL_list, B_list,
                                          props_list, model_part, geom);
    norm_res = np.linalg.norm(res, ord = 2)
                   print("RESIDUAL CRITERION :: norm res: {:.3e}".format(norm_res))
                       it += 1;
    }
*/
}

void RVELaw::CalculateResidual(Matrix &A, Vector &b, Vector &res, Parameters& rValues)
{
  const auto nr_points = mB_vec.size();
  const auto nr_modes = mB_vec[0].size2();
  const auto nr_comps = GetStrainSize();
  const auto dim = WorkingSpaceDimension();
/*
  Vector & strain_macro = rValues.GetStrainVector();
  Vector & stress = rValues.GetStressVector();
  Matrix & constit_matrix = rValues.GetConstitutiveMatrix();
  KRATOS_WATCH(strain_macro);
  const Properties& mat_props = rValues.GetMaterialProperties();

  for (auto i = 0; i < nr_points; i++)
  {
    Matrix B = mB_vec[i];
    Vector strain = strain_macro + prod(B, res);

    //TODO make this properly
    Vector N = ZeroVector(dim);
    Matrix F(dim,dim) = ZeroMatrix(dim, dim);
    F(0,0) = 1.0 + strain(0);
    F(0,1) = 0.5 * strain(3);
    F(0,2) = 0.5 * strain(5);
    F(1,0) = 0.5 * strain(3);
    F(1,1) = 1.0 + strain(1);
    F(1,2) = 0.5 * strain(4);
    F(2,0) = 0.5 * strain(5);
    F(2,1) = 0.5 * strain(4);
    F(2,2) = 1.0 + strain(2);
    //TODO compute det(F)
    double detF = 1.;

    Matrix DN_DX(3,2);

    mCL_vec[i].CalculateMaterialResponseCauchy(rValues);

    constit_matrix = rValues.GetConstitutiveMatrix();
    stress = rValues.GetStressVector();
  }
*/
}

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
int RVELaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo) {
     if (mB_vec[0].size1() != GetStrainSize())
         KRATOS_THROW_ERROR(std::invalid_argument, "Number of rows in modes matrix "
                 "rows differs from number of components of constitutive law", "");
     return 0;
    }

} /* namespace Kratos.*/
