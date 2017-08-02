#include <geometries/triangle_2d_3.h>
#include "rve_law.h"
#include "custom_utilities/qr_utility.h"

namespace Kratos
{
//Default constructor
RVELaw::RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param)
    : mpRVEModelPart(mpModelPart)
{
    // TODO why can't these three be const?
    auto w_list = param["w"];
    auto B_list = param["B"];
    auto prop_id_list = param["props_id"];
    const auto nr_points = B_list.size();
    const auto nr_modes = B_list[0][0].size();
    const auto nr_comps = GetStrainSize();

    for (auto i = 0; i < nr_points; i++)
    {
        Matrix BK(nr_comps, nr_modes);
        for (auto c = 0; c < nr_comps; c++)
            for (auto m = 0; m < nr_modes; m++)
                BK(c, m) = B_list[i][c][m].GetDouble();
        Properties::Pointer prop = mpRVEModelPart->pGetProperties(prop_id_list[i].GetInt());
        ConstitutiveLaw::Pointer pcl = prop->GetValue(CONSTITUTIVE_LAW)->Clone();
        mB_vec.push_back(BK);
        mIW_vec.push_back(w_list[i].GetDouble());
        mCL_vec.push_back(pcl);
        mPropId_vec.push_back(prop_id_list[i].GetInt());
    }
    // TODO why is ZeroVector 'const double' type, instead of 'Vector'?
    mModesWeights = ZeroVector(nr_modes);
}

// Constructor used by Clone()
RVELaw::RVELaw(ModelPart::Pointer mpModelPart, std::vector<Matrix> B_list,
               std::vector<double> IW_list, std::vector<ConstitutiveLaw::Pointer> CL_list,
               std::vector<int> prop_id_list)
    : mpRVEModelPart(mpModelPart), mB_vec(B_list), mIW_vec(IW_list),
      mCL_vec(CL_list), mPropId_vec(prop_id_list)
{
    KRATOS_WATCH("Inside Clone Constructor");
    KRATOS_WATCH(mB_vec[0]);
    KRATOS_WATCH(mIW_vec[0]);
    KRATOS_WATCH(mCL_vec[0]);
    KRATOS_WATCH(mPropId_vec[0]);
    const auto nr_modes = mB_vec[0].size2();
    mModesWeights = ZeroVector(nr_modes);
}

// Destructor
RVELaw::~RVELaw() {}

// Clone
ConstitutiveLaw::Pointer RVELaw::Clone() const
{
    RVELaw::Pointer pnewCL = boost::make_shared<RVELaw>(mpRVEModelPart,
                                                        mB_vec, mIW_vec,
                                                        mCL_vec, mPropId_vec);
    return pnewCL;
}

// Copy
RVELaw::RVELaw(const RVELaw& rOther)
        : ConstitutiveLaw(rOther)
{
}


void RVELaw::InitializeMaterial(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const Vector& rShapeFunctionsValues)
{
    for (auto i = 0; i < mCL_vec.size(); i++)
    {
        KRATOS_WATCH("initialized CL - fix argument. see TODO")
        // TODO: pass Property, not Id
        mCL_vec[i]->InitializeMaterial(mPropId_vec[i], rElementGeometry,
                                       rShapeFunctionsValues);
    }
}

void RVELaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_WATCH("inside calculate material response")

    const auto nr_points = mB_vec.size();
    const auto nr_modes = mB_vec[0].size2();
    const auto nr_comps = GetStrainSize();
    const auto dim = WorkingSpaceDimension();

    const Properties& mat_props = rValues.GetMaterialProperties();
    Vector & strain_macro = rValues.GetStrainVector();
    Vector & stress_rve = rValues.GetStressVector();
    Matrix & c_matrix_rve = rValues.GetConstitutiveMatrix();

    Vector res = ZeroVector(nr_modes);
    Matrix A = ZeroMatrix(nr_modes, nr_modes);
    double aux_qr_A[nr_modes][nr_modes];
    double aux_qr_res[nr_modes];
    double aux_qr_Dx[nr_modes];

    // row_major, col_mayor:order of the input matrix.
    // Should be col_major for the best performance.
    //enum storage_order {
    //    row_major,
    //    col_major
    //};
    //KRATOS_WATCH(storage_order::row_major);
    //KRATOS_WATCH(storage_order::col_major);
    QR<double, storage_order::row_major> QR_decomposition;
    /*
    // temporary validation od QR decomposition - solve
    auto tmp_nr_modes = 3;
    A.resize(3,3);
    A(0,0) = 1.; A(0,1) = 2.; A(0,2) = 3.;
    A(1,0) = 3.; A(1,1) = 2.; A(1,2) = 6.;
    A(2,0) = 2.; A(2,1) = 3.; A(2,2) = 2.;
    res(0) = 2.; res(1) = 3.; res(2) = 1.;
    for (auto ii = 0; ii < tmp_nr_modes; ii++){
        for (auto jj = 0; jj < tmp_nr_modes; jj++){
            aux_qr_A[ii][jj] = A(ii, jj);
        }
        aux_qr_res[ii] = res(ii);
    }
    KRATOS_WATCH(aux_qr_A[0][0]);
    KRATOS_WATCH(aux_qr_A[0][1]);
    KRATOS_WATCH(aux_qr_A[0][2]);
    KRATOS_WATCH(aux_qr_A[1][0]);
    KRATOS_WATCH(aux_qr_A[1][1]);
    KRATOS_WATCH(aux_qr_A[1][2]);
    KRATOS_WATCH(aux_qr_A[2][0]);
    KRATOS_WATCH(aux_qr_A[2][1]);
    KRATOS_WATCH(aux_qr_A[2][2]);
    QR_decomposition.compute(tmp_nr_modes, tmp_nr_modes, &(aux_qr_A[0][0]));
    KRATOS_WATCH(aux_qr_A[0][0]);
    KRATOS_WATCH(aux_qr_A[0][1]);
    KRATOS_WATCH(aux_qr_A[0][2]);
    KRATOS_WATCH(aux_qr_A[1][0]);
    KRATOS_WATCH(aux_qr_A[1][1]);
    KRATOS_WATCH(aux_qr_A[1][2]);
    KRATOS_WATCH(aux_qr_A[2][0]);
    KRATOS_WATCH(aux_qr_A[2][1]);
    KRATOS_WATCH(aux_qr_A[2][2]);
    QR_decomposition.solve(&(aux_qr_res[0]), &(aux_qr_Dx[0]));
    KRATOS_WATCH(aux_qr_Dx[0]);
    KRATOS_WATCH(aux_qr_Dx[1]);
    KRATOS_WATCH(aux_qr_Dx[2]);
    */

    // Caclulate material response of every CL
    for (auto i=0; i < nr_points; i++){
        // create and pass individual parameters
        Flags cl_flags;
        cl_flags.Set(COMPUTE_STRESS, true);
        cl_flags.Set(COMPUTE_CONSTITUTIVE_TENSOR, true);
        Vector stress;
        Matrix c_matrix;
        Vector strain;
        noalias(strain) = strain_macro + prod(mB_vec[i], res);
        Vector N = ZeroVector(dim);
        Matrix DN_DX(3,2);
        Matrix F(dim, dim);
        F(0,0) = 1.0 + strain(0);
        F(0,1) = 0.5 * strain(3);
        F(0,2) = 0.5 * strain(5);
        F(1,0) = 0.5 * strain(3);
        F(1,1) = 1.0 + strain(1);
        F(1,2) = 0.5 * strain(4);
        F(2,0) = 0.5 * strain(5);
        F(2,1) = 0.5 * strain(4);
        F(2,2) = 1.0 + strain(2);
        double detF = determinant(F);
        Parameters cl_params;
        cl_params.SetOptions(cl_flags);
        cl_params.SetDeformationGradientF(F);
        cl_params.SetDeterminantF(detF);
        cl_params.SetStrainVector(strain);
        cl_params.SetStressVector(stress);
        cl_params.SetConstitutiveMatrix(c_matrix);
        cl_params.SetShapeFunctionsValues(N);
        cl_params.SetShapeFunctionsDerivatives(DN_DX);
        cl_params.SetProcessInfo(mpRVEModelPart->GetProcessInfo());
        cl_params.SetMaterialProperties(mpRVEModelPart->pGetProperties(mPropId_vec[i]);
        // TODO: see how to pass proper geom and how it is used
        cl_params.SetElementGeometry(Triangle2D3);


        // get response
    }

   //const auto nr_modes = mB_vec[0].size2();


    /*
    int it = 1;
    double norm_res = 1.;
    while (norm_res > 1e-9 and it < 10)
    {
        // Solve
        for (auto ii = 0; ii < nr_modes; ii++){
            for (auto jj = 0; jj < nr_modes; jj++){
                aux_qr_A[ii][jj] = A(ii, jj);
            }
            aux_qr_res[ii] = res(ii);
        }
        QR_decomposition.compute(nr_modes, nr_modes, &(aux_qr_A[0][0]));
        QR_decomposition.solve(&(aux_qr_res[0]), &(aux_qr_Dx[0]));

        // Update
        for (auto ii = 0; ii < nr_modes; ii++) {
            mModesWeights[ii] -= aux_qr_Dx[ii];
        }

        //CalculateResidual(A, res, rValues);
        norm_res = norm_2(res);
        KRATOS_WATCH(norm_res);
        it++;
        break;
    }
    */
}

void RVELaw::CalculateResidual(Matrix &A, Vector &res, Parameters& rValues)
{
  const auto nr_points = mB_vec.size();
  const Vector & strain_macro = rValues.GetStrainVector();
  KRATOS_WATCH(strain_macro);
  const auto dim = WorkingSpaceDimension();
  const auto nr_comps = GetStrainSize();
    //const auto nr_modes = mB_vec[0].size2();
  //const Properties& mat_props = rValues.GetMaterialProperties();

  for (auto i = 0; i < nr_points; i++)
  {
    const Matrix B = mB_vec[i];
    const Vector strain = strain_macro + prod(B, res);

    //TODO make this properly
    Vector N = ZeroVector(dim);
    Matrix DN_DX(3,2);
    Matrix F(dim, dim);
    F(0,0) = 1.0 + strain(0);
    F(0,1) = 0.5 * strain(3);
    F(0,2) = 0.5 * strain(5);
    F(1,0) = 0.5 * strain(3);
    F(1,1) = 1.0 + strain(1);
    F(1,2) = 0.5 * strain(4);
    F(2,0) = 0.5 * strain(5);
    F(2,1) = 0.5 * strain(4);
    F(2,2) = 1.0 + strain(2);
    double detF = determinant(F);
    KRATOS_WATCH(F);
    KRATOS_WATCH(detF);

      /*
    Vector stress(nr_comps);
    FlagType cl_flags = Flags();
    //cl_flags.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    //cl_flags.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
    // Prepare parameters for CL
    Parameters cl_params;
    cl_params.SetOptions(cl_flags);
    cl_params.SetDeformationGradientF(F);
    cl_params.SetDeterminantF(detF);
    cl_params.SetStrainVector(strain);
    cl_params.SetStressVector(stress_vector)
    cl_params.SetConstitutiveMatrix(constitutive_matrix)
    cl_params.SetShapeFunctionsValues(N)
    cl_params.SetShapeFunctionsDerivatives(DN_DX)
    cl_params.SetProcessInfo(process_info)
    cl_params.SetMaterialProperties(properties)
    cl_params.SetElementGeometry(geom)

    // Compute RVE's point's response
    mCL_vec[i]->CalculateMaterialResponseCauchy(rValues);
    // Get RVE's point's response
    Vector & stress = rValues.GetStressVector();
      stress = rValues.GetStressVector();
    Matrix & constit_matrix = rValues.GetConstitutiveMatrix();
       */
  }
}


int RVELaw::determinant_sign(const permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}


double RVELaw::determinant(Matrix m) {
    permutation_matrix<std::size_t> pm(m.size1());
    double det = 1.0;
    if (lu_factorize(m, pm)) {
        det = 0.0;
    } else {
        for (int i = 0; i < m.size1(); i++)
            det *= m(i, i); // multiply by elements on diagonal
        det = det * determinant_sign(pm);
    }
    return det;
}


int RVELaw::Check(const Properties& rMaterialProperties,
                  const GeometryType& rElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo)
{
    if (mB_vec[0].size1() != GetStrainSize())
         KRATOS_THROW_ERROR(std::invalid_argument, "Number of rows in modes matrix "
                 "rows differs from number of components of constitutive law", "");

    //TODO: Implement call ->Check of every CL

    return 0;
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

} /* namespace Kratos.*/
