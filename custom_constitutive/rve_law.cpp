#include "rve_law.h"
#include "custom_utilities/qr_utility.h"
//#include <geometries/triangle_2d_3.h>
namespace Kratos
{



// Default constructor
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
        Properties::Pointer prop =
            mpRVEModelPart->pGetProperties(prop_id_list[i].GetInt());
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
RVELaw::RVELaw(ModelPart::Pointer mpModelPart,
               std::vector<Matrix> B_list,
               std::vector<double> IW_list,
               std::vector<ConstitutiveLaw::Pointer> CL_list,
               std::vector<int> prop_id_list)
    : mpRVEModelPart(mpModelPart),
      mB_vec(B_list),
      mIW_vec(IW_list),
      mCL_vec(CL_list),
      mPropId_vec(prop_id_list)
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
RVELaw::~RVELaw()
{
}



// Clone
ConstitutiveLaw::Pointer RVELaw::Clone() const
{
    RVELaw::Pointer pnewCL = boost::make_shared<RVELaw>(
        mpRVEModelPart, mB_vec, mIW_vec, mCL_vec, mPropId_vec);
    return pnewCL;
}



// Copy
RVELaw::RVELaw(const RVELaw& rOther) : ConstitutiveLaw(rOther)
{
}



void RVELaw::InitializeMaterial(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const Vector& rShapeFunctionsValues)
{
    KRATOS_WATCH("Initialize individual CLs")
    for (auto i = 0; i < mCL_vec.size(); i++)
    {
        const Properties& material_props =
            mpRVEModelPart->GetProperties(mPropId_vec[i]);

        // TODO We need geometry of the HF element, but we don have it.
        // Passing empty geometry, as is individual CL is not using it.
        const GeometryType dummy_element_geometry;
        mCL_vec[i]->InitializeMaterial(material_props, dummy_element_geometry,
                                       rShapeFunctionsValues);
    }
}



void RVELaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                  const GeometryType& rElementGeometry,
                                  const Vector& rShapeFunctionsValues,
                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_WATCH("Finalize solution step individual CL")
    for (auto i = 0; i < mCL_vec.size(); i++)
    {
        const Properties& material_props =
            mpRVEModelPart->GetProperties(mPropId_vec[i]);

        // TODO We need geometry of the HF element, but we don have it.
        // Passing empty geometry, as is individual CL is not using it.
        const GeometryType dummy_element_geometry;
        mCL_vec[i]->FinalizeSolutionStep(material_props, dummy_element_geometry,
                                       rShapeFunctionsValues, rCurrentProcessInfo);
    }
}



void RVELaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const auto nr_modes = mB_vec[0].size2();
    const Properties& mat_props = rValues.GetMaterialProperties();
    const Vector& strain_macro = rValues.GetStrainVector();

    Vector& stress_homog = rValues.GetStressVector();
    Matrix& c_matrix_homog = rValues.GetConstitutiveMatrix();

    Matrix A(nr_modes, nr_modes);
    Vector res(nr_modes);
    Vector Dx(nr_modes);

    accumulate(A, res, c_matrix_homog, stress_homog, strain_macro);
    double norm_res = 1.;
    int it = 1;
    while (norm_res > 1e-9 and it < 10)
    {
        solve(A, res, Dx);
        mModesWeights -= Dx;
        accumulate(A, res, c_matrix_homog, stress_homog, strain_macro);
        norm_res = norm_2(res);
        KRATOS_WATCH(norm_res);
        it++;
    }
}



void RVELaw::solve(const Matrix &A, const Vector &res, Vector &Dx){
    const auto nr_modes = mB_vec[0].size2();
    double aux_qr_A[nr_modes][nr_modes];
    double aux_qr_res[nr_modes];
    double aux_qr_Dx[nr_modes];

    // row_major, col_mayor:order of the input matrix.
    // Should be col_major for the best performance.
    // enum storage_order {
    //    row_major,
    //    col_major
    //};
    // KRATOS_WATCH(storage_order::row_major);
    // KRATOS_WATCH(storage_order::col_major);
    QR<double, storage_order::row_major> QR_decomposition;

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
            Dx[ii] = aux_qr_Dx[ii];
        }
    }



void
RVELaw::accumulate(Matrix &A, Vector &res, Matrix &c_matrix_homog, Vector &stress_homog, const Vector &strain_macro) {
    const auto nr_points = mB_vec.size();
    const auto nr_modes = mB_vec[0].size2();
    const auto nr_comps = GetStrainSize();
    res = ZeroVector(nr_modes);
    A = ZeroMatrix(nr_modes, nr_modes);
    for (auto i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);             // output
        Matrix c_matrix(nr_comps, nr_comps); // output
        Vector strain = strain_macro + prod(mB_vec[i], res);
        // TODO(marcelo): strain should be const
        calculate_individual_material_response(stress, c_matrix, strain, i);
        //TODO(marcelo): must use prod<temp_type>(...)
        // noalias(A) += mIW_vec[i] * prod(trans(mB_vec[i]), prod<temp_type>(c_matrix, mB_vec[i]));
        Matrix Aux1(nr_comps, nr_modes);
        A += mIW_vec[i] * Matrix(prod(trans(mB_vec[i]), prod(c_matrix, mB_vec[i], Aux1)));
        res += mIW_vec[i] * Vector(prod(trans(mB_vec[i]), stress));
        c_matrix_homog += mIW_vec[i] * c_matrix;
        stress_homog += mIW_vec[i] * stress;
    }
}



void RVELaw::calculate_individual_material_response(Vector &stress, Matrix &c_matrix, Vector &strain, std::size_t i)
{
    // create and pass individual parameters
    const auto dim = WorkingSpaceDimension();
    Flags cl_flags;
    cl_flags.Set(COMPUTE_STRESS, true);
    cl_flags.Set(COMPUTE_CONSTITUTIVE_TENSOR, true);

    Vector N = ZeroVector(dim);
    Matrix DN_DX(3, 2);
    Matrix F(dim, dim);
    F(0, 0) = 1.0 + strain(0); F(0, 1) = 0.5 * strain(3); F(0, 2) = 0.5 * strain(5);
    F(1, 0) = 0.5 * strain(3); F(1, 1) = 1.0 + strain(1); F(1, 2) = 0.5 * strain(4);
    F(2, 0) = 0.5 * strain(5); F(2, 1) = 0.5 * strain(4); F(2, 2) = 1.0 + strain(2);
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
    const Properties& material_props = mpRVEModelPart->GetProperties(mPropId_vec[i]);
    cl_params.SetMaterialProperties(material_props);
    // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
    // cl_params.SetElementGeometry();

    mCL_vec[i]->CalculateMaterialResponseCauchy(cl_params);
}



bool RVELaw::Has(const Variable<Vector>& rThisVariable)
{
    //if (rThisVariable == MODE_WEIGHT)
    //    return true;
    //else
        return false;
}



Vector& RVELaw::GetValue(const Variable<Vector>& rThisVariable, Vector&
rValue)
{
    //if (rThisVariable == MODE_WEIGHT)
    //    rValue = mModesWeights;
	return rValue;
}


int RVELaw::Check(const Properties& rMaterialProperties,
                  const GeometryType& rElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo)
{
    // Self check
    if (mB_vec[0].size1() != GetStrainSize())
        KRATOS_THROW_ERROR(
            std::invalid_argument,
            "Number of rows in modes matrix "
            "rows differs from number of components of constitutive law",
            "");

    // Individual CLs check
    const auto nr_points = mB_vec.size();
    for (auto i = 0; i < nr_points; i++)
    {
        const Properties& material_props =
            mpRVEModelPart->GetProperties(mPropId_vec[i]);
        const GeometryType dummy_element_geometry;
        mCL_vec[i]->Check(material_props, dummy_element_geometry, rCurrentProcessInfo);
    }

    return 0;
}

// bool RVELaw::Has(const Variable<double>& rThisVariable)
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
