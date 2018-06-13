#include "rve_law.h"
#include "custom_utilities/qr_utility.h"
#include <multiscale_rom_application_variables.h>
#include <utilities/math_utils.h>
namespace Kratos
{
// Constructor
RVELaw::RVELaw()
{
    //mpRVEModelPart = make_shared<ModelPart>("aaa");
    //mpRVEModelPart = Kernel::GetCurrentModel().CreateModelPart("aaa")
    //ReadMaterialsProcesn(param, mpRVEModelPart).Execute()
/*
    // TODO why can't these three be const?
    auto w_list = param["w"];
    auto B_list = param["B"];
    auto prop_id_list = param["props_id"];
    const auto nr_points = B_list.size();
    const auto nr_modes = B_list[0][0].size();
    const auto nr_comps = B_list[0].size();;
    // TODO(marcelo): Use GetStraubSize and Workingspacedimension properly
    //const auto nr_comps = GetStrainSize();
    //KRATOS_WATCH(nr_points);
    //KRATOS_WATCH(nr_modes);
    //KRATOS_WATCH(nr_comps);

    for (auto i = 0; i < nr_points; i++)
    {
        Matrix BK(nr_comps, nr_modes);
        for (auto c = 0; c < nr_comps; c++)
            for (auto m = 0; m < nr_modes; m++)
                BK(c, m) = B_list[i][c][m].GetDouble();
        Properties::Pointer prop =
            mpRVEModelPart->pGetProperties(prop_id_list[i].GetInt());
        ConstitutiveLaw::Pointer pcl = prop->GetValue(CONSTITUTIVE_LAW)->Clone();
        //KRATOS_WATCH(prop->GetValue(YOUNG_MODULUS));
        mB_vec.push_back(BK);
        mIW_vec.push_back(w_list[i].GetDouble());
        mCL_vec.push_back(pcl);
        mPropId_vec.push_back(prop_id_list[i].GetInt());
    }
    // TODO(marcelo): Discuss with Riccardo, following two != ZeroVector
    // mModesWeight.clear()
    // noalias(mModesWeights) = ZeroVector(nr_modes);
    mModesWeights = ZeroVector(nr_modes);
    */
}

// Default constructor
RVELaw::RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param)
    : mpRVEModelPart(mpModelPart)
{
    //mpRVEModelPart = make_shared<ModelPart>("aaa");
    //mpRVEModelPart = Kernel::GetCurrentModel().CreateModelPart("aaa")
    //ReadMaterialsProcesn(param, mpRVEModelPart).Execute()

    // TODO why can't these three be const?
    auto w_list = param["w"];
    auto B_list = param["B"];
    auto prop_id_list = param["props_id"];
    const auto nr_points = B_list.size();
    const auto nr_modes = B_list[0][0].size();
    const auto nr_comps = B_list[0].size();;
    // TODO(marcelo): Use GetStraubSize and Workingspacedimension properly
    //const auto nr_comps = GetStrainSize();
    //KRATOS_WATCH(nr_points);
    //KRATOS_WATCH(nr_modes);
    //KRATOS_WATCH(nr_comps);

    for (auto i = 0; i < nr_points; i++)
    {
        Matrix BK(nr_comps, nr_modes);
        for (auto c = 0; c < nr_comps; c++)
            for (auto m = 0; m < nr_modes; m++)
                BK(c, m) = B_list[i][c][m].GetDouble();
        Properties::Pointer prop =
            mpRVEModelPart->pGetProperties(prop_id_list[i].GetInt());
        ConstitutiveLaw::Pointer pcl = prop->GetValue(CONSTITUTIVE_LAW)->Clone();
        //KRATOS_WATCH(prop->GetValue(YOUNG_MODULUS));
        mB_vec.push_back(BK);
        mIW_vec.push_back(w_list[i].GetDouble());
        mCL_vec.push_back(pcl);
        mPropId_vec.push_back(prop_id_list[i].GetInt());
    }
    // TODO(marcelo): Discuss with Riccardo, following two != ZeroVector
    // mModesWeight.clear()
    // noalias(mModesWeights) = ZeroVector(nr_modes);
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
    const auto nr_modes = mB_vec[0].size2();
    mModesWeights = ZeroVector(nr_modes);
}

// Destructor
RVELaw::~RVELaw()
{
}

// Create
ConstitutiveLaw::Pointer RVELaw::Create(Kratos::Parameters Params) const
{
    //RVELaw::Pointer pnewCL = boost::make_shared<RVELaw>(
    //        mpRVEModelPart, mB_vec, mIW_vec, mCL_vec, mPropId_vec);
    RVELaw::Pointer p_clone(new RVELaw(mpRVEModelPart, mB_vec, mIW_vec, mCL_vec, mPropId_vec));
    return p_clone;
}

// Clone
ConstitutiveLaw::Pointer RVELaw::Clone() const
{
    //RVELaw::Pointer pnewCL = boost::make_shared<RVELaw>(
    //        mpRVEModelPart, mB_vec, mIW_vec, mCL_vec, mPropId_vec);
    RVELaw::Pointer p_clone(new RVELaw(mpRVEModelPart, mB_vec, mIW_vec, mCL_vec, mPropId_vec));
    return p_clone;
}

// Copy
RVELaw::RVELaw(const RVELaw& rOther) : ConstitutiveLaw(rOther)
{
}

void RVELaw::InitializeMaterial(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry,
                                const Vector& rShapeFunctionsValues)
{
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
    const auto nr_points = mB_vec.size();
    const auto nr_modes = mB_vec[0].size2();
    const auto nr_comps = mB_vec[0].size1();;
    // TODO(marcelo): Use GetStraubSize and Workingspacedimension properly
    //const auto nr_comps = GetStrainSize();
    //const Properties& mat_props = rValues.GetMaterialProperties();
    const Vector& strain_macro = rValues.GetStrainVector();

    Vector& homog_stress = rValues.GetStressVector();
    noalias(homog_stress) = ZeroVector(nr_comps);
    Matrix& homog_C = rValues.GetConstitutiveMatrix();
    noalias(homog_C) = ZeroMatrix(nr_comps, nr_comps);

    Matrix A(nr_modes, nr_modes);
    Vector res(nr_modes);
    Vector Dx(nr_modes);

    accumulate(A, res, strain_macro);
    double norm_res = 1.;
    int it = 1;
    while (norm_res > 1e-9 and it < 10)
    {
        solve(A, res, Dx);
        mModesWeights -= Dx;
        accumulate(A, res, strain_macro);
        norm_res = norm_2(res);
        KRATOS_WATCH(norm_res);
        it++;
    }

    // TODO(marcelo): Use flags for the computation of
    // Homogenize stress and constitutive tensor
    Matrix homog_C_taylor = ZeroMatrix(nr_comps, nr_comps);
    Matrix homog_C_fluct = ZeroMatrix(nr_comps, nr_comps);
    Matrix homog_Q = ZeroMatrix(nr_modes, nr_comps);
    Matrix homog_Op = ZeroMatrix(nr_modes, nr_comps);
    Matrix invA = ZeroMatrix(nr_modes, nr_modes);
    double vol_rve = 0.;
    double dummy_det;

    MathUtils<double>::InvertMatrix(A, invA, dummy_det);
    for (auto i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        // TODO(marcelo): strain argument should be const
        calculate_individual_material_response(stress, constit, strain, i);
        homog_stress += mIW_vec[i] * stress;
        homog_C_taylor += mIW_vec[i] * constit;
        homog_Q += mIW_vec[i] * prod(trans(mB_vec[i]), constit);
        vol_rve += mIW_vec[i];
    }
    homog_stress /= vol_rve;
    homog_Op = - prod(invA, homog_Q);
    Matrix homog_C_fluct_aux(nr_comps, nr_modes);
    for (auto i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        // TODO(marcelo): strain argument should be const
        calculate_individual_material_response(stress, constit, strain, i);
        homog_C_fluct_aux += mIW_vec[i] * prod(constit, mB_vec[i]);
    }
    noalias(homog_C_fluct) = prod(homog_C_fluct_aux, homog_Op);
    homog_C = homog_C_taylor + homog_C_fluct;
    homog_C /= vol_rve;
}

void RVELaw::solve(const Matrix& A, const Vector& res, Vector& Dx)
{
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

    //` Solve
    for (auto ii = 0; ii < nr_modes; ii++)
    {
        for (auto jj = 0; jj < nr_modes; jj++)
        {
            aux_qr_A[ii][jj] = A(ii, jj);
        }
        aux_qr_res[ii] = res(ii);
    }
    QR_decomposition.compute(nr_modes, nr_modes, &(aux_qr_A[0][0]));
    QR_decomposition.solve(&(aux_qr_res[0]), &(aux_qr_Dx[0]));

    // Update
    for (auto ii = 0; ii < nr_modes; ii++)
    {
        Dx[ii] = aux_qr_Dx[ii];
    }
}

void RVELaw::accumulate(Matrix& A, Vector& res, const Vector& strain_macro)
{
    const auto nr_points = mB_vec.size();
    const auto nr_modes = mB_vec[0].size2();
    const auto nr_comps = mB_vec[0].size1();;
    // TODO(marcelo): Use GetStraubSize and Workingspacedimension properly
    //const auto nr_comps = GetStrainSize();
    Matrix Aux1(nr_comps, nr_modes);

    noalias(A) = ZeroMatrix(nr_modes, nr_modes);
    noalias(res) = ZeroVector(nr_modes);
    for (auto i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);  // output
        Matrix constit(nr_comps, nr_comps);  // output
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        // TODO(marcelo): strain should be const
        calculate_individual_material_response(stress, constit, strain, i);
        // TODO(marcelo): explicitly write triple product for A
        // Dij = BTij Ckl Blj = for k for l for j for i
        noalias(Aux1) = prod(constit, mB_vec[i]);
        noalias(A) += mIW_vec[i] * prod(trans(mB_vec[i]), Aux1);
        noalias(res) += mIW_vec[i] * prod(trans(mB_vec[i]), stress);
    }
}

void RVELaw::calculate_individual_material_response(Vector& stress,
                                                    Matrix& constit,
                                                    Vector& strain,
                                                    std::size_t i)
{
    // create and pass individual parameters
    const auto dim = WorkingSpaceDimension();
    Flags cl_flags;
    cl_flags.Set(COMPUTE_STRESS, true);
    cl_flags.Set(COMPUTE_CONSTITUTIVE_TENSOR, true);

    Vector N(dim);
    Matrix DN_DX(dim, 2);
    Matrix F(dim, dim);
    F(0, 0) = 1.0 + strain(0);
    F(0, 1) = 0.5 * strain(3);
    F(0, 2) = 0.5 * strain(5);
    F(1, 0) = 0.5 * strain(3);
    F(1, 1) = 1.0 + strain(1);
    F(1, 2) = 0.5 * strain(4);
    F(2, 0) = 0.5 * strain(5);
    F(2, 1) = 0.5 * strain(4);
    F(2, 2) = 1.0 + strain(2);
    double detF = MathUtils<double>::Det(F);

    Parameters cl_params;
    cl_params.SetOptions(cl_flags);
    cl_params.SetDeformationGradientF(F);
    cl_params.SetDeterminantF(detF);
    cl_params.SetStrainVector(strain);
    cl_params.SetStressVector(stress);
    cl_params.SetConstitutiveMatrix(constit);
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
    if (rThisVariable == REDUCED_MODES_WEIGHTS)
        return true;
    return false;
}

Vector& RVELaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == REDUCED_MODES_WEIGHTS)
        rValue = mModesWeights;
    return rValue;
}

int RVELaw::Check(const Properties& rMaterialProperties,
                  const GeometryType& rElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo)
{
    // Self check
    // TODO(marcelo): check disable until proper handlo od 2D-3D
    //if (mB_vec[0].size1() != GetStrainSize())
    //    KRATOS_THROW_ERROR(
    //        std::invalid_argument,
    //        "Number of rows in modes matrix "
    //        "rows differs from number of components of constitutive law",
    //        "");

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

//************************************************************************************
//************************************************************************************

void RVELaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mpRVEModelPart", mpRVEModelPart);
    rSerializer.save("mB_vec", mB_vec);
    rSerializer.save("mIW_vec", mIW_vec);
    rSerializer.save("mCL_vec", mCL_vec);
    rSerializer.save("mPropId_vec", mPropId_vec);
    rSerializer.save("mModesWeights", mModesWeights);
}

//************************************************************************************
//************************************************************************************

void RVELaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mpRVEModelPart", mpRVEModelPart);
    rSerializer.load("mB_vec", mB_vec);
    rSerializer.load("mIW_vec", mIW_vec);
    rSerializer.load("mCL_vec", mCL_vec);
    rSerializer.load("mPropId_vec", mPropId_vec);
    rSerializer.load("mModesWeights", mModesWeights);
}

} /* namespace Kratos.*/
