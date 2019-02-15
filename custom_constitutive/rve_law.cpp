#include "rve_law.h"
#include "custom_utilities/qr_utility.h"
#include "multiscale_rom_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
/***********************************************************************************/
// Default constructor
/***********************************************************************************/
RVELaw::RVELaw()
{
}

/***********************************************************************************/
// Main constructor, used by Create
/***********************************************************************************/
RVELaw::RVELaw(Kratos::Parameters Params)
{
    // Parse RVE materials filename from Parameters
    Kratos::Parameters default_parameters(R"(
    {
        "name": "constitutive law name",
        "Parameters" : {
            "rve_materials_filename": "please specify the file to be opened",
            "rve_data_filename": "please specify the file to be opened",
            "convergence_criterion": "residual_criterion",
            "residual_relative_tolerance": 1e-4,
            "residual_absolute_tolerance": 1e-9,
            "max_iteration": 10,
            "verbose": 0
        }
    }  )"
    );
    Params.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Read json string in file, create parameters
    Kratos::Parameters materials_params(
            ReadFile(Params["Parameters"]["rve_materials_filename"].GetString()));
    Kratos::Parameters data_params(
            ReadFile(Params["Parameters"]["rve_data_filename"].GetString()));
    mRelativeTolerance = Params["Parameters"]["residual_relative_tolerance"].GetDouble();
    mAbsoluteTolerance = Params["Parameters"]["residual_absolute_tolerance"].GetDouble();
    mMaxIteration = Params["Parameters"]["max_iteration"].GetInt();
    mVerbose = Params["Parameters"]["verbose"].GetInt();

    // Parse material parameters and populate mpProperties
    GetPropertyBlock(materials_params);

    // Parse data parameters
    Kratos::Parameters B_list = data_params["B"];
    Kratos::Parameters w_list = data_params["w"];
    Kratos::Parameters prop_id_list = data_params["props_id"];
    const std::size_t nr_points = B_list.size();
    const std::size_t nr_modes = B_list[0][0].size();
    const std::size_t nr_comps = GetStrainSize();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Matrix BK(nr_comps, nr_modes);
        for (std::size_t c = 0; c < nr_comps; c++)
            for (std::size_t m = 0; m < nr_modes; m++)
                BK(c, m) = B_list[i][c][m].GetDouble();

        // Populate members
        mB_vec.push_back(BK);
        mIW_vec.push_back(w_list[i].GetDouble());
        mPropId_vec.push_back(prop_id_list[i].GetInt());
        Properties prop = mProperties_map[prop_id_list[i].GetInt()];
        ConstitutiveLaw::Pointer pcl = prop.GetValue(CONSTITUTIVE_LAW)->Clone();
        mCL_vec.push_back(pcl);

        //DEBUG
        //KRATOS_WATCH(prop.GetValue(YOUNG_MODULUS));
        //KRATOS_WATCH(prop.GetId());
        //KRATOS_WATCH(pcl->Info());
    }
    // preserve = false -> new elements (all of them in this case) not initialized
    mModesWeights.resize(nr_modes, false);
    mModesWeights.clear(); // values initialized to zero
}

/***********************************************************************************/
// Constructor used by Clone()
/***********************************************************************************/
RVELaw::RVELaw(PropertiesMap pProperties,
               std::vector<Matrix> B_list,
               std::vector<double> IW_list,
               std::vector<ConstitutiveLaw::Pointer> CL_list,
               std::vector<int> prop_id_list,
               double abs_tol, double rel_tol, int max_iter, int verbose)
    : mProperties_map(pProperties),
      mB_vec(B_list),
      mIW_vec(IW_list),
      mCL_vec(CL_list),
      mPropId_vec(prop_id_list),
      mAbsoluteTolerance(abs_tol),
      mRelativeTolerance(rel_tol),
      mMaxIteration(max_iter),
      mVerbose(verbose)
{
    const std::size_t nr_modes = mB_vec[0].size2();
    // preserve = false -> new elements (all of them in this case) not initialized
    mModesWeights.resize(nr_modes, false);
    mModesWeights.clear(); // values initialized to zero
}

/***********************************************************************************/
// Destructor
/***********************************************************************************/
RVELaw::~RVELaw()
{
}

/***********************************************************************************/
// Create
/***********************************************************************************/
ConstitutiveLaw::Pointer RVELaw::Create(Kratos::Parameters Params) const
{
    return Kratos::make_shared<RVELaw>(Params);
}

/***********************************************************************************/
// Clone
/***********************************************************************************/
ConstitutiveLaw::Pointer RVELaw::Clone() const
{
    //RVELaw::Pointer p_clone(new RVELaw(mProperties_map, mB_vec, mIW_vec, mCL_vec, mPropId_vec));
    RVELaw::Pointer p_clone(new RVELaw(mProperties_map, mB_vec, mIW_vec, mCL_vec, mPropId_vec,
            mAbsoluteTolerance, mRelativeTolerance, mMaxIteration, mVerbose));
    return p_clone;
}

/***********************************************************************************/
// Copy
/***********************************************************************************/
RVELaw::RVELaw(const RVELaw& rOther) : ConstitutiveLaw(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/
std::string RVELaw::ReadFile(const std::string &filename) const
{
    std::ifstream infile(filename);
    KRATOS_ERROR_IF_NOT(infile.good()) << "File " << filename << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::GetPropertyBlock(Kratos::Parameters Materials)
{
    for (auto i = 0; i < Materials["properties"].size(); ++i) {
        Kratos::Parameters material = Materials["properties"].GetArrayItem(i);
        AssignPropertyBlock(material);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void RVELaw::TrimComponentName(std::string& rLine){
    std::stringstream ss(rLine);
    std::size_t counter = 0;
    while (std::getline(ss, rLine, '.')){counter++;}
    if (counter > 1)
        KRATOS_WARNING("RVE Law") << "Ignoring module information for component " << rLine << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::AssignPropertyBlock(Kratos::Parameters Data)
{
    const std::size_t property_id = Data["properties_id"].GetInt();
    Properties property(property_id);

    //Set the CONSTITUTIVE_LAW for the current p_properties.
    if (Data["Material"].Has("constitutive_law")) {
        std::string constitutive_law_name = Data["Material"]["constitutive_law"]["name"].GetString();
        TrimComponentName(constitutive_law_name);
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(constitutive_law_name).Clone();
        property.SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
    } else {
        KRATOS_WARNING("RVE Law") << "No constitutive law defined for material ID: " << property_id << std::endl;
    }

    // Add / override the values of material parameters in the p_properties
    Kratos::Parameters variables = Data["Material"]["Variables"];
    for(auto iter = variables.begin(); iter != variables.end(); iter++) {
        const Kratos::Parameters value = variables.GetValue(iter.name());

        std::string variable_name = iter.name();
        TrimComponentName(variable_name);

        // TODO: Reuse this block from read_material_utility
        // We don't just copy the values, we do some transformation depending of the destination variable
        if(KratosComponents<Variable<double> >::Has(variable_name)) {
            const Variable<double>& variable = KratosComponents<Variable<double>>().Get(variable_name);
            if (value.IsDouble()) {
                property.SetValue(variable, value.GetDouble());
            } else if (value.IsInt()) {
                property.SetValue(variable, static_cast<double>(value.GetInt()));
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
            const Variable<bool>& variable = KratosComponents<Variable<bool>>().Get(variable_name);
            if (value.IsBool()) {
                property.SetValue(variable, value.GetBool());
            } else if (value.IsInt()) {
                property.SetValue(variable, static_cast<bool>(value.GetInt()));
            } else if (value.IsDouble()) {
                property.SetValue(variable, static_cast<bool>(value.GetDouble()));
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
            const Variable<int>& variable = KratosComponents<Variable<int>>().Get(variable_name);
            if (value.IsInt()) {
                property.SetValue(variable, value.GetInt());
            } else if (value.IsDouble()) {
                property.SetValue(variable, static_cast<int>(value.GetDouble()));
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
//        } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
//            const Variable<array_1d<double, 3>>& variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(variable_name);
//            if (value.IsVector()) {
//                array_1d<double, 3> temp(3, 0.0);
//                const Vector& value_variable = value.GetVector();
//                const std::size_t iter_number = (3 < value_variable.size()) ? 3 : value_variable.size();
//                for (std::size_t index = 0; index < iter_number; index++)
//                    temp[index] = value_variable[index];
//                property.SetValue(variable, temp);
//            } else {
//                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
//            }
//        } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
//            const Variable<array_1d<double, 6>>& variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(variable_name);
//            if (value.IsVector()) {
//                array_1d<double, 6> temp(6, 0.0);
//                const Vector& value_variable = value.GetVector();
//                const std::size_t iter_number = (6 < value_variable.size()) ? 6 : value_variable.size();
//                for (std::size_t index = 0; index < iter_number; index++)
//                    temp[index] = value_variable[index];
//                property.SetValue(variable, temp);
//            } else {
//                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
//            }
        } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
            const Variable<Vector>& variable = KratosComponents<Variable<Vector>>().Get(variable_name);
            if (value.IsVector()) {
                property.SetValue(variable, value.GetVector());
            } else if (value.IsMatrix()) {
                Vector temp;
                const Matrix& value_variable = value.GetMatrix();
                for (std::size_t index = 0; index < value_variable.size1(); index++)
                    temp[index] = value_variable(index, 0);
                property.SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
            const Variable<Matrix>& variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
            if (value.IsMatrix()) {
                property.SetValue(variable, value.GetMatrix());
            } else if (value.IsVector()) {
                Matrix temp;
                const Vector& value_variable = value.GetVector();
                for (std::size_t index = 0; index < value_variable.size(); index++)
                    temp(index, 0) = value_variable[index];
                property.SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
            const Variable<std::string>& variable = KratosComponents<Variable<std::string>>().Get(variable_name);
            if (value.IsString()) {
                property.SetValue(variable, value.GetString());
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else {
            KRATOS_ERROR << "Value type not defined";
        }
    }

    // Add / override tables in the p_properties
    Kratos::Parameters tables = Data["Material"]["Tables"];
    for(auto iter = tables.begin(); iter != tables.end(); iter++) {
        auto table_param = tables.GetValue(iter.name());
        // Case table is double, double. TODO(marandra): Does it make sense to consider other cases?
        Table<double> table;

        std::string input_var_name = table_param["input_variable"].GetString();
        TrimComponentName(input_var_name);
        std::string output_var_name = table_param["output_variable"].GetString();
        TrimComponentName(output_var_name);

        const auto input_var = KratosComponents<Variable<double>>().Get(input_var_name);
        const auto output_var = KratosComponents<Variable<double>>().Get(output_var_name);
        for (auto i = 0; i < table_param["data"].size(); i++) {
            table.insert(table_param["data"][i][0].GetDouble(),
                         table_param["data"][i][1].GetDouble());
        }
        property.SetTable(input_var, output_var, table);
    }
    mProperties_map[property_id] = property;
}
/***********************************************************************************/
/***********************************************************************************/
void RVELaw::InitializeMaterial(const Properties& rUnusedProperties,
                                const GeometryType& rUnusedElementGeometry,
                                const Vector& rUnusedShapeFunctionsValues)
{
    for (auto i = 0; i < mCL_vec.size(); i++)
    {
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        // Passing empty arguments, as individual CLs don't use them.
        const GeometryType dummy_element_geometry;
        const Vector dummy_shape_functions_value;

        mCL_vec[i]->InitializeMaterial(material_props,
                                       dummy_element_geometry,
                                       dummy_shape_functions_value);
    }
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::FinalizeSolutionStep(const Properties& rUnusedProperties,
                                  const GeometryType& rUnusedElementGeometry,
                                  const Vector& rUnusedShapeFunctionsValues,
                                  const ProcessInfo& rCurrentProcessInfo)
{
    for (auto i = 0; i < mCL_vec.size(); i++)
    {
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        // Passing empty arguments, as individual CLs don't use them.
        const GeometryType dummy_element_geometry;
        const Vector dummy_shape_functions_value;

        mCL_vec[i]->FinalizeSolutionStep(material_props,
                                         dummy_element_geometry,
                                         dummy_shape_functions_value,
                                         rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_modes = mB_vec[0].size2();
    const std::size_t nr_comps = GetStrainSize();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    Vector& homog_stress = rValues.GetStressVector(); // output
    homog_stress.clear();
    Matrix& homog_C = rValues.GetConstitutiveMatrix(); // output
    homog_C.clear();

    Matrix A(nr_modes, nr_modes);
    Vector res(nr_modes);
    Vector Dx(nr_modes);

    Accumulate(A, res, strain_macro, process_info);
    double residual = norm_2(res);
    double current_residual = residual;
    double ratio = 1.0;
    std::size_t it = 1;

    while (residual > mAbsoluteTolerance and ratio > mRelativeTolerance and it < mMaxIteration)
    {
        Solve(A, res, Dx);
        mModesWeights -= Dx;
        Accumulate(A, res, strain_macro, process_info);
        KRATOS_INFO_IF("RVE Law", mVerbose) << "Iteration " << it << " Residual: " << residual
                               << " Relative:" << ratio <<std::endl;
        current_residual = norm_2(res);
        ratio = current_residual / residual;
        residual = current_residual;
        it++;
    }
    KRATOS_INFO_IF("RVE Law", mVerbose) << std::endl;

    // Homogenize stress and constitutive tensor
    Matrix homog_C_taylor = ZeroMatrix(nr_comps, nr_comps);
    Matrix homog_C_fluct = ZeroMatrix(nr_comps, nr_comps);
    Matrix homog_C_fluct_aux = ZeroMatrix(nr_comps, nr_modes);
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
        CalculateIndividualMaterialResponse(stress, constit, strain, process_info, i);
        homog_stress += mIW_vec[i] * stress;
        homog_C_taylor += mIW_vec[i] * constit;
        homog_Q += mIW_vec[i] * prod(trans(mB_vec[i]), constit);
        vol_rve += mIW_vec[i];
    }
    homog_stress /= vol_rve;
    homog_Op = - prod(invA, homog_Q);
    for (auto i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        // TODO(marcelo): strain argument should be const
        CalculateIndividualMaterialResponse(stress, constit, strain, process_info, i);
        homog_C_fluct_aux += mIW_vec[i] * prod(constit, mB_vec[i]);
    }
    noalias(homog_C_fluct) = prod(homog_C_fluct_aux, homog_Op);
    homog_C = homog_C_taylor + homog_C_fluct;
    homog_C /= vol_rve;
}

void RVELaw::Solve(const Matrix &A, const Vector &res, Vector &Dx)
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

/***********************************************************************************/
/***********************************************************************************/
    void RVELaw::Accumulate(Matrix &A, Vector &res, const Vector &strain_macro, const ProcessInfo &process_info)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_modes = mB_vec[0].size2();
    const std::size_t nr_comps = GetStrainSize();
    Matrix Aux1(nr_comps, nr_modes);

    A.clear();
    res.clear();
    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);  // output
        Matrix constit(nr_comps, nr_comps);  // output
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        // TODO(marcelo): strain should be const
        CalculateIndividualMaterialResponse(stress, constit, strain, process_info, i);
        // TODO(marcelo): explicitly write triple product for A
        // Dij = BTij Ckl Blj = for k for l for j for i
        noalias(Aux1) = prod(constit, mB_vec[i]);
        noalias(A) += mIW_vec[i] * prod(trans(mB_vec[i]), Aux1);
        noalias(res) += mIW_vec[i] * prod(trans(mB_vec[i]), stress);
    }
}

/***********************************************************************************/
/***********************************************************************************/
    void RVELaw::CalculateIndividualMaterialResponse(Vector &stress,
                                                     Matrix &constit,
                                                     Vector &strain,
                                                     const ProcessInfo &process_info,
                                                     std::size_t ip_index)
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

    ConstitutiveLaw::Parameters cl_params;
    cl_params.SetOptions(cl_flags);
    cl_params.SetDeformationGradientF(F);
    cl_params.SetDeterminantF(detF);
    cl_params.SetStrainVector(strain);
    cl_params.SetStressVector(stress);
    cl_params.SetConstitutiveMatrix(constit);
    cl_params.SetShapeFunctionsValues(N);
    cl_params.SetShapeFunctionsDerivatives(DN_DX);
    const Properties material_props = mProperties_map[mPropId_vec[ip_index]];
    cl_params.SetMaterialProperties(material_props);
    cl_params.SetProcessInfo(process_info);
    // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
    // cl_params.SetElementGeometry();

    mCL_vec[ip_index]->CalculateMaterialResponseCauchy(cl_params);
}

/***********************************************************************************/
/***********************************************************************************/
bool RVELaw::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == REDUCED_MODES_WEIGHTS)
        return true;
    if (rThisVariable == INTERNAL_VARIABLES)
        return true;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/
Vector& RVELaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == REDUCED_MODES_WEIGHTS)
        rValue = mModesWeights;
    if (rThisVariable == INTERNAL_VARIABLES)
    {
        int count = 0;
        for (int i = 0; i < mCL_vec.size(); i++)
        {
            Vector rValue_i;
            mCL_vec[i]->GetValue(INTERNAL_VARIABLES, rValue_i);
            rValue.resize(count+rValue_i.size(), true);
            for (int j = 0; j < rValue_i.size(); j++) {
                rValue[count++] = rValue_i[j];
            }
        }
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::SetValue(
        const Variable<Vector>& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == INTERNAL_VARIABLES)
    {
        int count = 0;
        for (auto i = 0; i < mCL_vec.size(); i++)
        {
            int rsize = mCL_vec[i]->GetValue((const Variable<int>&)INTERNAL_VARIABLES, rsize);
            Vector rValue_i(rsize);
            for (auto j = 0; j < rsize; j++)
                rValue_i(j) = rValue(count++);
            mCL_vec[i]->SetValue(rThisVariable, rValue_i, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
int RVELaw::Check(const Properties& rUnusedProperties,
                  const GeometryType& rUnusedElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo)
{
    // Self check
    const size_t nr_comps = mB_vec.begin()->size1();
    KRATOS_ERROR_IF_NOT(nr_comps == GetStrainSize())
        << "Number of reduced base components ("<< nr_comps
        << ") differs from strain size " << GetStrainSize() << std::endl;

    // Individual CLs check
    const auto nr_points = mB_vec.size();
    for (auto i = 0; i < nr_points; i++)
    {
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        const GeometryType dummy_element_geometry;
        mCL_vec[i]->Check(material_props,
                          dummy_element_geometry,
                          rCurrentProcessInfo);
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
    rSerializer.save("mProperties_map", mProperties_map);
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
    rSerializer.load("mProperties_map", mProperties_map);
    rSerializer.load("mB_vec", mB_vec);
    rSerializer.load("mIW_vec", mIW_vec);
    rSerializer.load("mCL_vec", mCL_vec);
    rSerializer.load("mPropId_vec", mPropId_vec);
    rSerializer.load("mModesWeights", mModesWeights);
}

//************************************************************************************
//************************************************************************************

void RVELaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
    CalculateMaterialResponseCauchy(rValues);
}
void RVELaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    CalculateMaterialResponseCauchy(rValues);
}
void RVELaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) {
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void RVELaw::InitializeMaterialResponseCauchy(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
}

//************************************************************************************
//************************************************************************************

void RVELaw::InitializeMaterialResponsePK2(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void RVELaw::InitializeMaterialResponsePK1(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void RVELaw::InitializeMaterialResponseKirchhoff(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************
void RVELaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) { }
void RVELaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) { }
void RVELaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) { }
void RVELaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) { }

} /* namespace Kratos.*/
