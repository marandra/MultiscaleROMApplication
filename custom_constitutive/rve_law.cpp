#include "rve_law.h"
#include "custom_utilities/qr_utility.h"
#include "includes/checks.h"
#include "multiscale_rom_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/read_materials_utility.h"
#include "utilities/math_utils.h"

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
    const Kratos::Parameters default_parameters(R"(
    {
        "name": "constitutive law name",
        "Parameters" : {
            "rve_data_filename": "undefined_rve_data_file",
            "modified_material": [],
            "rotation": {
                "axis": [0, 0, 0],
                "angle_degree": 0.0
            },
            "convergence_criterion": "residual_criterion",
            "residual_relative_tolerance": 1e-4,
            "residual_absolute_tolerance": 1e-9,
            "max_iteration": 10,
            "verbose": 1
        }
    })");

    Params.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Read json string in file, create parameters
    const Kratos::Parameters data_params(
        ReadFile(Params["Parameters"]["rve_data_filename"].GetString()));
    mRelativeTolerance = Params["Parameters"]["residual_relative_tolerance"].GetDouble();
    mAbsoluteTolerance = Params["Parameters"]["residual_absolute_tolerance"].GetDouble();
    mMaxIteration = Params["Parameters"]["max_iteration"].GetInt();
    mVerbose = Params["Parameters"]["verbose"].GetInt();

    // Read rotation
    const double ax = Params["Parameters"]["rotation"]["axis"][0].GetDouble();
    const double ay = Params["Parameters"]["rotation"]["axis"][1].GetDouble();
    const double az  = Params["Parameters"]["rotation"]["axis"][2].GetDouble();
    const double angle = Params["Parameters"]["rotation"]["angle_degree"].GetDouble()
        * 3.141592653589793 / 180;
    mQ = Quaternion<double>::FromAxisAngle(ax, ay, az, angle);
    mQ.normalize();
    KRATOS_INFO("RVELaw") << "RVE rotation: " << angle
        << " rad on [" << ax << ", " << ay << ", " << az << "]"
        << " (q = " << mQ.W() << " + " << mQ.X() << "i + " << mQ.Y() << "j + " << mQ.Z() << "k)"
        << std::endl;
    //BoundedMatrix<double, 3, 3> Rq;
    //BoundedMatrix<double, 6, 6> Rm;
    //mQ.ToRotationMatrix(Rq);
    //const double l1 = Rq(0,0);
    //const double l2 = Rq(0,1);
    //const double l3 = Rq(0,2);
    //const double m1 = Rq(1,0);
    //const double m2 = Rq(1,1);
    //const double m3 = Rq(1,2);
    //const double n1 = Rq(2,0);
    //const double n2 = Rq(2,1);
    //const double n3 = Rq(2,2);
    // Rotation matrix (stress)
    //R[0,0]=l1*l1; R[0,1]=l2*l2; R[0,2]=l3*l3; R[0,3]=2*l2*l3;       R[0,4]=2*l1*l3;       R[0,5]=2*l1*l2;
    //R[1,0]=m1*m1; R[1,1]=m2*m2; R[1,2]=m3*m3; R[1,3]=2*m2*m3;       R[1,4]=2*m1*m3;       R[1,5]=2*m1*m2;
    //R[2,0]=n1*n1; R[2,1]=n2*n2; R[2,2]=n3*n3; R[2,3]=2*n2*n3;       R[2,4]=2*n1*n3;       R[2,5]=2*n1*n2;
    //R[3,0]=m1*n1; R[3,1]=m2*n2; R[3,2]=m3*n3; R[3,3]=(m2*n3+m3*n2); R[3,4]=(m1*n3+m3*n1); R[3,5]=(m1*n2+m2*n1);
    //R[4,0]=l1*n1; R[4,1]=l2*n2; R[4,2]=l3*n3; R[4,3]=(l2*n3+l3*n2); R[4,4]=(l1*n3+l3*n1); R[4,5]=(l1*n2+l2*n1);
    //R[5,0]=l1*m1; R[5,1]=l2*m2; R[5,2]=l3*m3; R[5,3]=(l2*m3+l3*m2); R[5,4]=(l1*m3+l3*m1); R[5,5]=(l1*m2+l2*m1);
    // Modified rotation matrix (strain)
    //Rm(0,0)=l1*l1;   Rm(0,1)=l2*l2;   Rm(0,2)=l3*l3;   Rm(0,3)=l2*l3;         Rm(0,4)=l1*l3;         Rm(0,5)=l1*l2;
    //Rm(1,0)=m1*m1;   Rm(1,1)=m2*m2;   Rm(1,2)=m3*m3;   Rm(1,3)=m2*m3;         Rm(1,4)=m1*m3;         Rm(1,5)=m1*m2;
    //Rm(2,0)=n1*n1;   Rm(2,1)=n2*n2;   Rm(2,2)=n3*n3;   Rm(2,3)=n2*n3;         Rm(2,4)=n1*n3;         Rm(2,5)=n1*n2;
    //Rm(3,0)=2*m1*n1; Rm(3,1)=2*m2*n2; Rm(3,2)=2*m3*n3; Rm(3,3)=(m2*n3+m3*n2); Rm(3,4)=(m1*n3+m3*n1); Rm(3,5)=(m1*n2+m2*n1);
    //Rm(4,0)=2*l1*n1; Rm(4,1)=2*l2*n2; Rm(4,2)=2*l3*n3; Rm(4,3)=(l2*n3+l3*n2); Rm(4,4)=(l1*n3+l3*n1); Rm(4,5)=(l1*n2+l2*n1);
    //Rm(5,0)=2*l1*m1; Rm(5,1)=2*l2*m2; Rm(5,2)=2*l3*m3; Rm(5,3)=(l2*m3+l3*m2); Rm(5,4)=(l1*m3+l3*m1); Rm(5,5)=(l1*m2+l2*m1);
    //KRATOS_INFO("RVELaw") << "RVE rotation matrix second: " << Rq << std::endl;
    //KRATOS_INFO("RVELaw") << "RVE rotation matrix fourth: " << Rm << std::endl;

    // Read material parameters:
    // material parameters are read from rve data file
    Kratos::Parameters material_parameters = data_params["material_parameters"];
    const Kratos::Parameters properties = material_parameters["properties"];
    // material parameters are modified if explicitly passed by user
    if (Params["Parameters"]["modified_material"].size() != 0)
    {
        Kratos::Parameters modified_material =
            Params["Parameters"]["modified_material"];
        for (std::size_t om = 0; om < properties.size(); ++om)
        {
            for (std::size_t mm = 0; mm < modified_material.size(); ++mm)
            {
                std::size_t op = properties.GetArrayItem(om)["properties_id"].GetInt();
                std::size_t mp =
                    modified_material.GetArrayItem(mm)["properties_id"].GetInt();
                if (op == mp)
                {
                    properties.GetArrayItem(om)["Material"].SetValue(
                        "Variables", modified_material.GetArrayItem(
                                         mm)["Material"]["Variables"]);
                    KRATOS_WARNING("RVE Law") << "WARNING: Material property " << op
                                              << " modified by user" << std::endl;
                    break;
                }
            }
        }
    }
    material_parameters.SetValue("properties", properties);

    // Create material properties
    std::string mp_name = properties.GetArrayItem(0)["model_part_name"].GetString();
    LSplit(mp_name);
    KRATOS_INFO("RVELaw") << "Modelpart name: " << mp_name << std::endl;
    Model aux_model;
    ModelPart& aux_modelpart = aux_model.CreateModelPart(mp_name);
    for (std::size_t om = 0; om < properties.size(); ++om)
    {
        std::string mp_name = properties.GetArrayItem(om)["model_part_name"].GetString();
        RSplit(mp_name);
        KRATOS_INFO("RVELaw") << "Submodelpart name: " << mp_name << std::endl;
        aux_modelpart.CreateSubModelPart(mp_name);
    }
    ReadMaterialsUtility(material_parameters.WriteJsonString(), aux_model);
    for (std::size_t om = 0; om < properties.size(); ++om)
    {
        const std::string mp_name = properties.GetArrayItem(om)["model_part_name"].GetString();
        const std::size_t property_id = properties.GetArrayItem(om)["properties_id"].GetInt();
        const ModelPart& r_aux_modelpart = aux_model.GetModelPart(mp_name);
        mProperties_map[property_id] = r_aux_modelpart.GetProperties(property_id);
    }

    // Parse data parameters
    Kratos::Parameters B_list = data_params["ip_strain_modes"];
    Kratos::Parameters w_list = data_params["ip_weight"];
    Kratos::Parameters prop_id_list = data_params["ip_property_id"];
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
        const Properties prop = mProperties_map[prop_id_list[i].GetInt()];
        const ConstitutiveLaw::Pointer pcl = prop.GetValue(CONSTITUTIVE_LAW)->Clone();
        mCL_vec.push_back(pcl);
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
               double abs_tol,
               double rel_tol,
               int max_iter,
               int verbose,
               Quaternion<double> q_rotation
               )
    : mProperties_map(pProperties),
      mB_vec(B_list),
      mIW_vec(IW_list),
      mCL_vec(CL_list),
      mPropId_vec(prop_id_list),
      mAbsoluteTolerance(abs_tol),
      mRelativeTolerance(rel_tol),
      mMaxIteration(max_iter),
      mVerbose(verbose),
      mQ(q_rotation)
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
    return Kratos::make_shared<RVELaw>(mProperties_map, mB_vec, mIW_vec,
                                       mCL_vec, mPropId_vec, mAbsoluteTolerance,
                                       mRelativeTolerance, mMaxIteration, mVerbose,
                                       mQ
                                       );
}

/***********************************************************************************/
// Copy
/***********************************************************************************/
RVELaw::RVELaw(const RVELaw& rOther) : ConstitutiveLaw(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/
//
//bool RVELaw::Has(const Variable<bool>& rThisVariable)
//{
//    if (rThisVariable == IS_INELASTIC)
//        // Here we should return "false", so the element
//        // know which function to use.
//        // GetValue when "true", CalculateValues when "false".
//        return false;
//    return false;
//}
//

//***********************************************************************************/
//***********************************************************************************/

bool RVELaw::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool RVELaw::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == REDUCED_MODES_WEIGHTS){
        return true;
    }

    if (rThisVariable == INTERNAL_VARIABLES){
        return true;
    }

    if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    // TODO: below measures are intercepted by BaseSolid element

    //if(rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR){
    //    // explicitly returning "false", so the element calls CalculateValue(...)
    //    return false;
    //}

    //if(rThisVariable == ALMANSI_STRAIN_VECTOR){
    //    // explicitly returning "false", so the element calls CalculateValue(...)
    //    return false;
    //}
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& RVELaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == REDUCED_MODES_WEIGHTS){
        // TODO: Check if it is missing a rValue.size()
        rValue = mModesWeights;
    }

    if (rThisVariable == INTERNAL_VARIABLES){
        std::size_t count = 0;
        for (std::size_t i = 0; i < mCL_vec.size(); i++)
        {
            Vector rValue_i;
            mCL_vec[i]->GetValue(INTERNAL_VARIABLES, rValue_i);
            rValue.resize(count + rValue_i.size(), true);
            for (std::size_t j = 0; j < rValue_i.size(); j++) {
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
    const ProcessInfo& rProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES){
        std::size_t count = 0;
        for (std::size_t i = 0; i < mCL_vec.size(); i++){
            Vector rValue_i;
            // call to have the right size of rValue_i
            mCL_vec[i]->GetValue(INTERNAL_VARIABLES, rValue_i);
            for (std::size_t j = 0; j < rValue_i.size(); j++)
                rValue_i(j) = rValue(count++);
            mCL_vec[i]->SetValue(INTERNAL_VARIABLES, rValue_i, rProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
void RVELaw::LSplit(std::string& rLine)
{
    std::stringstream ss(rLine);
    std::getline(ss, rLine, '.');
}

/***********************************************************************************/
/***********************************************************************************/

void RVELaw::RSplit(std::string& rLine)
{
    std::stringstream ss(rLine);
    std::size_t counter = 0;
    while (std::getline(ss, rLine, '.')){++counter;}
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

void RVELaw::InitializeMaterial(const Properties& rUnusedProperties,
                                const GeometryType& rUnusedElementGeometry,
                                const Vector& rUnusedShapeFunctionsValues)
{
    for (std::size_t i = 0; i < mCL_vec.size(); i++)
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

void RVELaw::CalculateStressResponse(Kratos::ConstitutiveLaw::Parameters &rValues,
                                     Kratos::Vector &rInternalVariables)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_modes = mB_vec[0].size2();
    const std::size_t nr_comps = GetStrainSize();
    const ProcessInfo& process_info = rValues.GetProcessInfo();
    Vector& strain_macro = rValues.GetStrainVector(); // input

    ///////////////////////////
    // Apply RVE rotation to second order tensor
    // We assume 3 and 3x3 components
    BoundedMatrix<double, 3, 3> aux_matrix;
    aux_matrix = MathUtils<double>::StrainVectorToTensor(strain_macro);

    // Test: compare rotations by quaternion and rotation matrix
    //BoundedMatrix<double, 3, 3> check_matrix;
    //BoundedMatrix<double, 3, 3> rot_matrix;
    //mQ.ToRotationMatrix(rot_matrix);
    //MathUtils<double>::BDBtProductOperation(check_matrix, aux_matrix, rot_matrix);
    //Vector check_vector(6);
    //check_vector = MathUtils<double>::StrainTensorToVector(check_matrix);
    //KRATOS_WATCH(check_vector);
    // End test

    // First, rotate columns
    for (std::size_t j = 0; j < 3; j++)
    {
        Vector aux_vector(3);
        for (std::size_t i = 0; i < 3; i++)
        {
            aux_vector[i] = aux_matrix(i, j);
        }
        mQ.RotateVector3(aux_vector);
        for (std::size_t i = 0; i < 3; i++)
        {
            aux_matrix(i, j) = aux_vector[i];
        }
    }
    // Second, rotate rows
    for (std::size_t i = 0; i < 3; i++)
    {
        Vector aux_vector(3);
        for (std::size_t j = 0; j < 3; j++)
        {
            aux_vector[j] = aux_matrix(i, j);
        }
        mQ.RotateVector3(aux_vector);
        for (std::size_t j = 0; j < 3; j++)
        {
            aux_matrix(i, j) = aux_vector[j];
        }
    }
    strain_macro = MathUtils<double>::StrainTensorToVector(aux_matrix);
    // End rotation strain

    Vector& homog_stress = rValues.GetStressVector(); // output
    homog_stress.clear();
    Matrix& homog_C = rValues.GetConstitutiveMatrix(); // output
    homog_C.clear();

    Matrix A(nr_modes, nr_modes);
    Vector res(nr_modes);
    Vector Dx(nr_modes);

    Accumulate(A, res, strain_macro, process_info);

    // Current criteria: relative displacement
    double ratio = 1.0;
    std::size_t it = 0;

    while (ratio > mRelativeTolerance and it < mMaxIteration)
    {
        Solve(A, res, Dx);
        mModesWeights -= Dx;
        Accumulate(A, res, strain_macro, process_info);
        KRATOS_INFO_IF("RVE Law", mVerbose) << "Iteration " << it
                                            << " Relative:" << ratio <<std::endl;
        const double norm_modes_weights = norm_2(mModesWeights);
        ratio = norm_2(Dx) / norm_modes_weights;
        it++;
    }
    // Previous criteria using residual.
    //double residual = norm_2(res);
    //double current_residual = residual;
    //double ratio = 1.0;
    //std::size_t it = 1;

    //while (residual > mAbsoluteTolerance and ratio > mRelativeTolerance and it < mMaxIteration)
    //{
    //    Solve(A, res, Dx);
    //    mModesWeights -= Dx;
    //    Accumulate(A, res, strain_macro, process_info);
    //    KRATOS_INFO_IF("RVE Law", mVerbose) << "Iteration " << it << " Residual: " << residual
    //                           << " Relative:" << ratio <<std::endl;
    //    current_residual = norm_2(res);
    //    ratio = current_residual / residual;
    //    residual = current_residual;
    //    it++;
    //}
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
    //for (std::size_t i = 0; i < nr_points; i++)
    //{
        //Vector stress(nr_comps);
        //Matrix constit(nr_comps, nr_comps);
        //Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        //CalculateIndividualMaterialResponse(stress, constit, strain, process_info, i);
        //homog_stress += mIW_vec[i] * stress;
        //homog_C_taylor += mIW_vec[i] * constit;
        //homog_Q += mIW_vec[i] * prod(trans(mB_vec[i]), constit);
        //vol_rve += mIW_vec[i];
    //}
    //homog_stress /= vol_rve;
    //homog_Op = - prod(invA, homog_Q);
    std::size_t count = 0;
    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
        Vector internal_variables;
        mCL_vec[i]->GetValue(INTERNAL_VARIABLES, internal_variables);
        // TODO(marcelo): strain argument should be const
        CalculateIndividualStressResponse(stress, constit, strain, internal_variables, process_info, i);
        for (std::size_t j = 0; j < internal_variables.size(); ++j)
        {
            rInternalVariables[count++] = internal_variables[j];
        }

        homog_stress += mIW_vec[i] * stress;
        homog_C_taylor += mIW_vec[i] * constit;
        homog_Q += mIW_vec[i] * prod(trans(mB_vec[i]), constit);
        vol_rve += mIW_vec[i];
        homog_C_fluct_aux += mIW_vec[i] * prod(constit, mB_vec[i]);
    }
    homog_stress /= vol_rve;
    homog_Op = - prod(invA, homog_Q);
    noalias(homog_C_fluct) = prod(homog_C_fluct_aux, homog_Op);
    homog_C = homog_C_taylor + homog_C_fluct;
    homog_C /= vol_rve;

    //////////////////////////////////
    // Apply inverse RVE rotation to stress and C
    Quaternion<double> iq = mQ.conjugate();
    aux_matrix = MathUtils<double>::StressVectorToTensor(homog_stress);

    // Test: compare rotations by quaternion and rotation matrix
    //BoundedMatrix<double, 3, 3> check_matrix;
    //BoundedMatrix<double, 3, 3> rot_matrix;
    //iq.ToRotationMatrix(rot_matrix);
    //MathUtils<double>::BDBtProductOperation(check_matrix, aux_matrix, rot_matrix);
    //Vector check_vector(6);
    //check_vector = MathUtils<double>::StressTensorToVector(check_matrix);
    //KRATOS_WATCH(check_vector);
    // End test

    // First, rotate columns
    for (std::size_t j = 0; j < 3; j++)
    {
        Vector aux_vector(3);
        for (std::size_t i = 0; i < 3; i++)
        {
            aux_vector[i] = aux_matrix(i, j);
        }
        iq.RotateVector3(aux_vector);
        for (std::size_t i = 0; i < 3; i++)
        {
            aux_matrix(i, j) = aux_vector[i];
        }
    }
    // Second, rotate rows
    for (std::size_t i = 0; i < 3; i++)
    {
        Vector aux_vector(3);
        for (std::size_t j = 0; j < 3; j++)
        {
            aux_vector[j] = aux_matrix(i, j);
        }
        iq.RotateVector3(aux_vector);
        for (std::size_t j = 0; j < 3; j++)
        {
            aux_matrix(i, j) = aux_vector[j];
        }
    }
    homog_stress = MathUtils<double>::StressTensorToVector(aux_matrix);

    // Rotate C
    BoundedMatrix<double, 3, 3> iQ;
    BoundedMatrix<double, 3, 3> Q;
    BoundedMatrix<double, 6, 6> aux_C;
    BoundedMatrix<double, 6, 6> iR;
    BoundedMatrix<double, 6, 6> Rm;
    double l1, l2, l3, m1, m2, m3, n1, n2, n3;

    // inverse Rotation matrix (stress)
    mQ.conjugate().ToRotationMatrix(iQ);
    l1 = iQ(0,0); l2 = iQ(0,1); l3 = iQ(0,2);
    m1 = iQ(1,0); m2 = iQ(1,1); m3 = iQ(1,2);
    n1 = iQ(2,0); n2 = iQ(2,1); n3 = iQ(2,2);
    iR(0,0)=l1*l1; iR(0,1)=l2*l2; iR(0,2)=l3*l3; iR(0,3)=2*l2*l3;       iR(0,4)=2*l1*l3;       iR(0,5)=2*l1*l2;
    iR(1,0)=m1*m1; iR(1,1)=m2*m2; iR(1,2)=m3*m3; iR(1,3)=2*m2*m3;       iR(1,4)=2*m1*m3;       iR(1,5)=2*m1*m2;
    iR(2,0)=n1*n1; iR(2,1)=n2*n2; iR(2,2)=n3*n3; iR(2,3)=2*n2*n3;       iR(2,4)=2*n1*n3;       iR(2,5)=2*n1*n2;
    iR(3,0)=m1*n1; iR(3,1)=m2*n2; iR(3,2)=m3*n3; iR(3,3)=(m2*n3+m3*n2); iR(3,4)=(m1*n3+m3*n1); iR(3,5)=(m1*n2+m2*n1);
    iR(4,0)=l1*n1; iR(4,1)=l2*n2; iR(4,2)=l3*n3; iR(4,3)=(l2*n3+l3*n2); iR(4,4)=(l1*n3+l3*n1); iR(4,5)=(l1*n2+l2*n1);
    iR(5,0)=l1*m1; iR(5,1)=l2*m2; iR(5,2)=l3*m3; iR(5,3)=(l2*m3+l3*m2); iR(5,4)=(l1*m3+l3*m1); iR(5,5)=(l1*m2+l2*m1);

    // modified rotation matrix (strain voigt)
    mQ.ToRotationMatrix(Q);
    l1 = Q(0,0); l2 = Q(0,1); l3 = Q(0,2);
    m1 = Q(1,0); m2 = Q(1,1); m3 = Q(1,2);
    n1 = Q(2,0); n2 = Q(2,1); n3 = Q(2,2);
    Rm(0,0)=l1*l1;   Rm(0,1)=l2*l2;   Rm(0,2)=l3*l3;   Rm(0,3)=l2*l3;         Rm(0,4)=l1*l3;         Rm(0,5)=l1*l2;
    Rm(1,0)=m1*m1;   Rm(1,1)=m2*m2;   Rm(1,2)=m3*m3;   Rm(1,3)=m2*m3;         Rm(1,4)=m1*m3;         Rm(1,5)=m1*m2;
    Rm(2,0)=n1*n1;   Rm(2,1)=n2*n2;   Rm(2,2)=n3*n3;   Rm(2,3)=n2*n3;         Rm(2,4)=n1*n3;         Rm(2,5)=n1*n2;
    Rm(3,0)=2*m1*n1; Rm(3,1)=2*m2*n2; Rm(3,2)=2*m3*n3; Rm(3,3)=(m2*n3+m3*n2); Rm(3,4)=(m1*n3+m3*n1); Rm(3,5)=(m1*n2+m2*n1);
    Rm(4,0)=2*l1*n1; Rm(4,1)=2*l2*n2; Rm(4,2)=2*l3*n3; Rm(4,3)=(l2*n3+l3*n2); Rm(4,4)=(l1*n3+l3*n1); Rm(4,5)=(l1*n2+l2*n1);
    Rm(5,0)=2*l1*m1; Rm(5,1)=2*l2*m2; Rm(5,2)=2*l3*m3; Rm(5,3)=(l2*m3+l3*m2); Rm(5,4)=(l1*m3+l3*m1); Rm(5,5)=(l1*m2+l2*m1);

    noalias(aux_C) = prod(homog_C, Rm);
    noalias(homog_C) = prod(iR, aux_C);

    // End rotation
    //////////////////////////////////
}

/***********************************************************************************/
/***********************************************************************************/

void RVELaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_modes = mB_vec[0].size2();
    const std::size_t nr_comps = GetStrainSize();
    const ProcessInfo& process_info = rValues.GetProcessInfo();
    Vector& strain_macro = rValues.GetStrainVector(); // input

    // Prepare rotation matices
    BoundedMatrix<double, 3, 3> iQ;
    BoundedMatrix<double, 3, 3> Q;
    BoundedMatrix<double, 6, 6> aux_C;
    BoundedMatrix<double, 6, 6> iR;
    BoundedMatrix<double, 6, 6> Rm;
    double l1, l2, l3, m1, m2, m3, n1, n2, n3;
    // inverse Rotation matrix (stress)
    mQ.conjugate().ToRotationMatrix(iQ);
    l1 = iQ(0,0); l2 = iQ(0,1); l3 = iQ(0,2);
    m1 = iQ(1,0); m2 = iQ(1,1); m3 = iQ(1,2);
    n1 = iQ(2,0); n2 = iQ(2,1); n3 = iQ(2,2);
    iR(0,0)=l1*l1; iR(0,1)=l2*l2; iR(0,2)=l3*l3; iR(0,3)=2*l2*l3;       iR(0,4)=2*l1*l3;       iR(0,5)=2*l1*l2;
    iR(1,0)=m1*m1; iR(1,1)=m2*m2; iR(1,2)=m3*m3; iR(1,3)=2*m2*m3;       iR(1,4)=2*m1*m3;       iR(1,5)=2*m1*m2;
    iR(2,0)=n1*n1; iR(2,1)=n2*n2; iR(2,2)=n3*n3; iR(2,3)=2*n2*n3;       iR(2,4)=2*n1*n3;       iR(2,5)=2*n1*n2;
    iR(3,0)=m1*n1; iR(3,1)=m2*n2; iR(3,2)=m3*n3; iR(3,3)=(m2*n3+m3*n2); iR(3,4)=(m1*n3+m3*n1); iR(3,5)=(m1*n2+m2*n1);
    iR(4,0)=l1*n1; iR(4,1)=l2*n2; iR(4,2)=l3*n3; iR(4,3)=(l2*n3+l3*n2); iR(4,4)=(l1*n3+l3*n1); iR(4,5)=(l1*n2+l2*n1);
    iR(5,0)=l1*m1; iR(5,1)=l2*m2; iR(5,2)=l3*m3; iR(5,3)=(l2*m3+l3*m2); iR(5,4)=(l1*m3+l3*m1); iR(5,5)=(l1*m2+l2*m1);
    // modified rotation matrix (strain voigt)
    mQ.ToRotationMatrix(Q);
    l1 = Q(0,0); l2 = Q(0,1); l3 = Q(0,2);
    m1 = Q(1,0); m2 = Q(1,1); m3 = Q(1,2);
    n1 = Q(2,0); n2 = Q(2,1); n3 = Q(2,2);
    Rm(0,0)=l1*l1;   Rm(0,1)=l2*l2;   Rm(0,2)=l3*l3;   Rm(0,3)=l2*l3;         Rm(0,4)=l1*l3;         Rm(0,5)=l1*l2;
    Rm(1,0)=m1*m1;   Rm(1,1)=m2*m2;   Rm(1,2)=m3*m3;   Rm(1,3)=m2*m3;         Rm(1,4)=m1*m3;         Rm(1,5)=m1*m2;
    Rm(2,0)=n1*n1;   Rm(2,1)=n2*n2;   Rm(2,2)=n3*n3;   Rm(2,3)=n2*n3;         Rm(2,4)=n1*n3;         Rm(2,5)=n1*n2;
    Rm(3,0)=2*m1*n1; Rm(3,1)=2*m2*n2; Rm(3,2)=2*m3*n3; Rm(3,3)=(m2*n3+m3*n2); Rm(3,4)=(m1*n3+m3*n1); Rm(3,5)=(m1*n2+m2*n1);
    Rm(4,0)=2*l1*n1; Rm(4,1)=2*l2*n2; Rm(4,2)=2*l3*n3; Rm(4,3)=(l2*n3+l3*n2); Rm(4,4)=(l1*n3+l3*n1); Rm(4,5)=(l1*n2+l2*n1);
    Rm(5,0)=2*l1*m1; Rm(5,1)=2*l2*m2; Rm(5,2)=2*l3*m3; Rm(5,3)=(l2*m3+l3*m2); Rm(5,4)=(l1*m3+l3*m1); Rm(5,5)=(l1*m2+l2*m1);

    // rotate strain
    strain_macro = prod(Rm, strain_macro);

    Vector& homog_stress = rValues.GetStressVector(); // output
    homog_stress.clear();
    Matrix& homog_C = rValues.GetConstitutiveMatrix(); // output
    homog_C.clear();

    Matrix A(nr_modes, nr_modes);
    Vector res(nr_modes);
    Vector Dx(nr_modes);

    Accumulate(A, res, strain_macro, process_info);

    // Current criteria: relative displacement
    double ratio = 1.0;
    std::size_t it = 0;

    while (ratio > mRelativeTolerance and it < mMaxIteration)
    {
        Solve(A, res, Dx);
        mModesWeights -= Dx;
        Accumulate(A, res, strain_macro, process_info);
        KRATOS_INFO_IF("RVE Law", mVerbose) << "Iteration " << it
                                            << " Relative:" << ratio <<std::endl;
        const double norm_modes_weights = norm_2(mModesWeights);
        ratio = norm_2(Dx) / norm_modes_weights;
        it++;
    }
    // Previous criteria using residual.
    //double residual = norm_2(res);
    //double current_residual = residual;
    //double ratio = 1.0;
    //std::size_t it = 1;

    //while (residual > mAbsoluteTolerance and ratio > mRelativeTolerance and it < mMaxIteration)
    //{
    //    Solve(A, res, Dx);
    //    mModesWeights -= Dx;
    //    Accumulate(A, res, strain_macro, process_info);
    //    KRATOS_INFO_IF("RVE Law", mVerbose) << "Iteration " << it << " Residual: " << residual
    //                           << " Relative:" << ratio <<std::endl;
    //    current_residual = norm_2(res);
    //    ratio = current_residual / residual;
    //    residual = current_residual;
    //    it++;
    //}
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
    for (std::size_t i = 0; i < nr_points; i++)
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
    for (std::size_t i = 0; i < nr_points; i++)
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

    // rotate stress to original base
    homog_stress = prod(iR, homog_stress);

    // rotate C to original base
    noalias(aux_C) = prod(homog_C, Rm);
    noalias(homog_C) = prod(iR, aux_C);
}

/***********************************************************************************/
/***********************************************************************************/

void RVELaw::Solve(const Matrix &A, const Vector &res, Vector &Dx)
{
    const std::size_t nr_modes = mB_vec[0].size2();
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
    for (std::size_t ii = 0; ii < nr_modes; ii++)
    {
        for (std::size_t jj = 0; jj < nr_modes; jj++)
        {
            aux_qr_A[ii][jj] = A(ii, jj);
        }
        aux_qr_res[ii] = res(ii);
    }
    QR_decomposition.compute(nr_modes, nr_modes, &(aux_qr_A[0][0]));
    QR_decomposition.solve(&(aux_qr_res[0]), &(aux_qr_Dx[0]));

    // Update
    for (std::size_t ii = 0; ii < nr_modes; ii++)
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

    void RVELaw::CalculateIndividualStressResponse(Vector &stress,
                                                     Matrix &constit,
                                                     Vector &strain,
                                                     Vector &rInternalVariables,
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

    mCL_vec[ip_index]->CalculateStressResponse(cl_params, rInternalVariables);
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

//************************************************************************************
//************************************************************************************
//
//bool& RVELaw::CalculateValue(
//        ConstitutiveLaw::Parameters& rValues,
//        const Variable<bool>& rThisVariable,
//        bool& rValue
//)
//{
//    if (rThisVariable == IS_INELASTIC) {
//
//        const std::size_t nr_points = mB_vec.size();
//
//        const ProcessInfo& process_info = rValues.GetProcessInfo();
//        const Vector& strain_macro = rValues.GetStrainVector(); // input
//
//        rValue = false;
//        for (std::size_t i = 0; i < nr_points; i++) {
//            const Properties material_props = mProperties_map[mPropId_vec[i]];
//            ConstitutiveLaw::Parameters cl_params;
//            cl_params.SetMaterialProperties(material_props);
//
//            Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
//            cl_params.SetStrainVector(strain);
//
//            cl_params.SetProcessInfo(process_info);
//
//            // IS_INELASTIC is retrieved with CalculateValue, so Has() returns "false"
//            if (not mCL_vec[i]->Has(IS_INELASTIC)) {
//                mCL_vec[i]->CalculateValue(cl_params, IS_INELASTIC, rValue);
//                if (rValue == true)
//                    return rValue;  // leave when the first "true" is found
//            }
//        }
//    }
//
//    return rValue;
//}

//************************************************************************************
//************************************************************************************

double& RVELaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == STRAIN_ENERGY) {
        const std::size_t nr_points = mB_vec.size();
        double volume = 0.0;
        Vector dummy;
        Vector& r_strain_energy_vector = dummy;

        CalculateValue(rParametersValues, STRAIN_ENERGY_VECTOR, r_strain_energy_vector);
        rValue = 0.0;
        for (std::size_t i = 0; i < nr_points; i++)
        {
            rValue += mIW_vec[i] * r_strain_energy_vector[i];
            volume += mIW_vec[i];
        }
        rValue /= volume;
     }

    return rValue;
}

//************************************************************************************
//************************************************************************************

Vector& RVELaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN_ENERGY_VECTOR) {
        const std::size_t nr_points = mB_vec.size();

        if (rValue.size() != nr_points)
            rValue.resize(nr_points, false);
        rValue.clear();

        const ProcessInfo& process_info = rParametersValues.GetProcessInfo();
        Vector& r_strain_vector = rParametersValues.GetStrainVector(); // input
        // In case there is an initial state (i.e., in sampling)
        //AddInitialStrainVectorContribution(r_strain_vector, rParametersValues);
        if (rParametersValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rParametersValues.GetProcessInfo()[INITIAL_STRAIN];
        }

        for (std::size_t i = 0; i < nr_points; i++) {
            double dummy;
            double& strain_energy = dummy;
            const Properties material_props = mProperties_map[mPropId_vec[i]];
            ConstitutiveLaw::Parameters cl_params;
            cl_params.SetMaterialProperties(material_props);
            Vector strain = r_strain_vector + prod(mB_vec[i], mModesWeights);
            cl_params.SetStrainVector(strain);
            cl_params.SetProcessInfo(process_info);
            mCL_vec[i]->CalculateValue(cl_params, STRAIN_ENERGY, strain_energy);
            rValue[i] = strain_energy;
        }
    }

    //if (rThisVariable == STRAIN ||
    //    rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
    //    rThisVariable == ALMANSI_STRAIN_VECTOR){
    if (rThisVariable == STRAIN){
        rValue = rParametersValues.GetStrainVector();
        // In case there is an initial state (i.e., in sampling)
        //AddInitialStrainVectorContribution(r_strain_vector, rParametersValues);
        if (rParametersValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(rValue) += rParametersValues.GetProcessInfo()[INITIAL_STRAIN];
        }
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

Matrix& RVELaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CAUCHY_STRESS_VECTOR_LIST) {
        const std::size_t nr_points = mB_vec.size();
        const std::size_t nr_comps = GetStrainSize();
        const Vector& strain_macro = rParametersValues.GetStrainVector();
        const ProcessInfo& process_info = rParametersValues.GetProcessInfo();

        if (rValue.size1() != nr_points || rValue.size2() != nr_comps)
            rValue.resize(nr_points, nr_comps, false);
        rValue.clear();

        for (std::size_t i = 0; i < nr_points; i++) {
            Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
            Vector stress(nr_comps);
            Matrix constit(nr_comps, nr_comps);  // unused
            CalculateIndividualMaterialResponse(stress, constit, strain, process_info, i);
            for (std::size_t j = 0; j < nr_comps; j++)
                rValue(i, j) = stress[j];
        }
    }
    return rValue;
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
    const std::size_t nr_points = mB_vec.size();
    for (std::size_t i = 0; i < nr_points; i++)
    {
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        const GeometryType dummy_element_geometry;
        mCL_vec[i]->Check(material_props,
                          dummy_element_geometry,
                          rCurrentProcessInfo);
    }

    return 0;
}

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

void RVELaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rParametersValues) {
    CalculateMaterialResponseCauchy(rParametersValues);
}
void RVELaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rParametersValues) {
    CalculateMaterialResponseCauchy(rParametersValues);
}
void RVELaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rParametersValues) {
    CalculateMaterialResponseCauchy(rParametersValues);
}

//************************************************************************************
//************************************************************************************

void RVELaw::InitializeMaterialResponseCauchy(
        Kratos::ConstitutiveLaw::Parameters &rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->InitializeMaterialResponseCauchy(cl_params);
    }
}

void RVELaw::InitializeMaterialResponsePK2(Kratos::ConstitutiveLaw::Parameters &rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const auto dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->InitializeMaterialResponsePK2(cl_params);
    }
}

void RVELaw::InitializeMaterialResponsePK1(Kratos::ConstitutiveLaw::Parameters &rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->InitializeMaterialResponsePK1(cl_params);
    }
}

void RVELaw::InitializeMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters &rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const auto dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->InitializeMaterialResponseKirchhoff(cl_params);
    }
}

//************************************************************************************
//************************************************************************************

void RVELaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->FinalizeMaterialResponsePK1(cl_params);
    }
}

void RVELaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->FinalizeMaterialResponsePK2(cl_params);
    }
}

void RVELaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->FinalizeMaterialResponseKirchhoff(cl_params);
    }
}

void RVELaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const std::size_t nr_points = mB_vec.size();
    const std::size_t nr_comps = GetStrainSize();
    const std::size_t dim = WorkingSpaceDimension();
    const Vector& strain_macro = rValues.GetStrainVector(); // input
    const ProcessInfo& process_info = rValues.GetProcessInfo();

    for (std::size_t i = 0; i < nr_points; i++)
    {
        Vector stress(nr_comps);
        Matrix constit(nr_comps, nr_comps);
        Vector strain = strain_macro + prod(mB_vec[i], mModesWeights);
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
        const Properties material_props = mProperties_map[mPropId_vec[i]];
        cl_params.SetMaterialProperties(material_props);
        cl_params.SetProcessInfo(process_info);
        // TODO(marcelo): needs HF elem geom. Currently not used in our iCL.
        // cl_params.SetElementGeometry();

        mCL_vec[i]->FinalizeMaterialResponseCauchy(cl_params);
    }
}

} /* namespace Kratos.*/
