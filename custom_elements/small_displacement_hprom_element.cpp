#include "custom_elements/small_displacement_hprom_element.hpp"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "multiscale_rom_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(SmallDisplacementHpromElement, COMPUTE_RHS_VECTOR, 0);
KRATOS_CREATE_LOCAL_FLAG(SmallDisplacementHpromElement, COMPUTE_LHS_MATRIX, 1);
KRATOS_CREATE_LOCAL_FLAG(SmallDisplacementHpromElement, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2);
KRATOS_CREATE_LOCAL_FLAG(SmallDisplacementHpromElement, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3);

SmallDisplacementHpromElement::SmallDisplacementHpromElement() : Element()
{
}

SmallDisplacementHpromElement::SmallDisplacementHpromElement(IndexType NewId,
                                                                 GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

SmallDisplacementHpromElement::SmallDisplacementHpromElement(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

SmallDisplacementHpromElement::SmallDisplacementHpromElement(SmallDisplacementHpromElement const& rOther)
    : Element(rOther),
      mThisIntegrationMethod(rOther.mThisIntegrationMethod),
      mConstitutiveLawVector(rOther.mConstitutiveLawVector)
{
}

SmallDisplacementHpromElement& SmallDisplacementHpromElement::operator=(
    SmallDisplacementHpromElement const& rOther)
{
    Element::operator=(rOther);
    mThisIntegrationMethod = rOther.mThisIntegrationMethod;
    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize(rOther.mConstitutiveLawVector.size());
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }
    return *this;
}

Element::Pointer SmallDisplacementHpromElement::Create(IndexType NewId,
                                                         NodesArrayType const& rThisNodes,
                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallDisplacementHpromElement(
        NewId, GetGeometry().Create(rThisNodes), pProperties));
}

Element::Pointer SmallDisplacementHpromElement::Clone(IndexType NewId,
                                                        NodesArrayType const& rThisNodes) const
{
    SmallDisplacementHpromElement NewElement(
        NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if (NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
    {
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

        if (NewElement.mConstitutiveLawVector.size() !=
            NewElement.GetGeometry().IntegrationPointsNumber())
            KRATOS_THROW_ERROR(std::logic_error,
                               "constitutive law not has the correct size ",
                               NewElement.mConstitutiveLawVector.size());
    }

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new SmallDisplacementHpromElement(NewElement));
}

SmallDisplacementHpromElement::~SmallDisplacementHpromElement()
{
}

// GETTING METHODS

SmallDisplacementHpromElement::IntegrationMethod SmallDisplacementHpromElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

// TODO check if this funcion es required in this case
/*
void SmallDisplacementHpromElement::GetDofList(
    DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo) {
  rElementalDofList.resize(0);
  const size_t dimension = GetGeometry().WorkingSpaceDimension();

  for (unsigned int i = 0; i < GetGeometry().size(); i++) {
    rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
    rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
    if (dimension == 3)
      rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
  }
}
*/

void SmallDisplacementHpromElement::EquationIdVector(EquationIdVectorType& rResult,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    // TODO mNumberOfModes is not initialized yet
    const size_t number_of_modes =
        static_cast<size_t>(rCurrentProcessInfo[NUMBER_REDUCED_MODES]);

    if (rResult.size() != number_of_modes)
        rResult.resize(number_of_modes, false);

    for (std::size_t i = 0; i < number_of_modes; i++)
        rResult[i] = static_cast<unsigned long>(i);
}

void SmallDisplacementHpromElement::SetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0;
         PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        mConstitutiveLawVector[PointNumber]->SetValue(
            rVariable, rValues[PointNumber], rCurrentProcessInfo);
    }
}

void SmallDisplacementHpromElement::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0;
         PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        mConstitutiveLawVector[PointNumber]->SetValue(
            rVariable, rValues[PointNumber], rCurrentProcessInfo);
    }
}

void SmallDisplacementHpromElement::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int PointNumber = 0;
         PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        mConstitutiveLawVector[PointNumber]->SetValue(
            rVariable, rValues[PointNumber], rCurrentProcessInfo);
    }
}

void SmallDisplacementHpromElement::SetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CONSTITUTIVE_LAW)
    {
        if (mConstitutiveLawVector.size() != rValues.size())
        {
            mConstitutiveLawVector.resize(rValues.size());

            if (mConstitutiveLawVector.size() !=
                GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod))
                KRATOS_THROW_ERROR(std::logic_error,
                                   "constitutive law not has the correct size ",
                                   mConstitutiveLawVector.size());
        }

        for (unsigned int i = 0; i < rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i];
        }
    }
}

//************************************************************************************

void SmallDisplacementHpromElement::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VON_MISES_STRESS)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else
    {
        const size_t& integration_points_number =
            GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

        if (rValues.size() != integration_points_number)
            rValues.resize(integration_points_number);

        for (unsigned int ii = 0; ii < integration_points_number; ii++)
        {
            rValues[ii] = 0.0;
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rValues[ii]);
        }
    }
}

void SmallDisplacementHpromElement::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    const size_t& integration_points_number = mConstitutiveLawVector.size();

    if (rValues.size() != integration_points_number)
        rValues.resize(integration_points_number);

    if (rVariable == PK2_STRESS_TENSOR || rVariable == CAUCHY_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else if (rVariable == PK2_STRESS_VECTOR || rVariable == CAUCHY_STRESS_VECTOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR || rVariable == ALMANSI_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else
    {
        for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(
                rVariable, rValues[PointNumber]);
        }
    }
}

void SmallDisplacementHpromElement::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    const size_t& integration_points_number = mConstitutiveLawVector.size();

    if (rValues.size() != integration_points_number)
        rValues.resize(integration_points_number);

    if (rVariable == PK2_STRESS_TENSOR || rVariable == CAUCHY_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR || rVariable == ALMANSI_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    else
    {
        for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(
                rVariable, rValues[PointNumber]);
        }
    }
}

void SmallDisplacementHpromElement::GetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CONSTITUTIVE_LAW)
    {
        if (rValues.size() != mConstitutiveLawVector.size())
        {
            rValues.resize(mConstitutiveLawVector.size());
        }

        for (unsigned int i = 0; i < rValues.size(); i++)
        {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }
}

// STARTING - ENDING  METHODS

void SmallDisplacementHpromElement::Initialize()
{
    KRATOS_TRY

    mVoigtSize = 3;
    // mVoigtSize = 4; // for 2D Bbar element + plasticiy
    if (GetGeometry().WorkingSpaceDimension() == 3)
    {
        mVoigtSize = 6;
    }

    InitializeMaterial();

    KRATOS_CATCH("")
}

void SmallDisplacementHpromElement::SetGeneralVariables(GeneralVariables& rVariables,
                                                          ConstitutiveLaw::Parameters& rValues,
                                                          const size_t& rPointNumber)
{
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

    if (rVariables.detJ < 0)
    {
        KRATOS_THROW_ERROR(
            std::invalid_argument,
            " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", rVariables.detJ)
    }

    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);
}

void SmallDisplacementHpromElement::InitializeGeneralVariables(GeneralVariables& rVariables,
                                                                 const ProcessInfo& rCurrentProcessInfo)
{
    const size_t number_of_nodes = GetGeometry().size();
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    rVariables.Initialize(mVoigtSize, dimension, number_of_nodes, rCurrentProcessInfo);

    // needed parameters for consistency with the general constitutive law:
    // small
    // displacements
    rVariables.detF = 1.0;
    rVariables.F = identity_matrix<double>(dimension);

    // set variables including all integration points values

    // reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod));

    // reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(
        GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod));

    // calculating the current jacobian from cartesian coordinates to parent
    // coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian(rVariables.j, mThisIntegrationMethod);

    // Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    // calculating the reference jacobian from cartesian coordinates to parent
    // coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian(rVariables.J, mThisIntegrationMethod,
                                          rVariables.DeltaPosition);
}

void SmallDisplacementHpromElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                                               VectorType& rRightHandSideVector,
                                                               Flags& rCalculationFlags)
{
    const size_t MatSize = mNumberOfModes;

    if (rCalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX))
    {
        if (rLeftHandSideMatrix.size1() != MatSize)
            rLeftHandSideMatrix.resize(MatSize, MatSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize);
    }

    if (rCalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR))
    {
        if (rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize, false);

        noalias(rRightHandSideVector) = ZeroVector(MatSize);
    }
}

void SmallDisplacementHpromElement::CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                                               ProcessInfo& rCurrentProcessInfo)
{
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    Flags& ConstitutiveLawOptions = Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // auxiliary terms
    const size_t dimension = GetGeometry().WorkingSpaceDimension();
    Vector VolumeForce(dimension);
    noalias(VolumeForce) = ZeroVector(dimension);

    Vector& mAssignedIntegrationWeights = this->GetValue(INTEGRATION_POINT_WEIGHT);

    for (size_t point_number = 0; point_number < integration_points.size(); point_number++)
    {
        double geometrical_integration_weight;
        geometrical_integration_weight = mAssignedIntegrationWeights[point_number];

        if (geometrical_integration_weight < 0)
        {
            continue;
        }

        this->CalculateKinematics(Variables, point_number);
        this->SetGeneralVariables(Variables, Values, point_number);
        mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy(Values);

        if (rLocalSystem.CalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX))
        {
            // contributions to stiffness matrix calculated on the reference
            // config
            this->CalculateAndAddLHS(rLocalSystem, Variables, geometrical_integration_weight);
        }
        if (rLocalSystem.CalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR))
        {
            // contribution to external forces
            VolumeForce = this->CalculateVolumeForce(VolumeForce, Variables);
            this->CalculateAndAddRHS(rLocalSystem, Variables, VolumeForce,
                                     geometrical_integration_weight);
        }
    }
}

//************************************************************************************

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                                         GeneralVariables& rVariables,
                                                         double& rIntegrationWeight)
{
    // contributions of the stiffness matrix calculated on the reference
    // configuration
    if (rLocalSystem.CalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX_WITH_COMPONENTS))
    {
        std::vector<MatrixType>& rLeftHandSideMatrices =
            rLocalSystem.GetLeftHandSideMatrices();
        const std::vector<Variable<MatrixType>>& rLeftHandSideVariables =
            rLocalSystem.GetLeftHandSideVariables();

        for (unsigned int i = 0; i < rLeftHandSideVariables.size(); i++)
        {
            bool calculated = false;
            if (rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX)
            {
                // operation performed: add Km to the rLefsHandSideMatrix
                this->CalculateAndAddKuum(rLeftHandSideMatrices[i], rVariables,
                                          rIntegrationWeight);
                calculated = true;
            }

            if (!calculated)
            {
                KRATOS_THROW_ERROR(std::logic_error,
                                   " ELEMENT can not supply the required local "
                                   "system variable: ",
                                   rLeftHandSideVariables[i])
            }
        }
    }
    else
    {
        MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        // operation performed: add Km to the rLefsHandSideMatrix
        this->CalculateAndAddKuum(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                                         GeneralVariables& rVariables,
                                                         Vector& rVolumeForce,
                                                         double& rIntegrationWeight)
{
    // contribution of the internal and external forces
    if (rLocalSystem.CalculationFlags.Is(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR_WITH_COMPONENTS))
    {
        std::vector<VectorType>& rRightHandSideVectors =
            rLocalSystem.GetRightHandSideVectors();
        const std::vector<Variable<VectorType>>& rRightHandSideVariables =
            rLocalSystem.GetRightHandSideVariables();
        for (unsigned int i = 0; i < rRightHandSideVariables.size(); i++)
        {
            bool calculated = false;

            /*if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ){
              // operation performed: rRightHandSideVector +=
            ExtForce*IntToReferenceWeight
              this->CalculateAndAddExternalForces( rRightHandSideVectors[i],
            rVariables, rVolumeForce, rIntegrationWeight );
              calculated = true;
            }*/

            if (rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR)
            {
                // operation performed: rRightHandSideVector -=
                // IntForce*IntToReferenceWeight
                this->CalculateAndAddInternalForces(
                    rRightHandSideVectors[i], rVariables, rIntegrationWeight);
                calculated = true;
            }

            if (!calculated)
            {
                KRATOS_THROW_ERROR(std::logic_error,
                                   " ELEMENT can not supply the required local "
                                   "system variable: ",
                                   rRightHandSideVariables[i])
            }
        }
    }
    else
    {
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

        // operation performed: rRightHandSideVector +=
        // ExtForce*IntToReferenceWeight
        // this->CalculateAndAddExternalForces( rRightHandSideVector,
        // rVariables,
        // rVolumeForce, rIntegrationWeight );

        // operation performed: rRightHandSideVector -=
        // IntForce*IntToReferenceWeight
        this->CalculateAndAddInternalForces(rRightHandSideVector, rVariables, rIntegrationWeight);
    }
}

//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

double& SmallDisplacementHpromElement::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension == 2)
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                             ProcessInfo& rCurrentProcessInfo)
{
    // create local system components
    LocalSystemComponents LocalSystem;

    // calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices(LeftHandSideMatrix, rRightHandSideVector,
                                   LocalSystem.CalculationFlags);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateRightHandSide(
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType>>& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo)
{
    // create local system components
    LocalSystemComponents LocalSystem;

    // calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(
        SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    // Initialize sizes for the system components:
    if (rRHSVariables.size() != rRightHandSideVectors.size())
        rRightHandSideVectors.resize(rRHSVariables.size());

    for (unsigned int i = 0; i < rRightHandSideVectors.size(); i++)
    {
        this->InitializeSystemMatrices(LeftHandSideMatrix, rRightHandSideVectors[i],
                                       LocalSystem.CalculationFlags);
    }

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    // Calculate elemental system
    CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                            ProcessInfo& rCurrentProcessInfo)
{
    // create local system components
    LocalSystemComponents LocalSystem;

    // calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices(rLeftHandSideMatrix, RightHandSideVector,
                                   LocalSystem.CalculationFlags);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                           VectorType& rRightHandSideVector,
                                                           ProcessInfo& rCurrentProcessInfo)
{
    // create local system components
    LocalSystemComponents LocalSystem;

    // calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR);

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices(rLeftHandSideMatrix, rRightHandSideVector,
                                   LocalSystem.CalculationFlags);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateLocalSystem(
    std::vector<MatrixType>& rLeftHandSideMatrices,
    const std::vector<Variable<MatrixType>>& rLHSVariables,
    std::vector<VectorType>& rRightHandSideVectors,
    const std::vector<Variable<VectorType>>& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo)
{
    // create local system components
    LocalSystemComponents LocalSystem;

    // calculation flags
    LocalSystem.CalculationFlags.Set(
        SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(
        SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    // Initialize sizes for the system components:
    if (rLHSVariables.size() != rLeftHandSideMatrices.size())
        rLeftHandSideMatrices.resize(rLHSVariables.size());

    if (rRHSVariables.size() != rRightHandSideVectors.size())
        rRightHandSideVectors.resize(rRHSVariables.size());

    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX);
    for (unsigned int i = 0; i < rLeftHandSideMatrices.size(); i++)
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices(
            rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags);
    }

    LocalSystem.CalculationFlags.Set(SmallDisplacementHpromElement::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(
        SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX, false);

    for (unsigned int i = 0; i < rRightHandSideVectors.size(); i++)
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices(
            rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags);
    }
    LocalSystem.CalculationFlags.Set(
        SmallDisplacementHpromElement::COMPUTE_LHS_MATRIX, true);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    // Calculate elemental system
    CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);
}

////************************************************************************************
////************************************************************************************

void SmallDisplacementHpromElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);
    size_t row_counter;

    // TODO: initialize only once
    {
        mModesWeights.resize(static_cast<size_t>(rCurrentProcessInfo[NUMBER_REDUCED_MODES]));
        Matrix& BMatrixImported = this->GetValue(REDUCED_MODES_MATRIX);
        mNumberOfModes =
            static_cast<std::size_t>(rCurrentProcessInfo[NUMBER_REDUCED_MODES]);
        mBMatrixVector.resize(integration_points.size());
        for (unsigned int PointNumber = 0;
             PointNumber < integration_points.size(); PointNumber++)
        {
            mBMatrixVector[PointNumber].resize(mVoigtSize, mNumberOfModes);
        }
        row_counter = 0;
        // KRATOS_WATCH(this->Id())
        for (size_t PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
        {
            for (size_t voigt_component = 0; voigt_component < mVoigtSize; voigt_component++)
            {
                for (size_t mode = 0; mode < mNumberOfModes; mode++)
                {
                    mBMatrixVector[PointNumber](voigt_component, mode) =
                        BMatrixImported(row_counter, mode);
                }
                row_counter++;
            }
            //KRATOS_WATCH(mBMatrixVector[PointNumber])
        }
    }

    mModesWeights = rCurrentProcessInfo[REDUCED_MODES_WEIGHTS];
    ClearNodalForces();

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        mConstitutiveLawVector[i]->InitializeSolutionStep(
            GetProperties(), GetGeometry(),
            row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i),
            rCurrentProcessInfo);
}

////************************************************************************************
////************************************************************************************
void SmallDisplacementHpromElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mModesWeights = rCurrentProcessInfo[REDUCED_MODES_WEIGHTS];
    ClearNodalForces();
}

////************************************************************************************
////************************************************************************************

void SmallDisplacementHpromElement::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
}

////************************************************************************************
////************************************************************************************

void SmallDisplacementHpromElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    // create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

    // create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // set constitutive law flags:
    Flags& ConstitutiveLawOptions = Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    for (unsigned int PointNumber = 0;
         PointNumber < mConstitutiveLawVector.size(); PointNumber++)
    {
        // compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, PointNumber);

        // set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, PointNumber);

        // call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy(Values);

        // call the constitutive law to finalize the solution step
        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep(
            GetProperties(), GetGeometry(), Variables.N, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::InitializeMaterial()
{
    KRATOS_TRY

    // NOTE:
    // This is the standard (previous) implementation:
    // If we are here, it means that no one already set up the constitutive law
    // vector
    // through the method SetValue<CONSTITUTIVE_LAW>

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // Constitutive Law initialization
    if (mConstitutiveLawVector.size() != integration_points.size())
    {
        mConstitutiveLawVector.resize(integration_points.size());
    }
    else
    {
        // check whether the constitutive law pointers have been already set up
        bool already_set_up = true;
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            if (mConstitutiveLawVector[i] == NULL)
                already_set_up = false;
        }
        if (already_set_up)
        {
            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
            {
                mConstitutiveLawVector[i]->InitializeMaterial(
                    GetProperties(), GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
            }
            return; // if so, we are done here!
        }
    }

    // NOTE:
    // This is the standard (previous) implementation:
    // If we are here, it means that no one already set up the constitutive law
    // vector
    // through the method SetValue<CONSTITUTIVE_LAW>

    if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(), GetGeometry(),
                row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
        }
    }
    else
    {
        KRATOS_THROW_ERROR(
            std::logic_error,
            "a constitutive law needs to be specified for the element with ID ",
            this->Id())
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
    {
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->ResetMaterial(
                GetProperties(), GetGeometry(),
                row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
                                                                    GeneralVariables& rVariables,
                                                                    Vector& rVolumeForce,
                                                                    double& rIntegrationWeight)
{
    KRATOS_TRY
    const size_t number_of_nodes = GetGeometry().PointsNumber();
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        size_t index = dimension * i;
        for (unsigned int j = 0; j < dimension; j++)
        {
            rRightHandSideVector[index + j] +=
                rIntegrationWeight * rVariables.N[i] * rVolumeForce[j];
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateAndAddInternalForces(
    VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType InternalForces =
        rIntegrationWeight * prod(trans(rVariables.B), rVariables.StressVector);
    noalias(rRightHandSideVector) -= InternalForces;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                                          GeneralVariables& rVariables,
                                                          double& rIntegrationWeight)
{
    KRATOS_TRY

    // contributions to stiffness matrix calculated on the reference config
    noalias(rLeftHandSideMatrix) += prod(
        trans(rVariables.B), rIntegrationWeight * Matrix(prod(rVariables.ConstitutiveMatrix,
                                                              rVariables.B))); // to be optimized to remove the temporary

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::ClearNodalForces()
{
    KRATOS_TRY

    const size_t number_of_nodes = GetGeometry().PointsNumber();
    for (size_t i = 0; i < number_of_nodes; i++)
    {
        if (GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) &&
            GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE))
        {
            array_1d<double, 3>& ExternalForce =
                GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3>& InternalForce =
                GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

            GetGeometry()[i].SetLock();
            ExternalForce.clear();
            InternalForce.clear();
            GetGeometry()[i].UnSetLock();
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SmallDisplacementHpromElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t number_of_nodes = GetGeometry().PointsNumber();
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    if (rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            size_t index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3>& ExternalForce =
                GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            for (unsigned int j = 0; j < dimension; j++)
            {
                ExternalForce[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }

    if (rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            size_t index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3>& InternalForce =
                GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
            for (unsigned int j = 0; j < dimension; j++)
            {
                InternalForce[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            size_t index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3>& ForceResidual =
                GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (unsigned int j = 0; j < dimension; j++)
            {
                ForceResidual[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }

    KRATOS_CATCH("")
}

//************* COMPUTING  METHODS

void SmallDisplacementHpromElement::CalculateKinematics(GeneralVariables& rVariables,
                                                          const size_t& rPointNumber)
{
    KRATOS_TRY

    // Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De =
        rVariables.GetShapeFunctionsGradients();
    // Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    // Calculating the inverse of the jacobian and the parameters needed
    // [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix(rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    // Compute cartesian derivatives  [dN/dx_n]
    noalias(rVariables.DN_DX) = prod(DN_De[rPointNumber], InvJ);

    // Displacement Gradient H  [dU/dx_n]
    // this->CalculateDisplacementGradient( rVariables.H, rVariables.DN_DX );

    // Set Shape Functions Values for this integration point
    rVariables.N = row(Ncontainer, rPointNumber);

    // Compute the deformation matrix B
    // this->CalculateDeformationMatrixBbar(rVariables.B, rVariables.Bh,
    // rVariables.DN_DX );
    rVariables.B = mBMatrixVector[rPointNumber];

    // Compute infinitessimal strain
    this->CalculateInfinitesimalStrain(rVariables.B, rVariables.StrainVector);

    KRATOS_CATCH("")
}

//*************************COMPUTE DELTA
// POSITION*************************************
//************************************************************************************

Matrix& SmallDisplacementHpromElement::CalculateDeltaPosition(Matrix& rDeltaPosition)
{
    KRATOS_TRY

    /*const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  =
    GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement =
    GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) =
    CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    return rDeltaPosition;*/

    GeometryType& geom = GetGeometry();
    const size_t number_of_nodes = geom.PointsNumber();
    const size_t dimension = geom.WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>(number_of_nodes, dimension);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const NodeType& iNode = geom[i];
        rDeltaPosition(i, 0) = iNode.X() - iNode.X0();
        rDeltaPosition(i, 1) = iNode.Y() - iNode.Y0();
        if (dimension == 3)
            rDeltaPosition(i, 2) = iNode.Z() - iNode.Z0();
    }

    return rDeltaPosition;

    KRATOS_CATCH("")
}

void SmallDisplacementHpromElement::CalculateInfinitesimalStrain(const Matrix& rB,
                                                                   Vector& rStrainVector)
{
    //KRATOS_TRY

    const std::size_t dimension =
        static_cast<std::size_t>(GetGeometry().WorkingSpaceDimension());

    rStrainVector.clear();
    noalias(rStrainVector) = ZeroVector(mVoigtSize);

    if (dimension == 3) // 3d, voigt size = 6
    {
        for (unsigned int i = 0; i < mNumberOfModes; i++)
        {
            rStrainVector[0] += mModesWeights[i] * rB(0, i); // xx
            rStrainVector[1] += mModesWeights[i] * rB(1, i); // yy
            rStrainVector[2] += mModesWeights[i] * rB(2, i); // zz
            rStrainVector[3] += mModesWeights[i] * rB(3, i); // xy
            rStrainVector[4] += mModesWeights[i] * rB(4, i); // yz
            rStrainVector[5] += mModesWeights[i] * rB(5, i); // xz
        }
    }
    else if (mVoigtSize == 3) // 2d, voigt size = 3
    {
        for (unsigned int i = 0; i < mNumberOfModes; i++)
        {
            rStrainVector[0] += mModesWeights[i] * rB(0, i); // xx
            rStrainVector[1] += mModesWeights[i] * rB(1, i); // yy
            rStrainVector[2] += mModesWeights[i] * rB(2, i); // xy
        }
    }
    else if (mVoigtSize == 4) // 2d, voigt size = 4
    {
        for (unsigned int i = 0; i < mNumberOfModes; i++)
        {
            rStrainVector[0] += mModesWeights[i] * rB(0, i); // xx
            rStrainVector[1] += mModesWeights[i] * rB(1, i); // yy
            rStrainVector[2] += mModesWeights[i] * rB(2, i); // zz
            rStrainVector[3] += mModesWeights[i] * rB(3, i); // xy
        }
    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "something is wrong with the dimension", "")
    }

    //KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void SmallDisplacementHpromElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
{
    KRATOS_TRY
    const size_t number_of_nodes = GetGeometry().PointsNumber();
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if (dimension == 3)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 3 * i;

            rB(0, index + 0) = rDN_DX(i, 0);
            rB(1, index + 1) = rDN_DX(i, 1);
            rB(2, index + 2) = rDN_DX(i, 2);
            rB(3, index + 0) = rDN_DX(i, 1);
            rB(3, index + 1) = rDN_DX(i, 0);
            rB(4, index + 1) = rDN_DX(i, 2);
            rB(4, index + 2) = rDN_DX(i, 1);
            rB(5, index + 0) = rDN_DX(i, 2);
            rB(5, index + 2) = rDN_DX(i, 0);
        }
    }
    else if (dimension == 2 && mVoigtSize == 3)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 2 * i;

            rB(0, index + 0) = rDN_DX(i, 0);
            rB(0, index + 1) = 0.0;
            rB(1, index + 0) = 0.0;
            rB(1, index + 1) = rDN_DX(i, 1);
            rB(2, index + 0) = rDN_DX(i, 1);
            rB(2, index + 1) = rDN_DX(i, 0);
        }
    }
    else if (dimension == 2 && mVoigtSize == 4)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 2 * i;

            rB(0, index + 0) = rDN_DX(i, 0);
            rB(0, index + 1) = 0.0;
            rB(1, index + 0) = 0.0;
            rB(1, index + 1) = rDN_DX(i, 1);
            rB(2, index + 0) = 0.0;
            rB(2, index + 1) = 0.0;
            rB(3, index + 0) = rDN_DX(i, 1);
            rB(3, index + 1) = rDN_DX(i, 0);
        }
    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "something is wrong with the dimension", "")
    }
    KRATOS_CATCH("")
}

//************************************CALCULATE TOTAL
// MASS****************************
//************************************************************************************

double& SmallDisplacementHpromElement::CalculateTotalMass(double& rTotalMass,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    // rTotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY]; //not
    // accurate

    // Compute the Volume Change acumulated:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // reading integration points
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // compute element kinematics
        this->CalculateKinematics(Variables, PointNumber);

        // getting informations for integration
        double IntegrationWeight =
            Variables.detJ * integration_points[PointNumber].Weight();

        // compute point volume changes
        rTotalMass += GetProperties()[DENSITY] * IntegrationWeight;
    }

    if (dimension == 2)
        rTotalMass *= GetProperties()[THICKNESS];

    return rTotalMass;

    KRATOS_CATCH("")
}

//************************************CALCULATE VOLUME
// ACCELERATION*******************
//************************************************************************************

Vector& SmallDisplacementHpromElement::CalculateVolumeForce(Vector& rVolumeForce,
                                                              GeneralVariables& rVariables)
{
    KRATOS_TRY

    const size_t number_of_nodes = GetGeometry().PointsNumber();
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    if (rVolumeForce.size() != dimension)
        rVolumeForce.resize(dimension, false);

    noalias(rVolumeForce) = ZeroVector(dimension);

    for (unsigned int j = 0; j < number_of_nodes; j++)
    {
        if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION))
        { // it must be checked once at the begining
            // only
            array_1d<double, 3>& VolumeAcceleration =
                GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
            for (unsigned int i = 0; i < dimension; i++)
                rVolumeForce[i] += rVariables.N[j] * VolumeAcceleration[i];
        }
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t& integration_points_number =
        GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

    if (rOutput.size() != integration_points_number)
        rOutput.resize(integration_points_number, false);

    if (rVariable == VON_MISES_STRESS)
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // set constitutive law flags:
        Flags& ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            // set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            // call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            ComparisonUtilities EquivalentStress;
            rOutput[PointNumber] =
                EquivalentStress.CalculateVonMises(Variables.StressVector);
        }
    }
    else if (rVariable == STRAIN_ENERGY)
    {
        double Thickness = 1.0;
        const size_t dimension = GetGeometry().WorkingSpaceDimension();

        if (dimension == 2)
        {
            Thickness = GetProperties()[THICKNESS];
        }

        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // set constitutive law flags:
        Flags& ConstitutiveLawOptions = Values.GetOptions();

        // ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN); it would
        // return 0.0 strain since in small def. F = Identity
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);

        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            // to take in account previous step writing
            // if( mFinalizedStep ){
            // this->GetHistoricalVariables(Variables,PointNumber);
            //}
            // set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            double StrainEnergy = 0.0;

            // compute stresses and constitutive parameters
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
            mConstitutiveLawVector[PointNumber]->GetValue(STRAIN_ENERGY, StrainEnergy);

            rOutput[PointNumber] =
                Variables.detJ * integration_points[PointNumber].Weight() *
                Thickness * StrainEnergy; // 1/2 * sigma * epsilon
            // rOutput[PointNumber] = Variables.detJ *
            // integration_points[PointNumber].Weight() * Thickness *
            // StrainEnergy;
            // // 1/2 * sigma * epsilon

        } // for each gauss_point
    }

    else if (rVariable == GAUSS_WEIGHTS)
    {
        Vector& mAssignedIntegrationWeights = this->GetValue(INTEGRATION_POINT_WEIGHT);
        for (size_t ii = 0; ii < integration_points_number; ii++)
        {
            rOutput[ii] = mAssignedIntegrationWeights[ii];
        }
    }
    else
    {
        for (unsigned int ii = 0; ii < integration_points_number; ii++)
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t& integration_points_number =
        GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

    if (rOutput.size() != integration_points_number)
        rOutput.resize(integration_points_number);

    for (size_t i = 0; i < integration_points_number; i++)
    {
        rOutput[i] = ZeroVector(this->mVoigtSize);
    }

    if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR)
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // set constitutive law flags:
        Flags& ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        Vector& mAssignedIntegrationWeights = this->GetValue(INTEGRATION_POINT_WEIGHT);
        for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
        {
            if (mAssignedIntegrationWeights[PointNumber] < 0)
            {
                continue;
            }

            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            // set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            // call the constitutive law to update material variables
            if (rVariable == CAUCHY_STRESS_VECTOR)
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
            else
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

            if (rOutput[PointNumber].size() != Variables.StressVector.size())
                rOutput[PointNumber].resize(Variables.StressVector.size(), false);

            rOutput[PointNumber] = Variables.StressVector;
        }
    }
    else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR || rVariable == ALMANSI_STRAIN_VECTOR)
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // reading integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            if (rOutput[PointNumber].size() != Variables.StrainVector.size())
                rOutput[PointNumber].resize(Variables.StrainVector.size(), false);

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else
    {
        for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementHpromElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t& integration_points_number =
        GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
    const size_t dimension = GetGeometry().WorkingSpaceDimension();

    if (rOutput.size() != integration_points_number)
        rOutput.resize(integration_points_number);

    if (rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR)
    {
        std::vector<Vector> StressVector;

        if (rVariable == CAUCHY_STRESS_TENSOR)
            this->CalculateOnIntegrationPoints(
                CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);
        else
            this->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, StressVector,
                                               rCurrentProcessInfo);

        // loop integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            if (rOutput[PointNumber].size2() != dimension)
                rOutput[PointNumber].resize(dimension, dimension, false);

            rOutput[PointNumber] =
                MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    }
    else if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR || rVariable == ALMANSI_STRAIN_TENSOR)
    {
        std::vector<Vector> StrainVector;
        if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
            CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,
                                         StrainVector, rCurrentProcessInfo);
        else
            CalculateOnIntegrationPoints(ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo);

        // loop integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            if (rOutput[PointNumber].size2() != dimension)
                rOutput[PointNumber].resize(dimension, dimension, false);

            rOutput[PointNumber] =
                MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else if (rVariable == CONSTITUTIVE_MATRIX)
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // set constitutive law flags:
        Flags& ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        // reading integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            // set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            // call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
            // mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

            if (rOutput[PointNumber].size2() != Variables.ConstitutiveMatrix.size2())
                rOutput[PointNumber].resize(Variables.ConstitutiveMatrix.size1(),
                                            Variables.ConstitutiveMatrix.size2(), false);

            rOutput[PointNumber] = Variables.ConstitutiveMatrix;
        }
    }
    else if (rVariable ==
             DEFORMATION_GRADIENT) // VARIABLE SET FOR TRANSFER PURPOUSES
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // reading integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            if (rOutput[PointNumber].size2() != Variables.F.size2())
                rOutput[PointNumber].resize(Variables.F.size1(), Variables.F.size2(), false);

            rOutput[PointNumber] = Variables.F;
        }
    }
    else if (rVariable == REDUCED_MODES_MATRIX)
    {
        // create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        // reading integration points
        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);

            if (rOutput[PointNumber].size2() != Variables.B.size2())
                rOutput[PointNumber].resize(Variables.B.size1(), Variables.B.size2(), false);

            rOutput[PointNumber] = Variables.B;
        }
    }
    else
    {
        for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
        }
    }

    KRATOS_CATCH("")
}

// DECIMAL CORRECTION OF STRAINS
void SmallDisplacementHpromElement::DecimalCorrection(Vector& rVector)
{
    KRATOS_TRY
    for (unsigned int i = 0; i < rVector.size(); i++)
    {
        if (rVector[i] * rVector[i] < 1e-24)
        {
            rVector[i] = 0;
        }
    }
    KRATOS_CATCH("")
}

/**
 * This function provides the place to perform checks on the completeness of the
 * input. It is designed to be called only once (or anyway, not often) typically
 * at the
 * beginning of the calculations, so to verify that nothing is missing from the
 * input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int SmallDisplacementHpromElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const size_t dimension = this->GetGeometry().WorkingSpaceDimension();

    if (VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "VELOCITY has Key zero! (check "
                           "if the application is correctly registered",
                           "")
    if (DISPLACEMENT.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "DISPLACEMENT has Key zero! "
                           "(check if the application is correctly registered",
                           "")
    if (ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "ACCELERATION has Key zero! "
                           "(check if the application is correctly registered",
                           "")
    if (DENSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "DENSITY has Key zero! (check if "
                           "the application is correctly "
                           "registered",
                           "")
    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        if (!this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing variable DISPLACEMENT on node ",
                               this->GetGeometry()[i].Id())
        if (!this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) ||
            !this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) ||
            !this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z))
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "missing one of the dofs for the variable "
                               "DISPLACEMENT on node ",
                               GetGeometry()[i].Id())
    }
    if (!this->GetProperties().Has(CONSTITUTIVE_LAW))
    {
        KRATOS_THROW_ERROR(std::logic_error,
                           "constitutive law not provided for property ",
                           this->GetProperties().Id())
    }
    // verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(LawFeatures);
    bool correct_strain_measure = false;
    for (unsigned int i = 0; i < LawFeatures.mStrainMeasures.size(); i++)
    {
        if (LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
            correct_strain_measure = true;
    }
    if (!correct_strain_measure)
        KRATOS_THROW_ERROR(
            std::logic_error,
            "constitutive law is not compatible with the element type ",
            " Small Displacements ");

    // Verify that the body force is defined
    // if ( this->GetProperties().Has( BODY_FORCE ) == false )
    // {
    //     KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for
    //     property ", this->GetProperties().Id() )
    // }

    // verify that the constitutive law has the correct dimension
    if (dimension == 2)
    {
        // if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW
        // )->GetStrainSize()
        // != 3 )
        //     KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law
        //     used.
        //     This is a 2D element! expected strain size is 3 (el id = ) ",
        //     this->Id() ) //fails in some 2D cases, i.e. axisymmetric

        // if ( THICKNESS.Key() == 0 )
        //     KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key
        //     zero!
        //     (check if the application is correctly registered", "" ) //if is
        //     not
        //     read from model part it will not exist

        if (!this->GetProperties().Has(THICKNESS))
        {
            if (LawFeatures.mOptions.Is(ConstitutiveLaw::PLANE_STRAIN_LAW) ||
                LawFeatures.mOptions.Is(ConstitutiveLaw::AXISYMMETRIC_LAW))
            {
                this->GetProperties().SetValue(THICKNESS, 1.0);
            }
            else
            {
                KRATOS_THROW_ERROR(std::logic_error,
                                   "THICKNESS not provided for element ", this->Id())
            }
        }
    }
    else
    {
        if (this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() != 6)
            KRATOS_THROW_ERROR(std::logic_error,
                               "wrong constitutive law used. This "
                               "is a 3D element! expected strain "
                               "size is 6 (el id = ) ",
                               this->Id())
    }

    // check constitutive law
    /*for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(),
    rCurrentProcessInfo );
    }*/
    // FIXED: At this point the constitutive law vector is not set yet.
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->Check(this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

    // check if it is in the XY plane for 2D case

    return 0;

    KRATOS_CATCH("");
}

void SmallDisplacementHpromElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod", IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

void SmallDisplacementHpromElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    int IntMethod;
    rSerializer.load("IntegrationMethod", IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
