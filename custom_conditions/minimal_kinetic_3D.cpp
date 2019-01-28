#include "minimal_kinetic_3D.hpp"
#include "multiscale_rom_application.h"
//#include "multiscale_rom_application_variables.h"
//
//#include "includes/constitutive_law.h"
//#include "includes/define.h"
//#include "includes/element.h"
//#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticCondition3D::MinimalKineticCondition3D() : Condition()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticCondition3D::MinimalKineticCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
MinimalKineticCondition3D::MinimalKineticCondition3D(IndexType NewId,
                                                     GeometryType::Pointer pGeometry,
                                                     PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

MinimalKineticCondition3D::MinimalKineticCondition3D(MinimalKineticCondition3D const& rOther)
    : Condition(rOther), mThisIntegrationMethod(rOther.mThisIntegrationMethod)
{
}

//****************** DESTRUCTOR ******************************************************
//************************************************************************************
MinimalKineticCondition3D::~MinimalKineticCondition3D()
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************
MinimalKineticCondition3D& MinimalKineticCondition3D::operator=(MinimalKineticCondition3D const& rOther)
{
    Condition::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
Condition::Pointer MinimalKineticCondition3D::Create(IndexType NewId,
                                                     NodesArrayType const& ThisNodes,
                                                     PropertiesType::Pointer pProperties) const
{
    return MinimalKineticCondition3D::Pointer(new MinimalKineticCondition3D(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
}

int MinimalKineticCondition3D::Check(const ProcessInfo& rCurrentProcessInfo)
{
    // return Condition::Check(rCurrentProcessInfo);
    return 0;
}

void MinimalKineticCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int StrainComp = 6;
    unsigned int n_dofs = number_of_nodes * dimension +
                          StrainComp; // 6 = Num components of (u_fl tens n)^s

    Vector currentValues(n_dofs, 0.0);
    Matrix rNintMatrix(dimension, dimension * number_of_nodes, 0.0);

    GeometryType& geom = GetGeometry();

    // resize system matrix and vector
    if (rLeftHandSideMatrix.size1() != n_dofs || rLeftHandSideMatrix.size2() != n_dofs)
        rLeftHandSideMatrix.resize(n_dofs, n_dofs, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_dofs, n_dofs);

    if (rRightHandSideVector.size() != n_dofs)
        rRightHandSideVector.resize(n_dofs, false);

    noalias(rRightHandSideVector) = ZeroVector(n_dofs);

    this->CalculateIntegralOfShapeFunctions(rNintMatrix, rCurrentProcessInfo);

    // Auxiliar magic node
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];


    // Displacements
    currentValues(0) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(1) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(2) = geom[0].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(3) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(4) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(5) = geom[1].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(6) = geom[2].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(7) = geom[2].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(8) = geom[2].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(9)  = geom[3].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(10) = geom[3].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(11) = geom[3].FastGetSolutionStepValue(DISPLACEMENT_Z);

    // Lagrange Multipliers
    currentValues(12) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
    currentValues(13) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
    currentValues(14) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);
    currentValues(15) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_4);
    currentValues(16) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_5);
    currentValues(17) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_6);

    // Vector V1(dimension, 0.0);
    // Vector V2(dimension, 0.0);
    // Vector V3(dimension, 0.0);

    array_1d<double, 3> V1(dimension);
    array_1d<double, 3> V2(dimension);
    array_1d<double, 3> V3(dimension);

    // Unitary normal vector V1
    V1(0) = geom[1].X0() - geom[0].X0();
    V1(1) = geom[1].Y0() - geom[0].Y0();
    V1(2) = geom[1].Z0() - geom[0].Z0();
    V1 /= MathUtils<double>::Norm3(V1);

    // Unitary normal vector V2
    V2(0) = geom[3].X0() - geom[0].X0();
    V2(1) = geom[3].Y0() - geom[0].Y0();
    V2(2) = geom[3].Z0() - geom[0].Z0();
    V2 /= MathUtils<double>::Norm3(V2);

    // Normal vector to the element surface
    MathUtils<double>::CrossProduct(V3, V1, V2);

    // compute the outward normal vector
    Matrix TensNormal(StrainComp, dimension, 0.0);
    // TensNormal = ZeroMatrix( 6, dimension );

    TensNormal(0, 0) = V3(0);  // xx
    TensNormal(1, 1) = V3(1);  // yy
    TensNormal(2, 2) = V3(2);  // zz
    TensNormal(3, 0) = V3(1);  // xy
    TensNormal(3, 1) = V3(0);
    TensNormal(4, 1) = V3(2);  // yz
    TensNormal(4, 2) = V3(1);
    TensNormal(5, 0) = V3(2);  // xz
    TensNormal(5, 2) = V3(0);

    // constraint matrix
    Matrix ElemConstraintMatrix = prod(TensNormal, rNintMatrix);

    unsigned int indexi = number_of_nodes * dimension;
    Matrix& K = rLeftHandSideMatrix;

    // Local assembly integral of shape functions matrix (for computing
    // constraint matrix - multiscale cases)
    for (unsigned int i = 0; i < StrainComp; i++)
    {
        for (unsigned int j = 0; j < number_of_nodes * dimension; j++)
        {
            K(indexi + i, j) = ElemConstraintMatrix(i, j);
            K(j, indexi + i) = ElemConstraintMatrix(i, j);
            // K(indexi, j) = ElemConstraintMatrix(i,j);
            // K(j, indexi) = ElemConstraintMatrix(i,j);
        }
    }

    // residual force
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, currentValues);

    ////KRATOS_WATCH(currentValues)
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH(rRightHandSideVector)
}

void MinimalKineticCondition3D::CalculateIntegralOfShapeFunctions(MatrixType& rNintMatrix,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

    // reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //this->GetGeometry().DeterminantOfJacobian(Variables.detJ,mThisIntegrationMethod);

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int MatSize = dimension * number_of_nodes;
    // MatrixType rNintMatrix( dimension, MatSize );

    if (rNintMatrix.size1() != dimension)
        rNintMatrix.resize(dimension, MatSize, false);

    noalias(rNintMatrix) = ZeroMatrix(dimension, MatSize);

    // Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = Variables.GetShapeFunctions();

    Matrix& KCond = rNintMatrix;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // Set Shape Functions Values for this integration point
        Variables.N = row(Ncontainer, PointNumber);

        // Calculating the inverse of the jacobian and the parameters needed
        // [d£/dx_n]
        Matrix InvJ;
        //MathUtils<double>::InvertMatrix(Variables.J[PointNumber], InvJ,
        //                                Variables.detJ);

        // DetJ using MathUtils (last version)
        Variables.detJ = MathUtils<double>::GeneralizedDet(Variables.J[PointNumber]);

        // Vicente's implementation (in the elemental routine)
        //Variables.detJ=GetGeometry().DeterminantOfJacobian(PointNumber,mThisIntegrationMethod);

        // calculating weights for integration on the "reference configuration"
        double IntegrationWeight =
            integration_points[PointNumber].Weight() * Variables.detJ;
        IntegrationWeight = this->CalculateIntegrationWeight(IntegrationWeight);

        // Local assembly integral of shape functions matrix (for computing
        // constraint matrix - multiscale cases)
        unsigned int indexi = 0;

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            KCond(0, indexi) += Variables.N[i] * IntegrationWeight;
            KCond(1, indexi + 1) += Variables.N[i] * IntegrationWeight;
            KCond(2, indexi + 2) += Variables.N[i] * IntegrationWeight;

            indexi += dimension;
        }
        // KRATOS_WATCH( rNintMatrix )
    }

    /*
    //noalias(rRightHandSideVector) = prod( MassMatrix,
       CurrentAccelerationVector );
    //KRATOS_WATCH( rRightHandSideVector )
    */
    KRATOS_CATCH("")
}

void MinimalKineticCondition3D::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0, 0, false);
}

void MinimalKineticCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
    CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}

void MinimalKineticCondition3D::EquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    // unsigned int n_lag_dofs = 3;
    // unsigned int n_dofs = 4 + n_lag_dofs;

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int n_dofs = dimension * number_of_nodes +
                          6; // 6 = Num components of (u_fl tens n)^s

    if (rResult.size() != n_dofs)
        rResult.resize(n_dofs, false);
    GeometryType& geom = GetGeometry();

    rResult[0] = geom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = geom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = geom[0].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[3] = geom[1].GetDof(DISPLACEMENT_X).EquationId();
    rResult[4] = geom[1].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[5] = geom[1].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[6] = geom[2].GetDof(DISPLACEMENT_X).EquationId();
    rResult[7] = geom[2].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[8] = geom[2].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[9] = geom[3].GetDof(DISPLACEMENT_X).EquationId();
    rResult[10] = geom[3].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[11] = geom[3].GetDof(DISPLACEMENT_Z).EquationId();

    // Auxiliar magic node
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];

    // KRATOS_WATCH(*pNode)

    // Lagrange Multipliers
    rResult[12] = pNode->GetDof(LAGRANGE_MULTIPLIER_1).EquationId();
    rResult[13] = pNode->GetDof(LAGRANGE_MULTIPLIER_2).EquationId();
    rResult[14] = pNode->GetDof(LAGRANGE_MULTIPLIER_3).EquationId();
    rResult[15] = pNode->GetDof(LAGRANGE_MULTIPLIER_4).EquationId();
    rResult[16] = pNode->GetDof(LAGRANGE_MULTIPLIER_5).EquationId();
    rResult[17] = pNode->GetDof(LAGRANGE_MULTIPLIER_6).EquationId();
}

void MinimalKineticCondition3D::GetDofList(DofsVectorType& rConditionDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& geom = GetGeometry();
    const unsigned int dimension = geom.WorkingSpaceDimension();
    unsigned int nr_of_nodes = 4;
    unsigned int nr_lagrang_dofs = 6;
    unsigned int nr_dofs = dimension * nr_of_nodes + nr_lagrang_dofs;

    // rConditionDofList.resize(0);
    if (rConditionDofList.size() != nr_dofs)
        rConditionDofList.resize(nr_dofs, false);

    // first node
    // rConditionDofList.push_back(geom[0].pGetDof(DISPLACEMENT_X));
    rConditionDofList[0] = geom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[1] = geom[0].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[2] = geom[0].pGetDof(DISPLACEMENT_Z);
    // second node
    rConditionDofList[3] = geom[1].pGetDof(DISPLACEMENT_X);
    rConditionDofList[4] = geom[1].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[5] = geom[1].pGetDof(DISPLACEMENT_Z);
    // third node
    rConditionDofList[6] = geom[2].pGetDof(DISPLACEMENT_X);
    rConditionDofList[7] = geom[2].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[8] = geom[2].pGetDof(DISPLACEMENT_Z);
    // fourth node
    rConditionDofList[9] = geom[3].pGetDof(DISPLACEMENT_X);
    rConditionDofList[10] = geom[3].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[11] = geom[3].pGetDof(DISPLACEMENT_Z);

    // Auxiliar magic node
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];

    // Lagrange Multipliers
    rConditionDofList[12] = pNode->pGetDof(LAGRANGE_MULTIPLIER_1);
    rConditionDofList[13] = pNode->pGetDof(LAGRANGE_MULTIPLIER_2);
    rConditionDofList[14] = pNode->pGetDof(LAGRANGE_MULTIPLIER_3);
    rConditionDofList[15] = pNode->pGetDof(LAGRANGE_MULTIPLIER_4);
    rConditionDofList[16] = pNode->pGetDof(LAGRANGE_MULTIPLIER_5);
    rConditionDofList[17] = pNode->pGetDof(LAGRANGE_MULTIPLIER_6);
}

//************************************************************************************
//************************************************************************************

void MinimalKineticCondition3D::InitializeGeneralVariables(GeneralVariables& rVariables,
                                                           const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

    rVariables.Initialize(voigt_size, dimension, number_of_nodes);

    // needed parameters for consistency with the general constitutive law:
    // small displacements
    // rVariables.detF  = 1.0;
    // rVariables.F     = identity_matrix<double>(dimension);

    // set variables including all integration points values

    // reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod));

    // reading shape functions local gradients
    // rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients(
    // mThisIntegrationMethod ));

    // calculating the current jacobian from cartesian coordinates to parent
    // coordinates for all integration points [dx_n+1/d£]
    //rVariables.j = GetGeometry().Jacobian(rVariables.j, mThisIntegrationMethod);

    // in this case, is not necessary to compute F, because the use of small
    // strain setting
    // Calculate Delta Position
    // rVariables.DeltaPosition =
    // CalculateDeltaPosition(rVariables.DeltaPosition);

    // Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    // calculating the reference jacobian from cartesian coordinates to parent
    // coordinates for all integration points [dx_n/d£]
    // compunting Jacobian using updated coordinates minus the increment of
    // displacement in order to get the reference coordinates, even if the
    //problem is under small strains setting
    rVariables.J = GetGeometry().Jacobian(rVariables.J, mThisIntegrationMethod,
                                          rVariables.DeltaPosition);

    //rVariables.J = GetGeometry().Jacobian(rVariables.J, mThisIntegrationMethod);

    //KRATOS_WATCH(rVariables.J)
}


//*************************COMPUTE DELTA POSITION******************************
//*****************************************************************************

Matrix& MinimalKineticCondition3D::CalculateDeltaPosition(Matrix& rDeltaPosition)
{
  KRATOS_TRY

    GeometryType& geom = GetGeometry();
    const unsigned int number_of_nodes = geom.PointsNumber();
    unsigned int dimension = geom.WorkingSpaceDimension();

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

//************************************************************************************
//************************************************************************************

double& MinimalKineticCondition3D::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension == 2)
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}

void MinimalKineticCondition3D::GetValuesVector(Vector& Values, int Step)
{
}

void MinimalKineticCondition3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void MinimalKineticCondition3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}
}
