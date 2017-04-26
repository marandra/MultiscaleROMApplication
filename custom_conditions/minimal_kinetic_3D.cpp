#include "minimal_kinetic_3D.hpp"
#include "multiscale_rom_application.h"
#include "multiscale_rom_application_variables.h"

#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"

/// Small Displacement Element for 3D and 2D geometries.

/**
 * Implements a Small Displacement Lagrangian definition for structural
 * analysis.
 * This works for arbitrary geometries in 3D and 2D
 */

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticCondition3D::MinimalKineticCondition3D() : Condition()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
// MinimalKineticCondition3D::MinimalKineticCondition3D(IndexType NewId, const
// NodesArrayType& ThisNodes):
//        Condition(NewId,ThisNodes)
//{
//}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticCondition3D::MinimalKineticCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

//******************************COPY
//CONSTRUCTOR**************************************
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

//****************** DESTRUCTOR
//******************************************************
//************************************************************************************
MinimalKineticCondition3D::~MinimalKineticCondition3D()
{
}

//**************** GETTING METHODS
//***************************************************
//************************************************************************************
MinimalKineticCondition3D::IntegrationMethod MinimalKineticCondition3D::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//*******************************ASSIGMENT
//OPERATOR***********************************
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
    KRATOS_TRY;

    // return Condition::Check(rCurrentProcessInfo);
    return 0;

    KRATOS_CATCH("");
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
    // MatrixType& rNintMatrix(dimension,dimension*number_of_nodes);
    // Matrix rNintMatrix;

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

    // Lagrange Multipliers
    currentValues[12] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
    currentValues[13] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
    currentValues[14] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);
    currentValues[15] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_4);
    currentValues[16] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_5);
    currentValues[17] = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_6);

    // currentValues(12) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
    // currentValues(13) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
    // currentValues(14) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);
    // currentValues(15) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_4);
    // currentValues(16) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_5);
    // currentValues(17) =
    // geom[4].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_6);

    // Vector V1(dimension, 0.0);
    // Vector V2(dimension, 0.0);
    // Vector V3(dimension, 0.0);

    array_1d<double, 3> V1(dimension, 0.0);
    array_1d<double, 3> V2(dimension, 0.0);
    array_1d<double, 3> V3(dimension, 0.0);

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

    /*
    //Matrix& K = rLeftHandSideMatrix;
    TensNormal(0, 0) = V3(0);       // xx
    TensNormal(1, 1) = V3(1);       // yy
    TensNormal(2, 2) = V3(2);       // zz
    TensNormal(3, 1) = 0.5*V3(0);   // xy
    TensNormal(3, 2) = 0.5*V3(1);
    TensNormal(4, 1) = 0.5*V3(1);   // yz
    TensNormal(4, 2) = 0.5*V3(2);
    TensNormal(5, 0) = 0.5*V3(0);   // xz
    TensNormal(5, 2) = 0.5*V3(2);
    */

    TensNormal(0, 0) = V3(0);       // xx
    TensNormal(1, 1) = V3(1);       // yy
    TensNormal(2, 2) = V3(2);       // zz
    TensNormal(3, 0) = 0.5 * V3(1); // xy
    TensNormal(3, 1) = 0.5 * V3(0);
    TensNormal(4, 1) = 0.5 * V3(2); // yz
    TensNormal(4, 2) = 0.5 * V3(1);
    TensNormal(5, 0) = 0.5 * V3(2); // xz
    TensNormal(5, 2) = 0.5 * V3(0);

    KRATOS_WATCH(TensNormal)
    KRATOS_WATCH(rNintMatrix)
    // KRATOS_WATCH(V3)

    // constraint matrix
    Matrix ElemConstraintMatrix = prod(TensNormal, rNintMatrix);

    KRATOS_WATCH(ElemConstraintMatrix)

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

    // form residual
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, currentValues);
    ////KRATOS_WATCH(currentValues)
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH(rRightHandSideVector)
}

//***********************************************************************************
// CrossProduct calculates the cross product of two 3d vectors. a x b.
// @return cross = a x b
//***********************************************************************************
// void MinimalKineticCondition3D::CrossProduct(Vector& cross,
//                                             const Vector& a,
//                                             const Vector& b)
//{
//    //array_1d<double, 3> cross;
//    cross[0] = a[1] * b[2] - a[2] * b[1];
//    cross[1] = a[2] * b[0] - a[0] * b[2];
//    cross[2] = a[0] * b[1] - a[1] * b[0];
//    //return cross;
//}

//************************************************************************************
// Compute the integral of the shape functions within the finite element (for 3D
// cases)
//************************************************************************************
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
        MathUtils<double>::InvertMatrix(Variables.J[PointNumber], InvJ,
                                        Variables.detJ);

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
        rResult.resize(n_dofs);
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
        rConditionDofList.resize(nr_dofs);

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
    rVariables.j = GetGeometry().Jacobian(rVariables.j, mThisIntegrationMethod);

    // in this case, is not necessary to compute F, because the use of small
    // strain setting
    // Calculate Delta Position
    // rVariables.DeltaPosition =
    // CalculateDeltaPosition(rVariables.DeltaPosition);

    // calculating the reference jacobian from cartesian coordinates to parent
    // coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian(rVariables.J, mThisIntegrationMethod,
                                          rVariables.DeltaPosition);
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
    /*    PeriodicVariablesContainer const& rPeriodicVariables =
       this->GetProperties().GetValue(PERIODIC_VARIABLES);
        const unsigned int BlockSize = rPeriodicVariables.size();
        const unsigned int LocalSize = 2 * BlockSize; // Total contribution size
       = 2 nodes * num dofs

        if (Values.size() != LocalSize)
            Values.resize(LocalSize,false);

        unsigned int LocalIndex = 0;

        for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar =
       rPeriodicVariables.DoubleVariablesBegin();
                itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
        {
            Values[LocalIndex] =
       this->GetGeometry()[0].FastGetSolutionStepValue(*itDVar,Step);
            Values[LocalIndex+BlockSize] =
       this->GetGeometry()[1].FastGetSolutionStepValue(*itDVar,Step);
            ++LocalIndex;
        }

        for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar =
       rPeriodicVariables.VariableComponentsBegin();
                itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
        {
            Values[LocalIndex] =
       this->GetGeometry()[0].FastGetSolutionStepValue(*itCVar,Step);
            Values[LocalIndex+BlockSize] =
       this->GetGeometry()[1].FastGetSolutionStepValue(*itCVar,Step);
            ++LocalIndex;
        }*/
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