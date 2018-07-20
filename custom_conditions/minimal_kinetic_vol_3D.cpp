#include "minimal_kinetic_vol_3D.hpp"
#include "multiscale_rom_application.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticVolCondition3D::MinimalKineticVolCondition3D() : Condition()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MinimalKineticVolCondition3D::MinimalKineticVolCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
MinimalKineticVolCondition3D::MinimalKineticVolCondition3D(IndexType NewId,
                                                     GeometryType::Pointer pGeometry,
                                                     PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

MinimalKineticVolCondition3D::MinimalKineticVolCondition3D(MinimalKineticVolCondition3D const& rOther)
    : Condition(rOther), mThisIntegrationMethod(rOther.mThisIntegrationMethod)
{
}

//****************** DESTRUCTOR ******************************************************
//************************************************************************************
MinimalKineticVolCondition3D::~MinimalKineticVolCondition3D()
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************
MinimalKineticVolCondition3D& MinimalKineticVolCondition3D::operator=(MinimalKineticVolCondition3D const& rOther)
{
    Condition::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
Condition::Pointer MinimalKineticVolCondition3D::Create(IndexType NewId,
                                                     NodesArrayType const& ThisNodes,
                                                     PropertiesType::Pointer pProperties) const
{
    return MinimalKineticVolCondition3D::Pointer(new MinimalKineticVolCondition3D(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
}

int MinimalKineticVolCondition3D::Check(const ProcessInfo& rCurrentProcessInfo)
{
    // return Condition::Check(rCurrentProcessInfo);
    return 0;
}


//************************************************************************************
//************************************************************************************
void MinimalKineticVolCondition3D::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
{
    const std::size_t nr_nodes = GetGeometry().PointsNumber();
    const std::size_t nr_dimensions = GetGeometry().WorkingSpaceDimension();
    const std::size_t nr_components = 6;
    const std::size_t nr_dofs = nr_nodes * nr_dimensions + nr_components;
    const std::size_t index = nr_nodes * nr_dimensions;
    Vector currentValues(nr_dofs, 0.0);
    GeometryType& geom = GetGeometry();

    if (rLeftHandSideMatrix.size1() != nr_dofs || rLeftHandSideMatrix.size2() != nr_dofs)
        rLeftHandSideMatrix.resize(nr_dofs, nr_dofs, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(nr_dofs, nr_dofs);

    if (rRightHandSideVector.size() != nr_dofs)
        rRightHandSideVector.resize(nr_dofs, false);
    noalias(rRightHandSideVector) = ZeroVector(nr_dofs);

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

    currentValues(12) = geom[4].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(13) = geom[4].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(14) = geom[4].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(15) = geom[5].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(16) = geom[5].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(17) = geom[5].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(18) = geom[6].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(19) = geom[6].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(20) = geom[6].FastGetSolutionStepValue(DISPLACEMENT_Z);

    currentValues(21) = geom[7].FastGetSolutionStepValue(DISPLACEMENT_X);
    currentValues(22) = geom[7].FastGetSolutionStepValue(DISPLACEMENT_Y);
    currentValues(23) = geom[7].FastGetSolutionStepValue(DISPLACEMENT_Z);

    // Lagrange Multipliers
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];
    currentValues(24) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
    currentValues(25) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
    currentValues(26) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);
    currentValues(27) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_4);
    currentValues(28) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_5);
    currentValues(29) = pNode->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_6);

    const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    Matrix J0(nr_dimensions, nr_dimensions, false);
    Matrix InvJ0(nr_dimensions, nr_dimensions, false);
    Matrix DN_DX(nr_nodes, nr_dimensions, false);
    Matrix B_matrix(nr_components, index, false);

    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        double detJ0 = this->CalculateDerivativesOnReferenceConfiguration(
                J0, InvJ0, DN_DX, point_number, this->GetIntegrationMethod());
        this->CalculateB(B_matrix, DN_DX);
        double integration_weight = integration_points[point_number].Weight() * detJ0;

        // local assembly
        for (std::size_t i = 0; i < nr_components; i++)
        {
            for (std::size_t j = 0; j < nr_nodes * nr_dimensions; j++)
            {
                rLeftHandSideMatrix(index + i, j) += integration_weight * (-B_matrix(i, j));
                rLeftHandSideMatrix(j, index + i) += integration_weight * (-B_matrix(i, j));
            }
        }
    }

    // residual force
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, currentValues);


}
/***********************************************************************************/
/***********************************************************************************/

double MinimalKineticVolCondition3D::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
{
    GeometryType& r_geom = GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geom,
        r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix& rDN_De =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}


/***********************************************************************************/
/***********************************************************************************/
void MinimalKineticVolCondition3D::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    const SizeType nr_nodes = GetGeometry().PointsNumber();

    rB.clear();

    for ( SizeType i = 0; i < nr_nodes; ++i ) {
        rB( 0, i*3     ) = rDN_DX( i, 0 );
        rB( 1, i*3 + 1 ) = rDN_DX( i, 1 );
        rB( 2, i*3 + 2 ) = rDN_DX( i, 2 );
        rB( 3, i*3     ) = rDN_DX( i, 1 );
        rB( 3, i*3 + 1 ) = rDN_DX( i, 0 );
        rB( 4, i*3 + 1 ) = rDN_DX( i, 2 );
        rB( 4, i*3 + 2 ) = rDN_DX( i, 1 );
        rB( 5, i*3     ) = rDN_DX( i, 2 );
        rB( 5, i*3 + 2 ) = rDN_DX( i, 0 );
    }
}


/***********************************************************************************/
/***********************************************************************************/
void MinimalKineticVolCondition3D::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0, 0, false);
}


/***********************************************************************************/
/***********************************************************************************/
void MinimalKineticVolCondition3D::CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType dummy;
    CalculateLocalSystem(dummy, rRightHandSideVector, rCurrentProcessInfo);
}


/***********************************************************************************/
/***********************************************************************************/
void MinimalKineticVolCondition3D::EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo
    )
{
    const std::size_t nr_nodes = GetGeometry().PointsNumber();
    const std::size_t nr_dimensions = GetGeometry().WorkingSpaceDimension();
    const std::size_t nr_components = 6;
    const std::size_t nr_dofs = nr_dimensions * nr_nodes + nr_components;

    if (rResult.size() != nr_dofs)
        rResult.resize(nr_dofs);

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

    rResult[12] = geom[4].GetDof(DISPLACEMENT_X).EquationId();
    rResult[13] = geom[4].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[14] = geom[4].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[15] = geom[5].GetDof(DISPLACEMENT_X).EquationId();
    rResult[16] = geom[5].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[17] = geom[5].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[18] = geom[6].GetDof(DISPLACEMENT_X).EquationId();
    rResult[19] = geom[6].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[20] = geom[6].GetDof(DISPLACEMENT_Z).EquationId();

    rResult[21] = geom[7].GetDof(DISPLACEMENT_X).EquationId();
    rResult[22] = geom[7].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[23] = geom[7].GetDof(DISPLACEMENT_Z).EquationId();

    // Lagrange Multipliers
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];
    rResult[24] = pNode->GetDof(LAGRANGE_MULTIPLIER_1).EquationId();
    rResult[25] = pNode->GetDof(LAGRANGE_MULTIPLIER_2).EquationId();
    rResult[26] = pNode->GetDof(LAGRANGE_MULTIPLIER_3).EquationId();
    rResult[27] = pNode->GetDof(LAGRANGE_MULTIPLIER_4).EquationId();
    rResult[28] = pNode->GetDof(LAGRANGE_MULTIPLIER_5).EquationId();
    rResult[29] = pNode->GetDof(LAGRANGE_MULTIPLIER_6).EquationId();
}

void MinimalKineticVolCondition3D::GetDofList(DofsVectorType& rConditionDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& geom = GetGeometry();
    const unsigned int dimension = geom.WorkingSpaceDimension();
    //unsigned int nr_of_nodes = 4;
    unsigned int nr_of_nodes = GetGeometry().PointsNumber();
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
    // fourth node
    rConditionDofList[12] = geom[4].pGetDof(DISPLACEMENT_X);
    rConditionDofList[13] = geom[4].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[14] = geom[4].pGetDof(DISPLACEMENT_Z);
    // fourth node
    rConditionDofList[15] = geom[5].pGetDof(DISPLACEMENT_X);
    rConditionDofList[16] = geom[5].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[17] = geom[5].pGetDof(DISPLACEMENT_Z);
    // fourth node
    rConditionDofList[18] = geom[6].pGetDof(DISPLACEMENT_X);
    rConditionDofList[19] = geom[6].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[20] = geom[6].pGetDof(DISPLACEMENT_Z);
    // fourth node
    rConditionDofList[21] = geom[7].pGetDof(DISPLACEMENT_X);
    rConditionDofList[22] = geom[7].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[23] = geom[7].pGetDof(DISPLACEMENT_Z);

    // Auxiliar magic node
    Node<3>::Pointer pNode = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_NODE];

    // Lagrange Multipliers
    rConditionDofList[24] = pNode->pGetDof(LAGRANGE_MULTIPLIER_1);
    rConditionDofList[25] = pNode->pGetDof(LAGRANGE_MULTIPLIER_2);
    rConditionDofList[26] = pNode->pGetDof(LAGRANGE_MULTIPLIER_3);
    rConditionDofList[27] = pNode->pGetDof(LAGRANGE_MULTIPLIER_4);
    rConditionDofList[28] = pNode->pGetDof(LAGRANGE_MULTIPLIER_5);
    rConditionDofList[29] = pNode->pGetDof(LAGRANGE_MULTIPLIER_6);
}

//************************************************************************************
//************************************************************************************

void MinimalKineticVolCondition3D::InitializeGeneralVariables(GeneralVariables& rVariables,
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

Matrix& MinimalKineticVolCondition3D::CalculateDeltaPosition(Matrix& rDeltaPosition)
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


void MinimalKineticVolCondition3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void MinimalKineticVolCondition3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}
}
