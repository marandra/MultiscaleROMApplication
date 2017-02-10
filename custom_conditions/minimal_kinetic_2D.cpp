#include "custom_conditions/minimal_kinetic_2D.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{

MinimalKineticCondition2D::MinimalKineticCondition2D(IndexType NewId):
    Condition(NewId)
{
}

MinimalKineticCondition2D::MinimalKineticCondition2D(IndexType NewId, const NodesArrayType& ThisNodes):
    Condition(NewId,ThisNodes)
{
}

MinimalKineticCondition2D::MinimalKineticCondition2D(IndexType NewId, GeometryType::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

MinimalKineticCondition2D::MinimalKineticCondition2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

MinimalKineticCondition2D::MinimalKineticCondition2D(MinimalKineticCondition2D const& rOther):
    Condition(rOther)
{
}

MinimalKineticCondition2D::~MinimalKineticCondition2D()
{
}

MinimalKineticCondition2D& MinimalKineticCondition2D::operator =(MinimalKineticCondition2D const& rOther)
{
    Condition::operator =(rOther);

    return *this;
}

Condition::Pointer MinimalKineticCondition2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return MinimalKineticCondition2D::Pointer(new MinimalKineticCondition2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

int MinimalKineticCondition2D::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    return Condition::Check(rCurrentProcessInfo);

    KRATOS_CATCH("");
}

void MinimalKineticCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    unsigned int n_dofs = 7;
    Vector currentValues(n_dofs, 0.0);
    GeometryType& geom = GetGeometry();

    // resize system matrix and vector
    if (rLeftHandSideMatrix.size1() != n_dofs ||
        rLeftHandSideMatrix.size2() != n_dofs)
        rLeftHandSideMatrix.resize(n_dofs, n_dofs, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_dofs, n_dofs);

    if (rRightHandSideVector.size() != n_dofs)
        rRightHandSideVector.resize(n_dofs, false);
    noalias(rRightHandSideVector) = ZeroVector(n_dofs);

    currentValues(4) = geom[2].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
    currentValues(5) = geom[2].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
    currentValues(6) = geom[2].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);

    // compute the outward normal vector
    double x0 = geom[0].X0();
    double y0 = geom[0].Y0();
    double x1 = geom[1].X0();
    double y1 = geom[1].Y0();
    double tx = (x1 - x0);
    double ty = (y1 - y0);
    //double L  = std::sqrt(tx * tx + ty * ty);
    double nx = 0.5 * ty;
    double ny = 0.5 * (-tx);

    Matrix& K = rLeftHandSideMatrix;
    K(4, 0) = nx;
    K(4, 2) = nx;
    K(0, 4) = nx;
    K(2, 4) = nx;
    K(5, 1) = ny;
    K(5, 3) = ny;
    K(1, 5) = ny;
    K(3, 5) = ny;
    K(6, 0) = ny;
    K(6, 1) = nx;
    K(6, 2) = ny;
    K(6, 3) = nx;
    K(0, 6) = ny;
    K(1, 6) = nx;
    K(2, 6) = ny;
    K(3, 6) = nx;
    // form residual
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, currentValues);

//KRATOS_WATCH("DEBUG CALCULATE LOCAL SYSTEM")
//KRATOS_WATCH(this->Id())
//KRATOS_WATCH(geom[0].Id())
//KRATOS_WATCH(geom[1].Id())
//KRATOS_WATCH(geom[2].Id())
//KRATOS_WATCH(K)
}

void MinimalKineticCondition2D::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0,0,false);
}

void MinimalKineticCondition2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType dummy;
    CalculateLocalSystem(dummy,rRightHandSideVector,rCurrentProcessInfo);
}

void MinimalKineticCondition2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    unsigned int n_lag_dofs = 3;
    unsigned int n_dofs = 4 + n_lag_dofs;

    if (rResult.size() != n_dofs) rResult.resize(n_dofs);
    GeometryType& geom = GetGeometry();
    rResult[0] = geom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = geom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = geom[1].GetDof(DISPLACEMENT_X).EquationId();
    rResult[3] = geom[1].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[4] = geom[2].GetDof(LAGRANGE_MULTIPLIER_1).EquationId();
    rResult[5] = geom[2].GetDof(LAGRANGE_MULTIPLIER_2).EquationId();
    rResult[6] = geom[2].GetDof(LAGRANGE_MULTIPLIER_3).EquationId();
}

void MinimalKineticCondition2D::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& geom = GetGeometry();
    const unsigned int dimension = geom.WorkingSpaceDimension();
    unsigned int nr_of_nodes = 2;
    unsigned int nr_lagrang_dofs = 3;
    unsigned int nr_dofs = 2 * nr_of_nodes + nr_lagrang_dofs;

    //rConditionDofList.resize(0);
    if (rConditionDofList.size() != nr_dofs) rConditionDofList.resize(nr_dofs);

    // first node
    //rConditionDofList.push_back(geom[0].pGetDof(DISPLACEMENT_X));
    rConditionDofList[0] = geom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[1] = geom[0].pGetDof(DISPLACEMENT_Y);
    if (dimension == 3){}
    // second node
    rConditionDofList[2] = geom[1].pGetDof(DISPLACEMENT_X);
    rConditionDofList[3] = geom[1].pGetDof(DISPLACEMENT_Y);
    if (dimension == 3){}
    // lagrangian node
    rConditionDofList[4] = geom[2].pGetDof(LAGRANGE_MULTIPLIER_1);
    rConditionDofList[5] = geom[2].pGetDof(LAGRANGE_MULTIPLIER_2);
    rConditionDofList[6] = geom[2].pGetDof(LAGRANGE_MULTIPLIER_3);
    if (dimension == 3){}
}

void MinimalKineticCondition2D::GetValuesVector(Vector& Values, int Step)
{
/*    PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
    const unsigned int BlockSize = rPeriodicVariables.size();
    const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs

    if (Values.size() != LocalSize)
        Values.resize(LocalSize,false);

    unsigned int LocalIndex = 0;

    for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar = rPeriodicVariables.DoubleVariablesBegin();
            itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
    {
        Values[LocalIndex] = this->GetGeometry()[0].FastGetSolutionStepValue(*itDVar,Step);
        Values[LocalIndex+BlockSize] = this->GetGeometry()[1].FastGetSolutionStepValue(*itDVar,Step);
        ++LocalIndex;
    }

    for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar = rPeriodicVariables.VariableComponentsBegin();
            itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
    {
        Values[LocalIndex] = this->GetGeometry()[0].FastGetSolutionStepValue(*itCVar,Step);
        Values[LocalIndex+BlockSize] = this->GetGeometry()[1].FastGetSolutionStepValue(*itCVar,Step);
        ++LocalIndex;
    }*/
}

void MinimalKineticCondition2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
}

void MinimalKineticCondition2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
}

}

