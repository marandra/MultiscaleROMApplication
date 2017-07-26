//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/total_lagrangian_fbar_element.hpp"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianFbarElement::TotalLagrangianFbarElement(IndexType NewId,
                                                       GeometryType::Pointer pGeometry)
    : LargeDisplacementElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianFbarElement::TotalLagrangianFbarElement(IndexType NewId,
                                                       GeometryType::Pointer pGeometry,
                                                       PropertiesType::Pointer pProperties)
    : LargeDisplacementElement(NewId, pGeometry, pProperties)
{
    // DO NOT ADD DOFS HERE!!!
}

//******************************COPY
//CONSTRUCTOR**************************************
//************************************************************************************

TotalLagrangianFbarElement::TotalLagrangianFbarElement(TotalLagrangianFbarElement const& rOther)
    : LargeDisplacementElement(rOther),
      mTotalDomainInitialSize(rOther.mTotalDomainInitialSize),
      mInvJ0(rOther.mInvJ0),
      mDetJ0(rOther.mDetJ0)
{
}

//*******************************ASSIGMENT
//OPERATOR***********************************
//************************************************************************************

TotalLagrangianFbarElement& TotalLagrangianFbarElement::operator=(TotalLagrangianFbarElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

    mInvJ0.clear();
    mInvJ0.resize(rOther.mInvJ0.size());

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mInvJ0[i] = rOther.mInvJ0[i];
    }

    mTotalDomainInitialSize = rOther.mTotalDomainInitialSize;
    mDetJ0 = rOther.mDetJ0;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer TotalLagrangianFbarElement::Create(IndexType NewId,
                                                    NodesArrayType const& rThisNodes,
                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TotalLagrangianFbarElement(
        NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer TotalLagrangianFbarElement::Clone(IndexType NewId,
                                                   NodesArrayType const& rThisNodes) const
{
    TotalLagrangianFbarElement NewElement(
        NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if (NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
    {
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

        if (NewElement.mConstitutiveLawVector.size() !=
            NewElement.GetGeometry().IntegrationPointsNumber())
            KRATOS_THROW_ERROR(std::logic_error,
                               "constitutive law not has the correct size ",
                               NewElement.mConstitutiveLawVector.size())
    }

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    //-----------//

    if (NewElement.mInvJ0.size() != mInvJ0.size())
        NewElement.mInvJ0.resize(mInvJ0.size());

    for (unsigned int i = 0; i < mInvJ0.size(); i++)
    {
        NewElement.mInvJ0[i] = mInvJ0[i];
    }

    NewElement.mTotalDomainInitialSize = mTotalDomainInitialSize;
    NewElement.mDetJ0 = mDetJ0;

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TotalLagrangianFbarElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianFbarElement::~TotalLagrangianFbarElement()
{
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementElement::Initialize();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    // Resizing jacobian inverses container
    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(), false);

    // Compute jacobian inverses and set the domain initial size:
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
    mTotalDomainInitialSize = 0.00;

    // calculating the inverse J0

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        // calculating and storing inverse of the jacobian and the parameters
        // needed
        MathUtils<double>::InvertMatrix(J0[PointNumber], mInvJ0[PointNumber],
                                        mDetJ0[PointNumber]);

        // calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    KRATOS_CATCH("")
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                                    GeneralVariables& rVariables,
                                                    double& rIntegrationWeight)
{
    KRATOS_WATCH("DEBUG ENTRA LHS2")
    KRATOS_TRY

    // int i = 0; //  bucle infinito
    // while (i == 0){
    //	i = i;
    //}

    // contributions of the stiffness matrix calculated on the reference
    // configuration
    if (rLocalSystem.CalculationFlags.Is(LargeDisplacementElement::COMPUTE_LHS_MATRIX_WITH_COMPONENTS))
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

            if (rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX)
            {
                // operation performed: add Kg to the rLefsHandSideMatrix
                this->CalculateAndAddKuug(rLeftHandSideMatrices[i], rVariables,
                                          rIntegrationWeight);
                calculated = true;
            }

            if (calculated == false)
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

        // operation performed: add Kg to the rLefsHandSideMatrix
        this->CalculateAndAddKuug(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }

    // KRATOS_WATCH( rLeftHandSideMatrix )
    // KRATOS_WATCH("DEBUG SALE LHS")
    KRATOS_CATCH("")
}

//************************************************************************************
//****************************************** JLM calcula la matriz de rigidez
//material

void TotalLagrangianFbarElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                                     GeneralVariables& rVariables,
                                                     double& rIntegrationWeight)
{
    KRATOS_TRY

    // noalias(rLeftHandSideMatrix) += prod(trans(rVariables.B),
    // rIntegrationWeight * Matrix(prod(rVariables.ConstitutiveMatrix,
    // rVariables.B))); //to be optimized to remove the temporary

    // Kmar = fact*B'*C*B
    double factor =
        std::sqrt(rVariables.detF0 / rVariables.detFT); // only for plain strain

    // contributions to stiffness matrix calculated on the reference config
    noalias(rLeftHandSideMatrix) +=
        factor * prod(trans(rVariables.B),
                      rIntegrationWeight *
                          Matrix(prod(rVariables.ConstitutiveMatrix, rVariables.B))); // to be optimized to remove the temporary

    // std::cout<<std::endl;
    // std::cout<<" Kmat "<<rLeftHandSideMatrix<<std::endl;
    // KRATOS_WATCH("DEBUG KUUM")

    // KRATOS_WATCH(rLeftHandSideMatrix)
    KRATOS_CATCH("")
}

//************************************************************************************
//***************************** JLM Modifico la forma en que se calcula Kgeo
//para FBAR

void TotalLagrangianFbarElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                                     GeneralVariables& rVariables,
                                                     double& rIntegrationWeight)

{
    KRATOS_TRY

    // Kgeom = fact*B'*C*F*(B'*F0^{-1} - B'*F^{-1})'/pot

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int element_size = number_of_nodes * dimension;
    unsigned int voigt_size = 3; // ver si no hace falta colocar la fila z

    if (dimension == 3)
    {
        KRATOS_THROW_ERROR(
            std::invalid_argument,
            "something is wrong with the dimension, solo esta hecho para 2D!!!",
            ""); //	voigt_size = 6;
    }

    double factor = 0.5; // only for plain strain

    Matrix InvF0;
    Matrix InvF;
    // InvF0.resize(dimension, dimension);
    // InvF.resize(dimension, dimension);

    MathUtils<double>::InvertMatrix(rVariables.F0, InvF0, rVariables.detF0);
    MathUtils<double>::InvertMatrix(rVariables.F, InvF, rVariables.detF);

    Vector v_InvF0 = MathUtils<double>::SymmetricTensorToVector(InvF0);
    Vector v_InvF = MathUtils<double>::SymmetricTensorToVector(InvF);
    Vector v_F = MathUtils<double>::SymmetricTensorToVector(rVariables.F);

    // Vector InvBF = factor * prod( trans(rVariables.B), Vector v_InvF0 ) ; //
    // JLM ver la resta de matrices
    Vector v_InvBF;
    v_InvBF.resize(element_size);

    v_InvF0.resize(voigt_size);
    v_InvF.resize(voigt_size);
    noalias(v_InvBF) = factor * Vector(prod(trans(rVariables.B), (v_InvF0 - v_InvF)));
    // KRATOS_WATCH(v_InvF0)
    //	KRATOS_WATCH(v_InvF)
    Matrix m_BFB;
    m_BFB.resize(voigt_size, element_size);
    noalias(m_BFB) = ZeroMatrix(voigt_size, element_size);

    for (unsigned int idof = 0; idof < element_size;
         idof++) // prod(v_F, trans(v_InvBF)) = Matrix
    {
        for (unsigned int ivoigt = 0; ivoigt < voigt_size; ivoigt++)
        {
            m_BFB(ivoigt, idof) = v_F(ivoigt) * v_InvBF(idof);
        }
    }

    // KRATOS_WATCH(v_F);
    // KRATOS_WATCH(v_InvBF);
    // KRATOS_WATCH(m_BFB);

    // int i = 0; //  bucle infinito
    // while (i == 0) {
    //	i = i;
    //}
    // Vector C(dimension);
    // noalias(C) = prod(trans(rVariables.B), v_InvF0);
    // MatrixType m_FBF = prod(trans(v_InvBF), v_F);
    // fact*B'*C*F*(B'*F0^{ -1 } -B'*F^{-1})' / pot
    // Matrix ReducedKg = prod(trans(rVariables.B), rIntegrationWeight *
    // Matrix(prod(rVariables.ConstitutiveMatrix, m_BFB ))); //to be optimized
    noalias(rLeftHandSideMatrix) += prod(
        trans(rVariables.B),
        rIntegrationWeight * Matrix(prod(rVariables.ConstitutiveMatrix, m_BFB))); // to be optimized
    // Matrix ReducedKg = ZeroMatrix(8, 8);
    // ReducedKg *= 0. ;
    // KRATOS_WATCH(ReducedKg.size2())

    // KRATOS_WATCH(ReducedKg.size1())
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH(ReducedKg)
    // MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix,
    // ReducedKg, dimension );

    // std::cout<<std::endl;
    // std::cout<<" Kmat + Kgeo "<<rLeftHandSideMatrix<<std::endl;

    // KRATOS_WATCH(ReducedKg.size2())

    // KRATOS_WATCH(ReducedKg)
    // int i = 0; //  bucle infinito
    // while (i == 0) {
    //	i = i;
    //}
    KRATOS_CATCH("")
}

//*********************************COMPUTE
//KINEMATICS*********************************
//************************************************************************************

void TotalLagrangianFbarElement::CalculateKinematics(GeneralVariables& rVariables,
                                                     const double& rPointNumber)

{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De =
        rVariables.GetShapeFunctionsGradients();

    // Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    // Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    // Jacobian Determinant for the isoparametric and numerical integration
    //
    rVariables.detJ = mDetJ0[rPointNumber];

    // Calculating the cartesian derivatives [dN/dx_n] = [dN/d£][d£/dx_0]
    noalias(rVariables.DN_DX) = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

    // Deformation Gradient F [dx_n+1/dx_0] = [dx_n+1/d£] [d£/dx_0]
    noalias(rVariables.F) = prod(rVariables.j[rPointNumber], mInvJ0[rPointNumber]);

    // Determinant of the deformation gradient F
    rVariables.detF = MathUtils<double>::Det(rVariables.F);

    //
    //
    //

    //
    //

    // Determinant of the Deformation Gradient F0
    // (in this element F = F0, then F0 is set to the identity for coherence in
    // the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0 = identity_matrix<double>(dimension);

    // Set Shape Functions Values for this integration point
    rVariables.N = row(Ncontainer, rPointNumber);

    // Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::CalculateDeformationMatrix(Matrix& rB, Matrix& rF, Matrix& rDN_DX)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear(); // set all components to zero

    if (dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 2 * i;

            rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
            rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
            rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
            rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
            rB(2, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
            rB(2, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        }
    }
    else if (dimension == 3)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 3 * i;

            rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
            rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
            rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
            rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
            rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
            rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
            rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
            rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
            rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
            rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
            rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
            rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
            rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
            rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
            rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
            rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
            rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
            rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
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
//MASS****************************
//************************************************************************************

double& TotalLagrangianFbarElement::CalculateTotalMass(double& rTotalMass,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rTotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

    if (dimension == 2)
        rTotalMass *= GetProperties()[THICKNESS];

    return rTotalMass;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                              std::vector<double>& rOutput,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int& integration_points_number =
        GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

    if (rOutput.size() != integration_points_number)
        rOutput.resize(integration_points_number);

    if (rVariable == STRAIN_ENERGY)
    {
        double Thickness = 1.0;
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

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

        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);

        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            // compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables, PointNumber);
            // to take in account previous step writing
            if (mFinalizedStep)
            {
                this->GetHistoricalVariables(Variables, PointNumber);
            }
            // set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables, Values, PointNumber);

            double StrainEnergy = 0.0;

            // compute stresses and constitutive parameters
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);
            mConstitutiveLawVector[PointNumber]->GetValue(STRAIN_ENERGY, StrainEnergy);

            rOutput[PointNumber] =
                Variables.detJ * integration_points[PointNumber].Weight() *
                Thickness * StrainEnergy; // 1/2 * sigma * epsilon

        } // for each gauss_point
    }
    else
    {
        LargeDisplacementElement::CalculateOnIntegrationPoints(
            rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::GetHistoricalVariables(GeneralVariables& rVariables,
                                                        const double& rPointNumber)
{
    // LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    // //Deformation Gradient F [dx_n/dx_0] = [dx_n/d£] [d£/dx_0]
    // noalias( rVariables.F0 ) = prod( rVariables.j[rPointNumber],
    // mInvJ0[rPointNumber] );

    // //Deformation Gradient F0
    // rVariables.detF0 = MathUtils<double>::Det(rVariables.F0);
}

//************************************CALCULATE VOLUME
//CHANGE*************************
//************************************************************************************

double& TotalLagrangianFbarElement::CalculateVolumeChange(double& rVolumeChange,
                                                          GeneralVariables& rVariables)
{
    KRATOS_TRY

    rVolumeChange = 1.0;

    return rVolumeChange;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TotalLagrangianFbarElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LargeDisplacementElement)
    rSerializer.save("mTotalDomainInitialSize", mTotalDomainInitialSize);
    rSerializer.save("InvJ0", mInvJ0);
    rSerializer.save("DetJ0", mDetJ0);
}

void TotalLagrangianFbarElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LargeDisplacementElement)
    rSerializer.load("mTotalDomainInitialSize", mTotalDomainInitialSize);
    rSerializer.load("InvJ0", mInvJ0);
    rSerializer.load("DetJ0", mDetJ0);
}

} // Namespace Kratos
