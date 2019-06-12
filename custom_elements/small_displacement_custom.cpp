// System includes

// External includes

// Project includes
#include "custom_elements/small_displacement_custom.h"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
SmallDisplacementCustom::SmallDisplacementCustom(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseSolidElement(NewId, pGeometry)
{
}

//************************************************************************************
//************************************************************************************

SmallDisplacementCustom::SmallDisplacementCustom(IndexType NewId,
                                                 GeometryType::Pointer pGeometry,
                                                 PropertiesType::Pointer pProperties)
    : BaseSolidElement(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************

Element::Pointer SmallDisplacementCustom::Create(IndexType NewId,
                                                 NodesArrayType const& ThisNodes,
                                                 PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementCustom>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementCustom::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementCustom>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementCustom::~SmallDisplacementCustom()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementCustom::UseElementProvidedStrain() const
{
    return true;
}


/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        // Contribution to external forces
        const Vector body_force = this->GetBodyForce(integration_points, point_number);

        const Matrix& r_global_modes = rCurrentProcessInfo.GetValue(GLOBAL_MODES_MATRIX);
        const int& r_mode_index = rCurrentProcessInfo.GetValue(MODE_INDEX);
        const IndexType element_index = this->Id();
        const int nr_components = this_constitutive_variables.StressVector.size();
        const int global_index = (element_index - 1) * nr_components * integration_points.size() + point_number * nr_components;
        Vector ip_mode(nr_components);
        for (int i_component = 0; i_component < nr_components; ++i_component) {
            ip_mode[i_component] = r_global_modes(global_index + i_component, r_mode_index);
        }
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        double int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, ip_mode, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY
    Matrix identity = IdentityMatrix(6);
    //noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod(trans(B), B);
    noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ), Matrix(prod(identity, B)));

    KRATOS_CATCH( "" )
}
/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& r_ip_mode,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    noalias( rRightHandSideVector ) += IntegrationWeight * prod( trans( rThisKinematicVariables.B ), r_ip_mode);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        GetGeometry().IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    Vector displacements;
    GetValuesVector(displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, displacements);
    rThisKinematicVariables.F = ComputeEquivalentF(strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Displacements vector
    Vector displacements;
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*2     ) = rDN_DX( i, 0 );
            rB( 1, i*2 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*2     ) = rDN_DX( i, 1 );
            rB( 2, i*2 + 1 ) = rDN_DX( i, 0 );
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
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

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

Matrix SmallDisplacementCustom::ComputeEquivalentF(const Vector& rStrainTensor) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    Matrix F(dim,dim);

    if(dim == 2) {
        F(0,0) = 1.0+rStrainTensor(0);
        F(0,1) = 0.5*rStrainTensor(2);
        F(1,0) = 0.5*rStrainTensor(2);
        F(1,1) = 1.0+rStrainTensor(1);
    } else {
        F(0,0) = 1.0+rStrainTensor(0);
        F(0,1) = 0.5*rStrainTensor(3);
        F(0,2) = 0.5*rStrainTensor(5);
        F(1,0) = 0.5*rStrainTensor(3);
        F(1,1) = 1.0+rStrainTensor(1);
        F(1,2) = 0.5*rStrainTensor(4);
        F(2,0) = 0.5*rStrainTensor(5);
        F(2,1) = 0.5*rStrainTensor(4);
        F(2,2) = 1.0+rStrainTensor(2);
    }

    return F;
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementCustom::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = BaseSolidElement::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementCustom::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseSolidElement);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementCustom::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseSolidElement);
}

} // Namespace Kratos
