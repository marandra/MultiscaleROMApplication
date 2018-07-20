// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//

//  Main authors:    Marcelo Raschi
//


// System includes

// External includes

// Project includes
#include "custom_elements/small_displacement_custom.hpp"

namespace Kratos
{
SmallDisplacementCustom::SmallDisplacementCustom(IndexType NewId, GeometryType::Pointer pGeometry)
    : SmallDisplacement(NewId, pGeometry)
{
}

//************************************************************************************
//************************************************************************************

SmallDisplacementCustom::SmallDisplacementCustom(IndexType NewId,
                                                 GeometryType::Pointer pGeometry,
                                                 PropertiesType::Pointer pProperties)
    : SmallDisplacement(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************

Element::Pointer SmallDisplacementCustom::Create(IndexType NewId,
                                                 NodesArrayType const& ThisNodes,
                                                 PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<SmallDisplacementCustom>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

SmallDisplacementCustom::~SmallDisplacementCustom()
{
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( GetProperties(),
            GetGeometry(),
            row( GetGeometry().ShapeFunctionsValues(  ), point_number )
            );
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementCustom::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    if ( rOutput.size() != GetGeometry().IntegrationPoints(  ).size() )
        rOutput.resize( GetGeometry().IntegrationPoints(  ).size() );

    if  (rVariable == DAMAGE_VARIABLE) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        for (unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number)
        {
            double damage = 0.0;
            mConstitutiveLawVector[point_number]->CalculateValue(Values, DAMAGE_VARIABLE, damage);
            rOutput[point_number] = damage;
        }
    } else {
        SmallDisplacement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    )
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ), Matrix(prod(D, B)));
        KRATOS_WATCH("DEBUG CUSTOM ELEMENT: OVERWRITE B")

    KRATOS_CATCH( "" )
}
/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    KRATOS_WATCH("DEBUG CUSTOM ELEMENT: CALCULATE ALL")

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
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_WATCH("DEBUG ******* MATRIX OUTPUT")
    KRATOS_WATCH(rLeftHandSideMatrix);
    KRATOS_WATCH(this_kinematic_variables.B)
    KRATOS_WATCH("DEBUG ******* ")
    KRATOS_WATCH("DEBUG ******* ")
    KRATOS_WATCH("DEBUG ******* ")

    KRATOS_CATCH( "" )
}
/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementCustom::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const Vector& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    )
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementCustom::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacement);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementCustom::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacement);
}

} // Namespace Kratos
