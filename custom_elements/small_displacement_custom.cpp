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
