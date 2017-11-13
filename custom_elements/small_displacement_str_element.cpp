// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_str_element.h"

namespace Kratos
{
SmallDisplacementStrElement::SmallDisplacementStrElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry)
    : SmallDisplacement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SmallDisplacementStrElement::SmallDisplacementStrElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry,
                                                         PropertiesType::Pointer pProperties)
    : SmallDisplacement(NewId, pGeometry, pProperties)
{
}

Element::Pointer SmallDisplacementStrElement::Create(IndexType NewId,
                                                     NodesArrayType const& ThisNodes,
                                                     PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallDisplacementStrElement(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
}

SmallDisplacementStrElement::~SmallDisplacementStrElement()
{
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementStrElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if ( rOutput.size() != GetGeometry().IntegrationPoints(  ).size() )
    {
        rOutput.resize( GetGeometry().IntegrationPoints(  ).size() );
    }

    if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
    {
        // Create and initialize element variables:
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        // Displacements vector
        Vector displacements;
        GetValuesVector(displacements);

        //reading integration points
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            if ( rOutput[point_number].size() != strain_size)
            {
                rOutput[point_number].resize( strain_size, false );
            }

            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_points);

            rOutput[point_number] = prod(this_kinematic_variables.B, displacements);
        }
    }
    else
    {
        SmallDisplacement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementStrElement::save(Serializer& rSerializer) const
{
    rSerializer.save("Name", "SmallDisplacementStrElement");
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacement);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementStrElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacement);
}

} // Namespace Kratos
