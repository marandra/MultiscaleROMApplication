// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//

// System includes

// External includes

// Project includes
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

    //if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
    //{
    //}
    //else
    //{
    //    SmallDisplacement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    //}
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementStrElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacement);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementStrElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacement);
}

} // Namespace Kratos
