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
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
    SmallDisplacementStrElement::SmallDisplacementStrElement( IndexType NewId, GeometryType::Pointer pGeometry )
            : KinematicLinear( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    SmallDisplacementStrElement::SmallDisplacementStrElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : KinematicLinear( NewId, pGeometry, pProperties )
    {
    }

    Element::Pointer SmallDisplacementStrElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new SmallDisplacementStrElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    SmallDisplacementStrElement::~SmallDisplacementStrElement()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void SmallDisplacementStrElement::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "SmallDisplacementStrElement" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, KinematicLinear );
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void SmallDisplacementStrElement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, KinematicLinear );
    }

} // Namespace Kratos


