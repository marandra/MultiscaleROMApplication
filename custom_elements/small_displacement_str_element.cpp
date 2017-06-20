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

    void SmallDisplacementStrElement::CalculateOnIntegrationPoints(
        const Variable<double >& rVariable,
        std::vector< double >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        if (rVariable == GAUSS_WEIGHTS)
        {
            const unsigned int number_of_nodes = GetGeometry().size();
            const unsigned int dim = GetGeometry().WorkingSpaceDimension();
            const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

            Matrix DN_DX( number_of_nodes, dim );
            Matrix J0(dim,dim), InvJ0(dim,dim);

            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
            if ( rOutput.size() != integration_points.size() )
            {
                rOutput.resize( integration_points.size() );
            }

            for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
            {
                const double detJ = CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod());

                //calculating weights for integration on the "reference configuration"
                double integration_weight = integration_points[point_number].Weight() * detJ;
                if( dim == 2 && this->GetProperties().Has( THICKNESS ) )
                {
                        integration_weight *= this->GetProperties()[THICKNESS];
                }
                rOutput[point_number] = integration_weight;
             }
        }

        //KinematicLinear::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo)
        BaseSolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

        KRATOS_CATCH( "" )
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


