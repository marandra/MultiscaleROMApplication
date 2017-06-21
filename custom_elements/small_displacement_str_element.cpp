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

    void SmallDisplacementStrElement::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Matrix F( dim, dim );
        Matrix D( strain_size, strain_size );
        Vector strain_vector( strain_size );
        Vector stress_vector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix B( strain_size, number_of_nodes * dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        Values.SetStrainVector(strain_vector); //this is the input  parameter
        Values.SetStressVector(stress_vector); //this is the output parameter

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        if ( rOutput.size() != integration_points.size() )
        {
            rOutput.resize( integration_points.size() );
        }

        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod());

            //Compute B and strain // TODO: MOVE THIS TO THE CL!!!!!!
            this->CalculateB( B, DN_DX, integration_points, point_number );
            noalias(strain_vector) = prod(B, displacements);
            F = ComputeEquivalentF(strain_vector);

            if ( rVariable == CAUCHY_STRESS_VECTOR)
            {
                if ( rOutput[point_number].size() != stress_vector.size() )
                {
                    rOutput[point_number].resize( stress_vector.size(), false );
                }

                // Here we essentially set the input parameters
                const double detF = MathUtils<double>::Det(F);
                Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
                Values.SetDeformationGradientF(F); //F computed somewhere else
                Values.SetConstitutiveMatrix(D); //assuming the determinant is computed somewhere else
                //Values.SetStressVector(stress_vector); //F computed somewhere else

                //actually do the computations in the ConstitutiveLaw
                mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy(Values); //here the calculations are actually done

                for ( unsigned int ii = 0; ii < stress_vector.size(); ii++ )
                {
                    rOutput[point_number]( ii ) = stress_vector[ii];
                }
            }
        }

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


