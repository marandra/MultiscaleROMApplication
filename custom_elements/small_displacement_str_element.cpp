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
#include "custom_elements/small_displacement_str_element.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "multiscale_rom_application_variables.h"
#include "utilities/math_utils.h"

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

void SmallDisplacementStrElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                               std::vector<Matrix>& rOutput,
                                                               const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rOutput.size() != GetGeometry().IntegrationPoints().size())
    {
        rOutput.resize(GetGeometry().IntegrationPoints().size());
    }

    if (rVariable == REDUCED_MODES_MATRIX)
    {
        // Create and initialize element variables:
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints();

        // Reading integration points
        for (unsigned int point_number = 0;
             point_number < mConstitutiveLawVector.size(); point_number++)
        {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number,
                                        integration_points);

            if (rOutput[point_number].size2() != this_kinematic_variables.B.size2())
            {
                rOutput[point_number].resize(this_kinematic_variables.B.size1(),
                                             this_kinematic_variables.B.size2(), false);
            }

            rOutput[point_number] = this_kinematic_variables.B;
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
