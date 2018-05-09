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

#if !defined(KRATOS_SMALL_DISPLACEMENT_CUSTOM_H_INCLUDED )
#define KRATOS_SMALL_DISPLACEMENT_CUSTOM_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SmallDisplacementCustom
 * @ingroup MultiscaleROMApplication
 * @brief Small displacement element for 2D and 3D geometries.
 * @details Auxiliary small displacement element, to add custom methods.
 * @author Marcelo Raschi
 */

class KRATOS_API(MULTISCALE_ROM_APPLICATION) SmallDisplacementCustom
    : public SmallDisplacement
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of SmallDisplacementCustom
    KRATOS_CLASS_POINTER_DEFINITION(SmallDisplacementCustom);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallDisplacementCustom(IndexType NewId, GeometryType::Pointer pGeometry);
    SmallDisplacementCustom(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~SmallDisplacementCustom() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * @brief Returns the currently selected integration method
     * @return Current integration method selected
     * @todo ADD THE OTHER CREATE FUNCTION
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

	SmallDisplacementCustom() : SmallDisplacement()
    {
    }

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //SmallDisplacement& operator=(const SmallDisplacement& rOther);
    /// Copy constructor.
    //SmallDisplacement(const SmallDisplacement& rOther);
    ///@}

}; // Class SmallDisplacementCustom

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_CUSTOM_H_INCLUDED defined
