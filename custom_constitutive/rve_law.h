#if !defined(KRATOS_RVE_IDENTIFIER_LAW_H_INCLUDED)
#define KRATOS_RVE_IDENTIFIER_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
class KRATOS_API(MULTISCALE_ROM_APPLICATION) RVELaw : public ConstitutiveLaw
{
protected:

public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    /**
     * Counted pointer of RVELaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(RVELaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    //    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    //    RVELaw (const RVELaw& rOther);

    /**
     * Assignment operator.
     */

    /**
     * Destructor.
     */
    virtual ~RVELaw(){};

    /**
     * Operators
     */
    /**
    * This is to be called at the very beginning of the calculation
    * (e.g. from InitializeElement) in order to initialize all relevant
    * attributes of the constitutive law
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    */
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);

    /**
     * Operations needed by the base class:
     */

protected:
private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{
    std::vector<ConstitutiveLaw::Pointer> mCL_list;
    std::vector<double> mIW_list;
    std::vector<Matrix> mB_list;
    std::vector<Properties::Pointer> mprop_list;
    ModelPart::Pointer mpRVEModelPart;
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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    ///@}

}; // Class RVELaw
} // namespace Kratos.
#endif // KRATOS_RVE_LAW_H_INCLUDED  defined
