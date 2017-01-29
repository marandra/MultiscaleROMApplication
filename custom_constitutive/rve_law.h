#if !defined (KRATOS_RVE_LAW_H_INCLUDED)
#define  KRATOS_RVE_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{

template<class TStrategyType >
class KRATOS_API(MULTISCALE_ROM_APPLICATION) RVELaw : public ConstitutiveLaw
{
protected:

    typename TStrategyType::Pointer mpSolvingStrategy;

public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /**
     * Counted pointer of RVELaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( RVELaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    RVELaw(){};
    
    RVELaw(typename TStrategyType::Pointer pSolvingStrategy)
        : mpSolvingStrategy(pSolvingStrategy)
    {
        //check if we can correctly get the model part from the strategy
        ModelPart& mr_model_part = mpSolvingStrategy->GetModelPart();
        std::cout << mr_model_part << std::endl;
    };

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
     * Operations needed by the base class:
     */



protected:



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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }


    ///@}

}; // Class RVELaw
}  // namespace Kratos.
#endif // KRATOS_RVE_LAW_H_INCLUDED  defined 
