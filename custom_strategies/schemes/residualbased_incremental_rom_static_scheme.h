#if !defined(KRATOS_RESIDUALBASED_INCREMENTAL_ROM_STATIC_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUALBASED_INCREMENTAL_ROM_STATIC_SCHEME_H_INCLUDED

#include "boost/smart_ptr.hpp"

#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

namespace Kratos
{
// TODO Description of the scheme here
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedIncrementalROMStaticScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalROMStaticScheme);
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef ModelPart::ElementsContainerType             ElementsArrayType;
    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    /** Constructor. */
    ResidualBasedIncrementalROMStaticScheme()
        : Scheme<TSparseSpace,TDenseSpace>()
    {}

    /** Destructor. */
    virtual ~ResidualBasedIncrementalROMStaticScheme() {}

    virtual void Initialize(ModelPart& r_model_part) override
    {
        BaseType::Initialize(r_model_part);
        std::size_t nr_modes =
            r_model_part.GetProcessInfo()[REDUCED_MODES_WEIGHTS].size();
        TSystemVectorType modes_weights = ZeroVector(nr_modes);
        r_model_part.GetProcessInfo().SetValue(REDUCED_MODES_WEIGHTS, modes_weights);
    }


    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        TSystemVectorType& modes_weights =
            r_model_part.GetProcessInfo().GetValue(REDUCED_MODES_WEIGHTS);
        noalias(modes_weights) += Dx;
        r_model_part.GetProcessInfo().SetValue(REDUCED_MODES_WEIGHTS, modes_weights);

        KRATOS_CATCH("")
    }


    virtual void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }


    /**
     * It initializes a non-linear iteration (for the element)
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY;

        // Initializes the non-linear iteration for all the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd   = rElements.begin() + ElementPartition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->InitializeNonLinearIteration(CurrentProcessInfo);
            }
        }
        
        // Initializes the non-linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();
        
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd   = rConditions.begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                itCond->InitializeNonLinearIteration(CurrentProcessInfo);
            }
        }

        KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont: The condition to compute
     * @param CurrentProcessInfo: The current process info instance
     */


    void InitializeNonLinearIteration(
        Condition::Pointer rCurrentCondition,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentElement: The element to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
        Element::Pointer rCurrentElement,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }


    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

    this function calculates at the same time the contribution to the LHS and to the RHS
    of the system
    */
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    //***************************************************************************
    //***************************************************************************
    virtual void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************
    virtual void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }


    /** functions totally analogous to the precedent but applied to
    the "condition" objects
    */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

protected:
private:

}; /* Class Scheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_STANDARD_STATIC_SCHEME  defined */

