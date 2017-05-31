#if !defined(KRATOS_RESIDUALBASED_INCREMENTAL_ROM_STATIC_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUALBASED_INCREMENTAL_ROM_STATIC_SCHEME_H_INCLUDED

#include <iostream>
#include <string>

#include "includes/define.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

namespace Kratos
{
// TODO Description of the scheme here
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedIncrementalROMStaticScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    /// Pointer definition of ResidualBasedIncrementalROMStaticScheme
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalROMStaticScheme);
    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> BaseType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /// Default constructor.
    ResidualBasedIncrementalROMStaticScheme()
    {
    }

    /// Destructor.
    virtual ~ResidualBasedIncrementalROMStaticScheme()
    {
    }

    virtual void Initialize(ModelPart& r_model_part)
    {
        BaseType::Initialize(r_model_part);
        std::size_t nr_modes =
            r_model_part.GetProcessInfo()[REDUCED_MODES_WEIGHTS].size();
        TSystemVectorType modes_weights = ZeroVector(nr_modes);
        r_model_part.GetProcessInfo().SetValue(REDUCED_MODES_WEIGHTS, modes_weights);
    }

    virtual void Predict(
            ModelPart& r_model_part,
            DofsArrayType& rDofSet,
            LocalSystemMatrixType & A,
            LocalSystemVectorType & Dx,
            LocalSystemVectorType & b
    )
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    virtual void Update(ModelPart& r_model_part,
                        DofsArrayType& rDofSet,
                        TSystemMatrixType& A,
                        TSystemVectorType& Dx,
                        TSystemVectorType& b)
    {
        TSystemVectorType& modes_weights =
            r_model_part.GetProcessInfo().GetValue(REDUCED_MODES_WEIGHTS);
        noalias(modes_weights) += Dx;
        r_model_part.GetProcessInfo().SetValue(REDUCED_MODES_WEIGHTS, modes_weights);
    }

protected:
private:
    /// Assignment operator.
    ResidualBasedIncrementalROMStaticScheme& operator=(ResidualBasedIncrementalROMStaticScheme const& rOther)
    {
    }

    /// Copy constructor.
    ResidualBasedIncrementalROMStaticScheme(ResidualBasedIncrementalROMStaticScheme const& rOther)
    {
    }

}; // Class ResidualBasedIncrementalROMStaticScheme
} // namespace Kratos.

#endif // KRATOS_RESIDUALBASED_INCREMENTAL_AITKEN_STATIC_SCHEME_H_INCLUDED
       // defined
