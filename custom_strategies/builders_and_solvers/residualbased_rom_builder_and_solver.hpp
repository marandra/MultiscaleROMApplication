#if !defined(KRATOS_ROM_BUILDER_AND_SOLVER )
#define  KRATOS_ROM_BUILDER_AND_SOLVER


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedROMBuilderAndSolver
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedROMBuilderAndSolver);

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {

        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }

    /** Destructor.
     */
    virtual ~ResidualBasedROMBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    //**************************************************************************
    //**************************************************************************

    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY
			if (!pScheme)
				KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

		const int nelements = static_cast<int>(r_model_part.Elements().size());

		ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
		ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
		ModelPart::ConditionsContainerType::iterator cond_begin = r_model_part.ConditionsBegin();

		LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
		LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
		Element::EquationIdVectorType EquationId;

		for (int k = 0; k < nelements; k++)
		{
			ModelPart::ElementsContainerType::iterator it = el_begin + k;

				pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);


                A += LHS_Contribution;
                b += RHS_Contribution;

				pScheme->CleanMemory(*(it.base()));

		}

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY;

		if (this->GetEchoLevel() > 1 && r_model_part.GetCommunicator().MyPID() == 0) {
			std::cout << "Setting up the dofs" << std::endl;
		}

        BaseType::mDofSet = DofsArrayType();

        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo
    )
    {
        std::size_t number_of_modes = CurrentProcessInfo[NUMBER_REDUCED_MODES];

        KRATOS_TRY

        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(number_of_modes, number_of_modes, false);
            ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
        }
        else
        {
            if (A.size1() != number_of_modes || A.size2() != number_of_modes)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                A.resize(number_of_modes, number_of_modes, true);
                ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
            }
        }
        if (Dx.size() != number_of_modes)
            Dx.resize(number_of_modes, false);
        if (b.size() != number_of_modes)
            b.resize(number_of_modes, false);

        //TODO fix or remove this block
        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
        }

        KRATOS_CATCH("")

    }


protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    

    //**************************************************************************
	virtual void ConstructMatrixStructure(
		TSystemMatrixType& A,
		ElementsContainerType& rElements,
		ConditionsArrayType& rConditions,
		ProcessInfo& CurrentProcessInfo)
	{
		//filling with zero the matrix (creating the structure)
		Timer::Start("MatrixStructure");

		const std::size_t equation_size = CurrentProcessInfo[NUMBER_REDUCED_MODES];

#ifdef USE_GOOGLE_HASH
		std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
		const std::size_t empty_key = 2 * equation_size + 10;
#else
		std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
		for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
		{
#ifdef USE_GOOGLE_HASH
			indices[iii].set_empty_key(empty_key);
#else
			indices[iii].reserve(40);
#endif
		}

		Element::EquationIdVectorType ids(3, 0);

		const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel for firstprivate(nelements, ids)
		for (int iii = 0; iii<nelements; iii++)
		{
			typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
			(i_element)->EquationIdVector(ids, CurrentProcessInfo);

			for (std::size_t i = 0; i < ids.size(); i++)
			{
				if (ids[i] < equation_size)
				{
#ifdef _OPENMP
                    //TODO remove comments
                    // omp_set_lock(&mlock_array[ids[i]]);
#endif
					auto& row_indices = indices[ids[i]];
					for (auto it = ids.begin(); it != ids.end(); it++)
					{
						if (*it < equation_size)
							row_indices.insert(*it);
					}
#ifdef _OPENMP
                    //TODO remove comments
				//	omp_unset_lock(&mlock_array[ids[i]]);
#endif
				}
			}

		}

		const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
		for (int iii = 0; iii<nconditions; iii++)
		{
			typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
			(i_condition)->EquationIdVector(ids, CurrentProcessInfo);
			for (std::size_t i = 0; i < ids.size(); i++)
			{
				if (ids[i] < equation_size)
				{
#ifdef _OPENMP
					omp_set_lock(&mlock_array[ids[i]]);
#endif
					auto& row_indices = indices[ids[i]];
					for (auto it = ids.begin(); it != ids.end(); it++)
					{
						if (*it < equation_size)
							row_indices.insert(*it);
					}
#ifdef _OPENMP
					omp_unset_lock(&mlock_array[ids[i]]);
#endif
				}
			}
		}

		//count the row sizes
		unsigned int nnz = 0;
		for (unsigned int i = 0; i < indices.size(); i++)
			nnz += indices[i].size();

		A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

		double* Avalues = A.value_data().begin();
		std::size_t* Arow_indices = A.index1_data().begin();
		std::size_t* Acol_indices = A.index2_data().begin();

		//filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
		Arow_indices[0] = 0;
		for (int i = 0; i < static_cast<int>(A.size1()); i++)
			Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();



#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(A.size1()); i++)
		{
			const unsigned int row_begin = Arow_indices[i];
			const unsigned int row_end = Arow_indices[i + 1];
			unsigned int k = row_begin;
			for (auto it = indices[i].begin(); it != indices[i].end(); it++)
			{
				Acol_indices[k] = *it;
				Avalues[k] = 0.0;
				k++;
			}

			std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

		}

		A.set_filled(indices.size() + 1, nnz);

		Timer::Stop("MatrixStructure");
	}


private:

#ifdef _OPENMP
	std::vector< omp_lock_t > mlock_array;
#endif

}; /* Class ResidualBasedROMBuilderAndSolver */



} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

