//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

// Project includes
#include "add_custom_strategies_to_python.h"
#include "containers/flags.h"
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "custom_strategies/builders_and_solvers/residualbased_rom_builder_and_solver.hpp"
#include "custom_strategies/schemes/residualbased_incremental_rom_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;
//typedef UblasSpace<double, CompressedMatrix, Vector> DenseSparseSpaceType;
//typedef UblasSpace<double, Matrix, Vector> DenseSparseSpaceType;
//typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

void AddCustomStrategiesToPython()
{
    //typedef UblasSpace<double, CompressedMatrix, Vector> DenseSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // base types
    typedef LinearSolver<DenseSparseSpaceType, LocalSpaceType> DenseLinearSolverType;
    typedef BuilderAndSolver<DenseSparseSpaceType, LocalSpaceType, DenseLinearSolverType> BuilderAndSolverType;
    typedef Scheme<DenseSparseSpaceType, LocalSpaceType> BaseSchemeType;

    // custom builder_and_solver types
    typedef ResidualBasedROMBuilderAndSolver<DenseSparseSpaceType, LocalSpaceType, DenseLinearSolverType> ResidualBasedROMBuilderAndSolverType;
    typedef ResidualBasedIncrementalROMStaticScheme<DenseSparseSpaceType, LocalSpaceType> ResidualBasedIncrementalROMStaticSchemeType;
    typedef ResidualBasedNewtonRaphsonStrategy< DenseSparseSpaceType, LocalSpaceType, DenseLinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

    typedef ConvergenceCriteria< DenseSparseSpaceType, LocalSpaceType > DenseConvergenceCriteriaType;


    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************
    class_<DenseLinearSolverType, DenseLinearSolverType::Pointer, boost::noncopyable>(
        "DenseLinearSolver", init<>());


    class_<ResidualBasedROMBuilderAndSolverType,  ResidualBasedROMBuilderAndSolverType::Pointer, boost::noncopyable>(
        "ResidualBasedROMBuilderAndSolver", init<typename DenseLinearSolverType::Pointer>());

    class_<ResidualBasedIncrementalROMStaticSchemeType/*, bases<BaseSchemeType>*/, boost::noncopyable>(
        "ResidualBasedIncrementalROMStaticScheme")
        .def("Initialize",
             &ResidualBasedIncrementalROMStaticScheme<DenseSparseSpaceType, LocalSpaceType>::Initialize);



    class_< ConvergenceCriteria< DenseSparseSpaceType, LocalSpaceType >, boost::noncopyable > ("DenseConvergenceCriteria", init<>())
        .def("SetActualizeRHSFlag", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
        .def("GetActualizeRHSflag", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::GetActualizeRHSflag)
        .def("PreCriteria", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::PreCriteria)
        .def("PostCriteria", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::PostCriteria)
        .def("Initialize", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::Initialize)
        .def("InitializeNonLinearIteration", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::InitializeNonLinearIteration)
        .def("InitializeSolutionStep", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::InitializeSolutionStep)
        .def("FinalizeNonLinearIteration", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::FinalizeSolutionStep)
        .def("Check", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::Check)
        .def("SetEchoLevel", &ConvergenceCriteria<DenseSparseSpaceType, LocalSpaceType >::SetEchoLevel)
        ;

    class_< DisplacementCriteria<DenseSparseSpaceType, LocalSpaceType >,
        bases<ConvergenceCriteria< DenseSparseSpaceType, LocalSpaceType > >,
        boost::noncopyable >
            ("DenseDisplacementCriteria", init< double, double>())
                .def("SetEchoLevel",&ResidualCriteria<DenseSparseSpaceType, LocalSpaceType >::SetEchoLevel)
                .def("SetActualizeRHSFlag",&ResidualCriteria<DenseSparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
        ;

    class_<ResidualCriteria<DenseSparseSpaceType, LocalSpaceType >,
        bases<ConvergenceCriteria< DenseSparseSpaceType, LocalSpaceType > >,
        boost::noncopyable >
            ("DenseResidualCriteria", init< double, double>())
                .def("SetEchoLevel",&ResidualCriteria<DenseSparseSpaceType, LocalSpaceType >::SetEchoLevel)
                .def("SetActualizeRHSFlag",&ResidualCriteria<DenseSparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
        ;

  // Residual Based Newton-Raphson Line Search Strategy
  class_< ResidualBasedNewtonRaphsonStrategyType, boost::noncopyable >
      (
          "ResidualBasedNewtonRaphsonROMStrategy", init < ModelPart&,
                                                          BaseSchemeType::Pointer,
                                                          DenseLinearSolverType::Pointer,
                                                          DenseConvergenceCriteriaType::Pointer,
                                                          ResidualBasedROMBuilderAndSolverType::Pointer, int, bool, bool, bool >())
      .def("Solve", &ResidualBasedNewtonRaphsonStrategyType::Solve)
      .def("Initialize", &ResidualBasedNewtonRaphsonStrategyType::Initialize)
      ;
}

} // namespace Python.

} // Namespace Kratos
