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

namespace Kratos
{
namespace Python
{
using namespace boost::python;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

void AddCustomStrategiesToPython()
{
    // base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    // typedef SolvingStrategy< SparseSpaceType, LocalSpaceType,
    // LinearSolverType > BaseSolvingStrategyType;
    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;
    // typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType >
    // ConvergenceCriteriaType;

    // custom builder_and_solver types
    typedef ResidualBasedROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedROMBuilderAndSolverType;
    typedef ResidualBasedIncrementalROMStaticScheme<SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalROMStaticSchemeType;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    class_<ResidualBasedROMBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable>(
        "ResidualBasedROMBuilderAndSolver", init<LinearSolverType::Pointer>());

    class_<ResidualBasedIncrementalROMStaticSchemeType, bases<BaseSchemeType>, boost::noncopyable>(
        "ResidualBasedIncrementalROMStaticScheme")
        .def("Initialize",
             &ResidualBasedIncrementalROMStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize);
}

} // namespace Python.

} // Namespace Kratos
