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
#include "includes/define.h"
#include "containers/flags.h"
#include "add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_strategies/builders_and_solvers/rom_builder_and_solver.hpp"

namespace Kratos
{

namespace Python
{
using namespace boost::python;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

void  AddCustomStrategiesToPython()
{

    //base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    //typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    //typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

    //custom builder_and_solver types
    typedef ROMBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ROMBuilderAndSolverType;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************


    // Component Wise Builder and Solver
    class_< ROMBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
            (
              "ROMBuilderAndSolver", init< LinearSolverType::Pointer > ()
            );

}

}  // namespace Python.

} // Namespace Kratos

