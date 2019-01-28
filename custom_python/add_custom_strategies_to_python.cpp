//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "containers/flags.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "custom_strategies/builders_and_solvers/residualbased_block_builder_and_solver_custom.hpp"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    // base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    //typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;
    //typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // custom builder_and_solver types
    typedef ResidualBasedBlockBuilderAndSolverCustom<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverCustomType;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

   pybind11::class_<ResidualBasedBlockBuilderAndSolverCustomType,
        typename ResidualBasedBlockBuilderAndSolverCustomType::Pointer,
       //BuilderAndSolverType> (m, "ResidualBasedBlockBuilderAndSolverCustom")
        ResidualBasedBlockBuilderAndSolverType> (m, "ResidualBasedBlockBuilderAndSolverCustom")
        .def(pybind11::init<LinearSolverType::Pointer>());
}

} // namespace Python.

} // Namespace Kratos
