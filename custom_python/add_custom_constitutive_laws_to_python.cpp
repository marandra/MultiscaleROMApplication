// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
//#include "includes/define_python.h":
#include "includes/constitutive_law.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/rve_law.h"

// For RVELaw
//#include "linear_solvers/linear_solver.h"
//#include "solving_strategies/strategies/solving_strategy.h"
//#include "spaces/ublas_space.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;


void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    class_<RVELaw, typename RVELaw::Pointer, ConstitutiveLaw>
    (m, "RVELaw").def(init<>() )
    ;

   // typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
   // typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
   // typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
   // typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
   // class_<RVELaw<SolvingStrategyType>, bases<ConstitutiveLawBaseType>, boost::noncopyable>(
   //     "RVELaw", init<typename SolvingStrategyType::Pointer>());
}

} // namespace Python.
} // namespace Kratos.
