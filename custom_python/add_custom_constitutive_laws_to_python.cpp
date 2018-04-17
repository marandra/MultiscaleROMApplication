// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
//#include "includes/define_python.h":
#include "includes/constitutive_law.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/exponential_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_isotropic_damage_3D_law.hpp"
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
    //class_<LinearIsotropicDamagePlaneStrain2DLaw, typename LinearIsotropicDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw>
    //        (m, "LinearIsotropicDamagePlaneStrain2DLaw").def(init<>() )
    //        ;

    //class_<ExponentialIsotropicDamagePlaneStrain2DLaw, typename ExponentialIsotropicDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw>
    //        (m, "ExponentialIsotropicDamagePlaneStrain2DLaw").def(init<>() )
    //        ;

    //class_<LinearIsotropicDamage3DLaw, typename LinearIsotropicDamage3DLaw::Pointer, ConstitutiveLaw>
    //        (m, "LinearIsotropicDamage3DLaw").def(init<>() )
    //        ;

    class_<RVELaw, typename RVELaw::Pointer, ConstitutiveLaw>
    (m, "RVELaw").def(init<ModelPart::Pointer, Parameters>() )
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
