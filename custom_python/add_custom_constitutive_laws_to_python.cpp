// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/exponential_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/small_displacement_elasto_plastic_J2_3D_law.hpp"
#include "custom_constitutive/rve_law.h"

//For RVELaw
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef Properties::Pointer                    PropertiesPointer;
typedef Mesh<Node<3>, Properties, Element, Condition>   MeshType;
typedef ConstitutiveLaw                  ConstitutiveLawBaseType;
typedef ConstitutiveLaw::Pointer          ConstitutiveLawPointer;
typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;


void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
                                  ConstitutiveLawPointer ThisConstitutiveLaw )
{
    ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
}

void  AddCustomConstitutiveLawsToPython()
{
    class_< LinearIsotropicDamagePlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearIsotropicDamagePlaneStrain2DLaw",
      init<>() )
    ;

    class_< ExponentialIsotropicDamagePlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "ExponentialIsotropicDamagePlaneStrain2DLaw",
      init<>() )
    ;

    class_< LinearElasticPlasticJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlasticJ2PlaneStrain2DLaw",
      init<>() )
    ;

    class_< SmallDisplacementElastoPlasticJ23DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "SmallDisplacementElastoPlasticJ23DLaw",
      init<>() )
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > SolvingStrategyType;
    class_< RVELaw<SolvingStrategyType>, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "RVELaw",
      init<typename SolvingStrategyType::Pointer>() )
    ;

}

}  // namespace Python.
}  // namespace Kratos.
