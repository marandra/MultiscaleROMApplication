#include "multiscale_rom_application_variables.h"

namespace Kratos {
  ///@name Type Definitions
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;

  //Create Variables
  KRATOS_CREATE_VARIABLE( double, ISOTROPIC_DAMAGE_MODULUS )
  KRATOS_CREATE_VARIABLE( double, INFINITY_YIELD_STRESS )
  KRATOS_CREATE_VARIABLE( int, FLOW_RULE_IS_TRACTION_ONLY )
  KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_1)
  KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_2)
  KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_3)
  KRATOS_CREATE_VARIABLE(Vector, INITIAL_STRAIN_VECTOR)

}
