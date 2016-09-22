#include "solid_mechanics_application_variables.h"

namespace Kratos {
  ///@name Type Definitions
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;

  //Create Variables
  KRATOS_CREATE_VARIABLE( double, ISOTROPIC_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, INFINITY_YIELD_STRESS )
}
