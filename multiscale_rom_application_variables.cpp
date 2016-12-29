#include "multiscale_rom_application_variables.h"

namespace Kratos {
  ///@name Type Definitions
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;

  //Create Variables
  KRATOS_CREATE_VARIABLE( double, ISOTROPIC_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, INFINITY_YIELD_STRESS )
  KRATOS_CREATE_VARIABLE( int, FLOW_RULE_IS_TRACTION_ONLY )
//  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( RVE_FULL_DISPLACEMENT )
  KRATOS_CREATE_VARIABLE( double, REFERENCE_TEMPERATURE )
  KRATOS_CREATE_VARIABLE( double, DETERMINANT_F )
}
