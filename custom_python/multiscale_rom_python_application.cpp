// System includes
#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>
// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_elements/small_displacement_bbar_element.hpp"
#include "custom_python/add_custom_strategies_to_python.h"
#include "multiscale_rom_application.h"

namespace Kratos {
    namespace Python {
        using namespace boost::python;
        BOOST_PYTHON_MODULE(KratosMultiscaleROMApplication) {
            class_<KratosMultiscaleROMApplication, KratosMultiscaleROMApplication::Pointer,
                   bases<KratosApplication>, boost::noncopyable >("KratosMultiscaleROMApplication");

            AddCustomConstitutiveLawsToPython();
            AddCustomStrategiesToPython();

            //registering variables in python ( if must to be seen from python )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(ISOTROPIC_DAMAGE_MODULUS)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INFINITY_YIELD_STRESS)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(FLOW_RULE_IS_TRACTION_ONLY)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_1)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_2)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_3)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INITIAL_STRAIN_VECTOR)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(NUMBER_REDUCED_MODES)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(REDUCED_MODES_MATRIX)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(REDUCED_MODES_WEIGHTS)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INTEGRATION_POINT_WEIGHT)
        }
    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
