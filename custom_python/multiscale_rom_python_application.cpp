// System includes
#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>
// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "multiscale_rom_application.h"

namespace Kratos {
    namespace Python {
        using namespace boost::python;
        BOOST_PYTHON_MODULE(KratosMultiscaleROMApplication) {
            class_<KratosMultiscaleROMApplication, KratosMultiscaleROMApplication::Pointer,
                   bases<KratosApplication>, boost::noncopyable >("KratosMultiscaleROMApplication");
        
            AddCustomConstitutiveLawsToPython();
            AddCustomUtilitiesToPython();
            
            //registering variables in python ( if must to be seen from python )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(ISOTROPIC_DAMAGE_MODULUS)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INFINITY_YIELD_STRESS)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(FLOW_RULE_IS_TRACTION_ONLY)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_1)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_2)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_3)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_4)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_5)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_6)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INITIAL_STRAIN_VECTOR)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(INELASTICITY_FLAG)
            //KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_MULTIPLIER_NODE)

        }
    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
