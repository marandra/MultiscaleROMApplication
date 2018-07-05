#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "multiscale_rom_application_variables.h"
#include "multiscale_rom_application.h"
//#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;



PYBIND11_MODULE(KratosMultiscaleROMApplication, m)
{

    class_<KratosMultiscaleROMApplication,
           KratosMultiscaleROMApplication::Pointer,
           KratosApplication >(m, "KratosMultiscaleROMApplication")
            .def(init<>())
            ;

    //AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);

    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ISOTROPIC_DAMAGE_MODULUS)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INFINITY_YIELD_STRESS)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLOW_RULE_IS_TRACTION_ONLY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_4)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_5)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_MULTIPLIER_6)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_REDUCED_MODES)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REDUCED_MODES_MATRIX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REDUCED_MODES_WEIGHTS)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTEGRATION_POINT_WEIGHT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GREEN_LAGRANGE_STRAIN_VECTOR)
}
} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
