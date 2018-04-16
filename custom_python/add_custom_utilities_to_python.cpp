//  License:		 BSD License
//					 license:
//
//  Main authors:
//

// System includes

// External includes

// Project includes
//#include "includes/define_python.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
// Utilities
#include "custom_utilities/lagrange_multiplier_utility.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<LagrangeMultiplierUtility>(m, "LagrangeMultiplierUtility")
    .def(init<ModelPart&>())
    .def("Execute", &LagrangeMultiplierUtility::Execute)
    ;
}

} // namespace Python.

}  // namespace Kratos

