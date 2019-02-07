//  License:		 BSD License
//					 license:
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/lagrange_multiplier_utility.h"
#include "custom_utilities/modes_to_nodes_matrix_utility.hpp"

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

    class_<ModesToNodesMatrixUtility>(m, "ModeToNodeMatrixUtility")
    .def(init<ModelPart&>())
    .def("Execute", &ModesToNodesMatrixUtility::Execute)
    ;
}

} // namespace Python.

}  // namespace Kratos

