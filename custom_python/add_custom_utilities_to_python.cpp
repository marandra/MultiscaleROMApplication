//  License:		 BSD License
//					 license:
//
//  Main authors:
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/lagrange_multiplier_utility.h"

namespace Kratos
{
    namespace Python
    {
        void  AddCustomUtilitiesToPython()
        {
            using namespace boost::python;

            class_<LagrangeMultiplierUtility>("LagrangeMultiplierUtility", init<ModelPart&>())
                    .def("Execute",&LagrangeMultiplierUtility::Execute)
                    ;
        }

    }  // namespace Python.

} // Namespace Kratos