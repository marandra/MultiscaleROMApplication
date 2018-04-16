//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//


#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"


namespace Kratos
{
namespace Python
{

  void  AddCustomStrategiesToPython(pybind11::module& m);

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
