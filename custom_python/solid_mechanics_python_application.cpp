// System includes
#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>
// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "solid_mechanics_application.h"

namespace Kratos {
    namespace Python {
        using namespace boost::python;
        BOOST_PYTHON_MODULE(KratosMultiscaleROMApplication) {
            class_<KratosMultiscaleROMApplication, KratosMultiscaleROMApplication::Pointer,
                   bases<KratosApplication>, boost::noncopyable >("KratosMultiscaleROMApplication");
        
            AddCustomConstitutiveLawsToPython();   
            
            //registering variables in python ( if must to be seen from python )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( INFINITY_YIELD_STRESS )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( RVE_FULL_DISPLACEMENT )
        }
    }  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
