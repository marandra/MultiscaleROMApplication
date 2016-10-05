// System includes
// External includes
// Project includes
#include "includes/variables.h"
#include "includes/serializer.h"
#include "multiscale_rom_application.h"
namespace Kratos {
    //Application variables creation: (see solid_mechanics_application_variables.cpp)
    //Application Constructor:
    KratosMultiscaleROMApplication::KratosMultiscaleROMApplication()//:
    {}
    
    void KratosMultiscaleROMApplication::Register() {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosMultiscaleROMApplication...  " << std::endl;
        
        //Register Variables (variables created in solid_mechanics_application_variables.cpp)
        KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
        KRATOS_REGISTER_VARIABLE( INFINITY_YIELD_STRESS )
        
        //Register Constitutive Laws
        Serializer::Register( "LinearIsotropicDamagePlaneStrain2DLaw", mLinearIsotropicDamagePlaneStrain2DLaw );
        Serializer::Register( "ExponentialIsotropicDamagePlaneStrain2DLaw", mExponentialIsotropicDamagePlaneStrain2DLaw );
    }
}

