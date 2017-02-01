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
        KRATOS_REGISTER_VARIABLE(ISOTROPIC_DAMAGE_MODULUS)
        KRATOS_REGISTER_VARIABLE(INFINITY_YIELD_STRESS)
        KRATOS_REGISTER_VARIABLE(FLOW_RULE_IS_TRACTION_ONLY)
        KRATOS_REGISTER_VARIABLE(REFERENCE_TEMPERATURE)
        KRATOS_REGISTER_VARIABLE(DETERMINANT_F)
        KRATOS_REGISTER_VARIABLE(LAGRANGIAN_DOF_1)
        KRATOS_REGISTER_VARIABLE(LAGRANGIAN_DOF_2)
        KRATOS_REGISTER_VARIABLE(LAGRANGIAN_DOF_3)

        KRATOS_REGISTER_CONDITION( "MinimalKineticCondition2D", mMinimalKineticCondition2D)

        //Register Constitutive Laws
        Serializer::Register("LinearIsotropicDamagePlaneStrain2DLaw", mLinearIsotropicDamagePlaneStrain2DLaw);
        Serializer::Register("ExponentialIsotropicDamagePlaneStrain2DLaw", mExponentialIsotropicDamagePlaneStrain2DLaw);


    }
}

