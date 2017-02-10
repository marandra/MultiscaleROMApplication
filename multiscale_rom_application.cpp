// System includes
// External includes
// Project includes
#include "includes/variables.h"
#include "includes/serializer.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "multiscale_rom_application.h"
namespace Kratos {
    //Application variables creation: (see solid_mechanics_application_variables.cpp)
    //Application Constructor:
    KratosMultiscaleROMApplication::KratosMultiscaleROMApplication():
        mSmallDisplacementBbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mUpdatedLagrangianFbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMinimalKineticCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
        mMinimalKineticCondition3D3N( 0, Condition::GeometryType::Pointer( new Line3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) )
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
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_1)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_2)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_3)

        KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D4N", mSmallDisplacementBbarElement2D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianFbarElement2D4N", mUpdatedLagrangianFbarElement2D4N )

        KRATOS_REGISTER_CONDITION( "MinimalKineticCondition2D3N", mMinimalKineticCondition2D3N)
        KRATOS_REGISTER_CONDITION( "MinimalKineticCondition3D3N", mMinimalKineticCondition3D3N)

        //KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
        //KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR)
        //KRATOS_REGISTER_VARIABLE(PRESTRESS)

        //duplicated from SolidMechanics, necessary for SmallDisplacementBbarElement
        KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS)
        KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_TENSOR)
        KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_VECTOR)
        KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_VECTOR)
        KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
        KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

        //Register Constitutive Laws
        Serializer::Register("LinearIsotropicDamagePlaneStrain2DLaw", mLinearIsotropicDamagePlaneStrain2DLaw);
        Serializer::Register("ExponentialIsotropicDamagePlaneStrain2DLaw", mExponentialIsotropicDamagePlaneStrain2DLaw);


    }
}

