#include "includes/variables.h"
#include "includes/serializer.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "multiscale_rom_application.h"
namespace Kratos {
    //Application variables creation: (see solid_mechanics_application_variables.cpp)
    //Application Constructor:
    KratosMultiscaleROMApplication::KratosMultiscaleROMApplication():
        mSmallDisplacementBbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mSmallDisplacementBbarElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
        mUpdatedLagrangianFbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mTotalLagrangianFbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMinimalKineticCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
        mMinimalKineticCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ))
    {}

    void KratosMultiscaleROMApplication::Register() {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosMultiscaleROMApplication...  " << std::endl;
        
        //Register Variables (variables created in solid_mechanics_application_variables.cpp)
        KRATOS_REGISTER_VARIABLE(ISOTROPIC_DAMAGE_MODULUS)
        KRATOS_REGISTER_VARIABLE(INFINITY_YIELD_STRESS)
        KRATOS_REGISTER_VARIABLE(FLOW_RULE_IS_TRACTION_ONLY)
        // For 2D cases
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_1)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_2)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_3)
        // For 3D cases
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_4)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_5)
        KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_6)
        // Multiscale macro
        KRATOS_REGISTER_VARIABLE(INITIAL_STRAIN_VECTOR)

        //KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_NODE)

        KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D4N", mSmallDisplacementBbarElement2D4N )
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D8N", mSmallDisplacementBbarElement3D8N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianFbarElement2D4N", mUpdatedLagrangianFbarElement2D4N )
        KRATOS_REGISTER_ELEMENT( "TotalLagrangianFbarElement2D4N", mTotalLagrangianFbarElement2D4N )

        KRATOS_REGISTER_CONDITION( "MinimalKineticCondition2D3N", mMinimalKineticCondition2D3N)
        KRATOS_REGISTER_CONDITION( "MinimalKineticCondition3D4N", mMinimalKineticCondition3D4N)

        //Register Constitutive Laws
        Serializer::Register("LinearIsotropicDamagePlaneStrain2DLaw", mLinearIsotropicDamagePlaneStrain2DLaw);
        Serializer::Register("SmallDisplacementIsotropicDamage3DLaw", mSmallDisplacementIsotropicDamage3DLaw);
        Serializer::Register("ExponentialIsotropicDamagePlaneStrain2DLaw", mExponentialIsotropicDamagePlaneStrain2DLaw);
        Serializer::Register("LinearElasticPlasticJ2PlaneStrain2DLaw", mLinearElasticPlasticJ2PlaneStrain2DLaw);
        Serializer::Register("SmallDisplacementElastoPlasticJ23DLaw", mSmallDisplacementElastoPlasticJ23DLaw);

    }
}

