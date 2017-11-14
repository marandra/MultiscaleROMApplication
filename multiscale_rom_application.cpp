#include "multiscale_rom_application.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"
#include "includes/variables.h"
namespace Kratos
{
KratosMultiscaleROMApplication::KratosMultiscaleROMApplication() :
      mSmallDisplacementStrElement2D4N( 0,
                Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
                        Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementStrElement3D8N( 0,
                Element::GeometryType::Pointer(
                        new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mSmallDisplacementStrBbarElement2D4N( 0,
//                Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//                        Element::GeometryType::PointsArrayType(4)))),
//      mSmallDisplacementStrBbarElement3D8N( 0,
//                Element::GeometryType::Pointer(
//                        new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mSmallDisplacementStdElement2D4N( 0,
//                Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//                        Element::GeometryType::PointsArrayType(4)))),
//      mSmallDisplacementStdElement3D8N( 0,
//                Element::GeometryType::Pointer(
//                        new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mSmallDisplacementBbarElement2D4N( 0,
//          Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//              Element::GeometryType::PointsArrayType(4)))),
//      mSmallDisplacementBbarElement3D8N( 0,
//          Element::GeometryType::Pointer(
//              new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mSmallDisplacementHpromElement2D4N( 0,
//          Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//              Element::GeometryType::PointsArrayType(4)))),
//      mSmallDisplacementHpromElement3D8N( 0,
//          Element::GeometryType::Pointer(
//              new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mSmallDisplacementHpromJ2Element2D4N( 0,
//          Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//              Element::GeometryType::PointsArrayType(4)))),
//      mSmallDisplacementHpromJ2Element3D8N( 0,
//          Element::GeometryType::Pointer(
//              new Hexahedra3D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
//      mUpdatedLagrangianFbarElement2D4N( 0,
//          Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//              Element::GeometryType::PointsArrayType(4)))),
//      mTotalLagrangianFbarElement2D4N( 0,
//          Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(
//              Element::GeometryType::PointsArrayType(4)))),
      mMinimalKineticCondition2D3N( 0,
          Condition::GeometryType::Pointer(
              new Line2D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
      mMinimalKineticCondition3D4N( 0,
          Condition::GeometryType::Pointer(new Quadrilateral3D4<Node<3>>(
              Condition::GeometryType::PointsArrayType(4))))
{
}

void KratosMultiscaleROMApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosMultiscaleROMApplication...  " << std::endl;

    KRATOS_REGISTER_VARIABLE(ISOTROPIC_DAMAGE_MODULUS)
    KRATOS_REGISTER_VARIABLE(INFINITY_YIELD_STRESS)
    KRATOS_REGISTER_VARIABLE(FLOW_RULE_IS_TRACTION_ONLY)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_1)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_2)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_3)
    KRATOS_REGISTER_VARIABLE(NUMBER_REDUCED_MODES)
    KRATOS_REGISTER_VARIABLE(REDUCED_MODES_MATRIX)
    KRATOS_REGISTER_VARIABLE(REDUCED_MODES_WEIGHTS)
    KRATOS_REGISTER_VARIABLE(INTEGRATION_POINT_WEIGHT)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_4)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_5)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_6)
    KRATOS_REGISTER_VARIABLE(INELASTIC_FLAG)
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_NODE)
    KRATOS_REGISTER_VARIABLE(GAUSS_WEIGHTS)


    KRATOS_REGISTER_ELEMENT("SmallDisplacementStrElement2D4N", mSmallDisplacementStrElement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementStrElement3D8N", mSmallDisplacementStrElement3D8N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementStrBbarElement2D4N", mSmallDisplacementStrBbarElement2D4N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementStrBbarElement3D8N", mSmallDisplacementStrBbarElement3D8N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementStdElement2D4N", mSmallDisplacementStdElement2D4N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementStdElement3D8N", mSmallDisplacementStdElement3D8N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement2D4N", mSmallDisplacementBbarElement2D4N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement3D8N", mSmallDisplacementBbarElement3D8N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementHpromElement2D4N", mSmallDisplacementHpromElement2D4N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementHpromElement3D8N", mSmallDisplacementHpromElement3D8N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementHpromJ2Element2D4N", mSmallDisplacementHpromJ2Element2D4N)
//    KRATOS_REGISTER_ELEMENT("SmallDisplacementHpromJ2Element3D8N", mSmallDisplacementHpromJ2Element3D8N)
//    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFbarElement2D4N", mUpdatedLagrangianFbarElement2D4N)
//    KRATOS_REGISTER_ELEMENT("TotalLagrangianFbarElement2D4N", mTotalLagrangianFbarElement2D4N)
    KRATOS_REGISTER_CONDITION("MinimalKineticCondition2D3N", mMinimalKineticCondition2D3N)
    KRATOS_REGISTER_CONDITION("MinimalKineticCondition3D4N", mMinimalKineticCondition3D4N)

    // Register Constitutive Laws
    Serializer::Register("LinearIsotropicDamagePlaneStrain2DLaw",
                         mLinearIsotropicDamagePlaneStrain2DLaw);
    Serializer::Register("LinearIsotropicDamage3DLaw",
                         mLinearIsotropicDamage3DLaw);
    Serializer::Register("ExponentialIsotropicDamagePlaneStrain2DLaw",
                         mExponentialIsotropicDamagePlaneStrain2DLaw);
//    Serializer::Register("LinearJ2PlasticityPlaneStrain2DLaw",
//                         mLinearJ2PlasticityPlaneStrain2DLaw);
//    Serializer::Register("LinearJ2Plasticity3DLaw",
//                         mLinearJ2Plasticity3DLaw);
}
}
