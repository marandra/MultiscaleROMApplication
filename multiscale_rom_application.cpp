// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "multiscale_rom_application.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/prism_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"

namespace Kratos
{
KratosMultiscaleROMApplication::KratosMultiscaleROMApplication():
     KratosApplication("MultiscaleROMApplication"),
      mSmallDisplacementCustom2D4N(0,
          Element::GeometryType::Pointer(
                  new Quadrilateral2D4<Node<3>>(
                          Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementCustom3D8N(0,
          Element::GeometryType::Pointer(
                  new Hexahedra3D8<Node<3>>(
                          Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacementCustom3D6N(0,
          Element::GeometryType::Pointer(
                  new Prism3D6<Node<3>>(
                          Element::GeometryType::PointsArrayType(6)))),
      //mMinimalKineticCondition2D3N(0,
      //    Condition::GeometryType::Pointer(
      //            new Line2D3<Node<3>>(
      //                    Condition::GeometryType::PointsArrayType(3)))),
      //mMinimalKineticCondition3D4N(0,
      //    Condition::GeometryType::Pointer(
      //            new Quadrilateral3D4<Node<3>>(
      //                    Condition::GeometryType::PointsArrayType(4)))),
      mMinimalKineticCondition3D8N(0,
          Condition::GeometryType::Pointer(
                  new Hexahedra3D8<Node<3>>(
                          Condition::GeometryType::PointsArrayType(8)))),
      mRVELaw()
{
}

void KratosMultiscaleROMApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    KRATOS_INFO("") << "               ___    __   __             ___          __   __   __       " << std::endl
                    << "|\\/| |  | |     |  | /__` /  `  /\\  |    |__     |__| |__) |__) /  \\  |\\/|" << std::endl
                    << "|  | \\__/ |___  |  | .__/ \\__, /~~\\ |___ |___    |  | |    |  \\ \\__/  |  |" << std::endl
                    << "                                                                TECHNOLOGY" << std::endl;

    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_1);
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_2);
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_3);
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_4);
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_5);
    KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_6);
    KRATOS_REGISTER_VARIABLE(REDUCED_MODES_WEIGHTS);
    KRATOS_REGISTER_VARIABLE(MODE_INDEX)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_VECTOR_LIST)
    KRATOS_REGISTER_VARIABLE(GLOBAL_MODES_MATRIX)
    KRATOS_REGISTER_VARIABLE(RHS_MATRIX)
    KRATOS_REGISTER_VARIABLE(LHS_MATRIX)
   // KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_NODE);

    KRATOS_REGISTER_ELEMENT("SmallDisplacementCustomElement2D4N", mSmallDisplacementCustom2D4N);
    KRATOS_REGISTER_ELEMENT("SmallDisplacementCustomElement3D8N", mSmallDisplacementCustom3D8N);
    KRATOS_REGISTER_ELEMENT("SmallDisplacementCustomElement3D6N", mSmallDisplacementCustom3D6N);

    //KRATOS_REGISTER_CONDITION("MinimalKineticCondition2D3N", mMinimalKineticCondition2D3N);
    //KRATOS_REGISTER_CONDITION("MinimalKineticCondition3D4N", mMinimalKineticCondition3D4N);
    KRATOS_REGISTER_CONDITION("MinimalKineticCondition3D8N", mMinimalKineticCondition3D8N);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("RVELaw", mRVELaw);
}
}
