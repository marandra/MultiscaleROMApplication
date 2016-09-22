// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //Application variables creation: (see solid_mechanics_application_variables.cpp)

  //Application Constructor:

  KratosMultiscaleROMApplication::KratosMultiscaleROMApplication():
    mSmallDisplacementElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSmallDisplacementElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mSmallDisplacementElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mSmallDisplacementElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mSmallDisplacementElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mSmallDisplacementElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )

  {}

  void KratosMultiscaleROMApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "Initializing KratosMultiscaleROMApplication...  " << std::endl;

    //Register Variables (variables created in solid_mechanics_application_variables.cpp)

    //explicit schemes
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

    //solution   
    KRATOS_REGISTER_VARIABLE( WRITE_ID )
    KRATOS_REGISTER_VARIABLE( PREVIOUS_DELTA_TIME )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_ALPHA )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_BETA )


    //constitutive law
    KRATOS_REGISTER_VARIABLE( IMPLEX )
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NAME )
       
    //condition nodal load variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_SURFACE_LOAD )

    
    //material orientation
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DX )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DY )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DZ )
    
    //othotropic/anisotropic constants
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )
    
    //material : hyperelastic_plastic
    KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS )
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( HARDENING_EXPONENT )
    KRATOS_REGISTER_VARIABLE( REFERENCE_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( INFINITY_HARDENING_MODULUS )

    //material : isotropic damage
    KRATOS_REGISTER_VARIABLE( DAMAGE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( DAMAGE_THRESHOLD )
    KRATOS_REGISTER_VARIABLE( STRENGTH_RATIO )
    KRATOS_REGISTER_VARIABLE( FRACTURE_ENERGY )

    //thermal
    KRATOS_REGISTER_VARIABLE( THERMAL_EXPANSION_COEFFICIENT );
    KRATOS_REGISTER_VARIABLE( REFERENCE_TEMPERATURE );
    KRATOS_REGISTER_VARIABLE( PLASTIC_DISSIPATION );
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_DISSIPATION );

    //element
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_VECTOR )

    KRATOS_REGISTER_VARIABLE( MATERIAL_STIFFNESS_MATRIX )
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS_MATRIX )

    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS )

    //nodal dofs
    KRATOS_REGISTER_VARIABLE( PRESSURE_REACTION )  
     
    //Register Elements

    //Register small displacement elements
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D3N", mSmallDisplacementElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D4N", mSmallDisplacementElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D6N", mSmallDisplacementElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D8N", mSmallDisplacementElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D9N", mSmallDisplacementElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D4N", mSmallDisplacementElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D6N", mSmallDisplacementElement3D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D8N", mSmallDisplacementElement3D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D10N", mSmallDisplacementElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D15N", mSmallDisplacementElement3D15N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D20N", mSmallDisplacementElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D27N", mSmallDisplacementElement3D27N )

    KRATOS_REGISTER_VARIABLE( INFINITY_YIELD_STRESS )
    //Register Constitutive Laws
    Serializer::Register( "LinearIsotropicDamagePlaneStrain2DLaw", mLinearIsotropicDamagePlaneStrain2DLaw );
   }

}  // namespace Kratos.


