#if !defined(KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"

#include "containers/flags.h"

#include "custom_elements/small_displacement_element.hpp"
#include "solid_mechanics_application_variables.h"
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"

namespace Kratos
{
///@name Type Definitions
///@{
typedef array_1d<double,3> Vector3;
typedef array_1d<double,6> Vector6;
///@}

///@name Kratos Globals
///@{

//Application variables definition:  (see solid_mechanics_application_variables.h)

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
 class KratosMultiscaleROMApplication : public KratosApplication
 {
 public:


   ///@name Type Definitions
   ///@{


   /// Pointer definition of KratosSolidMechanicsApplication
   KRATOS_CLASS_POINTER_DEFINITION(KratosMultiscaleROMApplication);


   ///@}
   ///@name Life Cycle
   ///@{

   /// Default constructor.
   KratosMultiscaleROMApplication();

   /// Destructor.
   virtual ~KratosMultiscaleROMApplication() {}


   ///@}
   ///@name Operators
   ///@{


   ///@}
   ///@name Operations
   ///@{

   virtual void Register();



   ///@}
   ///@name Access
   ///@{


   ///@}
   ///@name Inquiry
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   /// Turn back information as a string.
   virtual std::string Info() const
   {
     return "KratosMultiscaleROMApplication";
   }

   /// Print information about this object.
   virtual void PrintInfo(std::ostream& rOStream) const
   {
     rOStream << Info();
     PrintData(rOStream);
   }

   ///// Print object's data.
   virtual void PrintData(std::ostream& rOStream) const
   {
     KRATOS_WATCH( "in KratosMultiscaleROMApplication" )
     KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
     rOStream << "Variables:" << std::endl;
     KratosComponents<VariableData>().PrintData(rOStream);
     rOStream << std::endl;
     rOStream << "Elements:" << std::endl;
     KratosComponents<Element>().PrintData(rOStream);
     rOStream << std::endl;
     rOStream << "Conditions:" << std::endl;
     KratosComponents<Condition>().PrintData(rOStream);
   }


   ///@}
   ///@name Friends
   ///@{


   ///@}

 protected:
   ///@name Protected static Member Variables
   ///@{


   ///@}
   ///@name Protected member Variables
   ///@{


   ///@}
   ///@name Protected Operators
   ///@{


   ///@}
   ///@name Protected Operations
   ///@{


   ///@}
   ///@name Protected  Access
   ///@{


   ///@}
   ///@name Protected Inquiry
   ///@{


   ///@}
   ///@name Protected LifeCycle
   ///@{


   ///@}

 private:

   ///@name Static Member Variables
   ///@{


   ///@}
   ///@name Member Variables
   ///@{


   //small displacement
   const SmallDisplacementElement mSmallDisplacementElement2D3N;
   const SmallDisplacementElement mSmallDisplacementElement2D4N;
   const SmallDisplacementElement mSmallDisplacementElement2D6N;
   const SmallDisplacementElement mSmallDisplacementElement2D8N;
   const SmallDisplacementElement mSmallDisplacementElement2D9N;

   const SmallDisplacementElement mSmallDisplacementElement3D4N;
   const SmallDisplacementElement mSmallDisplacementElement3D6N;
   const SmallDisplacementElement mSmallDisplacementElement3D8N;
   const SmallDisplacementElement mSmallDisplacementElement3D10N;
   const SmallDisplacementElement mSmallDisplacementElement3D15N;
   const SmallDisplacementElement mSmallDisplacementElement3D20N;
   const SmallDisplacementElement mSmallDisplacementElement3D27N;

   //constitutive laws
   const LinearIsotropicDamagePlaneStrain2DLaw   mLinearIsotropicDamagePlaneStrain2DLaw;

   ///@}
   ///@name Private Operators
   ///@{


   ///@}
   ///@name Private Operations
   ///@{


   ///@}
   ///@name Private  Access
   ///@{


   ///@}
   ///@name Private Inquiry
   ///@{


   ///@}
   ///@name Un accessible methods
   ///@{

   /// Assignment operator.
   KratosMultiscaleROMApplication& operator=(KratosMultiscaleROMApplication const& rOther);

   /// Copy constructor.
   KratosMultiscaleROMApplication(KratosMultiscaleROMApplication const& rOther);


   ///@}

 }; // Class KratosSolidMechanicsApplication

 ///@}


 ///@name Type Definitions
 ///@{


 ///@}
 ///@name Input and output
 ///@{

 ///@}


}  // namespace Kratos.

#endif // KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED  defined 


