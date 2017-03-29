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

#include "multiscale_rom_application_variables.h"
#include "custom_conditions/minimal_kinetic_2D.hpp"
#include "custom_conditions/minimal_kinetic_3D.hpp"
#include "custom_elements/small_displacement_bbar_element.hpp"
#include "custom_elements/updated_lagrangian_fbar_element.hpp"
#include "custom_elements/total_lagrangian_fbar_element.hpp"
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/small_displacement_isotropic_damage_3D_law.hpp"
#include "custom_constitutive/exponential_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/small_displacement_elasto_plastic_J2_3D_law.hpp"

namespace Kratos {
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
    class KratosMultiscaleROMApplication : public KratosApplication {
       public:
           /// Pointer definition of KratosSolidMechanicsApplication
           KRATOS_CLASS_POINTER_DEFINITION(KratosMultiscaleROMApplication);
           
           /// Default constructor.
           KratosMultiscaleROMApplication();
           
           /// Destructor.
           virtual ~KratosMultiscaleROMApplication() {}
           
           virtual void Register();
           
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
           
       protected:
       
       private:
           //elements
           const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D4N;
           const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D8N;
           const UpdatedLagrangianFbarElement mUpdatedLagrangianFbarElement2D4N;
           const TotalLagrangianFbarElement mTotalLagrangianFbarElement2D4N;
           //conditions
           const MinimalKineticCondition2D mMinimalKineticCondition2D3N;
           const MinimalKineticCondition3D mMinimalKineticCondition3D4N;
           //constitutive laws
           const LinearIsotropicDamagePlaneStrain2DLaw mLinearIsotropicDamagePlaneStrain2DLaw;
           const SmallDisplacementIsotropicDamage3DLaw mSmallDisplacementIsotropicDamage3DLaw;
           const ExponentialIsotropicDamagePlaneStrain2DLaw mExponentialIsotropicDamagePlaneStrain2DLaw;
           const LinearElasticPlasticJ2PlaneStrain2DLaw mLinearElasticPlasticJ2PlaneStrain2DLaw;
           const SmallDisplacementElastoPlasticJ23DLaw mSmallDisplacementElastoPlasticJ23DLaw;
           /// Assignment operator.
           KratosMultiscaleROMApplication& operator=(KratosMultiscaleROMApplication const& rOther);
           /// Copy constructor.
           KratosMultiscaleROMApplication(KratosMultiscaleROMApplication const& rOther);
    };
}

#endif // KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED  defined 


