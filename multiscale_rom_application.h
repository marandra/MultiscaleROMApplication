#if !defined(KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED)
#define KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED

#include <iostream>
#include <string>

#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"

#include "containers/flags.h"

#include "custom_conditions/minimal_kinetic_2D.hpp"
#include "custom_conditions/minimal_kinetic_3D.hpp"
#include "custom_constitutive/exponential_isotropic_damage_plane_strain_2D_law.hpp"
//#include "custom_constitutive/linear_J2_plasticity_plane_strain_2D_law.hpp"
//#include "custom_constitutive/linear_J2_plasticity_3D_law.hpp"
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_isotropic_damage_3D_law.hpp"
#include "custom_elements/small_displacement_str_element.h"
//#include "custom_elements/small_displacement_str_bbar_element.h"
//#include "custom_elements/small_displacement_std_element.hpp"
//#include "custom_elements/small_displacement_bbar_element.hpp"
//#include "custom_elements/small_displacement_hprom_element.hpp"
//#include "custom_elements/small_displacement_hprom_J2_element.hpp"
//#include "custom_elements/total_lagrangian_fbar_element.hpp"
//#include "custom_elements/updated_lagrangian_fbar_element.hpp"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
typedef array_1d<double, 3> Vector3;
class KratosMultiscaleROMApplication : public KratosApplication
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultiscaleROMApplication);
    KratosMultiscaleROMApplication();
    virtual ~KratosMultiscaleROMApplication()
    {
    }
    virtual void Register();
    virtual std::string Info() const
    {
        return "KratosMultiscaleROMApplication";
    }
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in KratosMultiscaleROMApplication")
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size())
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
    // elements
    const SmallDisplacementStrElement mSmallDisplacementStrElement2D4N;
    const SmallDisplacementStrElement mSmallDisplacementStrElement3D8N;
//    const SmallDisplacementStrBbarElement mSmallDisplacementStrBbarElement2D4N;
//    const SmallDisplacementStrBbarElement mSmallDisplacementStrBbarElement3D8N;
//    const SmallDisplacementStdElement mSmallDisplacementStdElement2D4N;
//    const SmallDisplacementStdElement mSmallDisplacementStdElement3D8N;
//    const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D4N;
//    const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D8N;
//    const SmallDisplacementHpromElement mSmallDisplacementHpromElement2D4N;
//    const SmallDisplacementHpromElement mSmallDisplacementHpromElement3D8N;
//    const SmallDisplacementHpromJ2Element mSmallDisplacementHpromJ2Element2D4N;
//    const SmallDisplacementHpromJ2Element mSmallDisplacementHpromJ2Element3D8N;
//    const UpdatedLagrangianFbarElement mUpdatedLagrangianFbarElement2D4N;
//    const TotalLagrangianFbarElement mTotalLagrangianFbarElement2D4N;
    // restrictions
    const MinimalKineticCondition2D mMinimalKineticCondition2D3N;
    const MinimalKineticCondition3D mMinimalKineticCondition3D4N;
    // constitutive laws
    const LinearIsotropicDamagePlaneStrain2DLaw mLinearIsotropicDamagePlaneStrain2DLaw;
    const LinearIsotropicDamage3DLaw mLinearIsotropicDamage3DLaw;
    const ExponentialIsotropicDamagePlaneStrain2DLaw mExponentialIsotropicDamagePlaneStrain2DLaw;
//    const LinearJ2PlasticityPlaneStrain2DLaw mLinearJ2PlasticityPlaneStrain2DLaw;
//    const LinearJ2Plasticity3DLaw mLinearJ2Plasticity3DLaw;
    // others
    KratosMultiscaleROMApplication& operator=(KratosMultiscaleROMApplication const& rOther);
    KratosMultiscaleROMApplication(KratosMultiscaleROMApplication const& rOther);
};
}
#endif
