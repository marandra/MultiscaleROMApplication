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
#include "custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_isotropic_damage_3D_law.hpp"
//#include "custom_elements/small_displacement_str_element.h"
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
    //const SmallDisplacementStrElement mSmallDisplacementStrElement2D4N;
    //const SmallDisplacementStrElement mSmallDisplacementStrElement3D8N;
    // restrictions
    const MinimalKineticCondition2D mMinimalKineticCondition2D3N;
    const MinimalKineticCondition3D mMinimalKineticCondition3D4N;
    // constitutive laws
    const LinearIsotropicDamagePlaneStrain2DLaw mLinearIsotropicDamagePlaneStrain2DLaw;
    const LinearIsotropicDamage3DLaw mLinearIsotropicDamage3DLaw;
    const ExponentialIsotropicDamagePlaneStrain2DLaw mExponentialIsotropicDamagePlaneStrain2DLaw;
    // others
    KratosMultiscaleROMApplication& operator=(KratosMultiscaleROMApplication const& rOther);
    KratosMultiscaleROMApplication(KratosMultiscaleROMApplication const& rOther);
};
}
#endif
