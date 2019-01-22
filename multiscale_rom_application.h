#if !defined(KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED)
#define KRATOS_MULTISCALE_ROM_APPLICATION_H_INCLUDED

#include <iostream>
#include <string>

#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"

#include "containers/flags.h"

//#include "custom_conditions/minimal_kinetic_2D.hpp"
//#include "custom_conditions/minimal_kinetic_3D.hpp"
#include "custom_conditions/minimal_kinetic_vol_3D.hpp"
#include "custom_constitutive/rve_law.h"
#include "custom_elements/small_displacement_custom.h"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
typedef array_1d<double, 3> Vector3;
class KratosMultiscaleROMApplication : public KratosApplication
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultiscaleROMApplication);
    KratosMultiscaleROMApplication();
    ~KratosMultiscaleROMApplication() override
    {
    }
    void Register() override;
    std::string Info() const override
    {
        return "KratosMultiscaleROMApplication";
    }
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }
    void PrintData(std::ostream& rOStream) const override
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
    const SmallDisplacementCustom mSmallDisplacementCustom2D4N;
    const SmallDisplacementCustom mSmallDisplacementCustom3D8N;

    // restrictions
    //const MinimalKineticCondition2D mMinimalKineticCondition2D3N;
    //const MinimalKineticCondition3D mMinimalKineticCondition3D4N;
    const MinimalKineticVolCondition3D mMinimalKineticCondition3D8N;

    // others
    KratosMultiscaleROMApplication& operator=(KratosMultiscaleROMApplication const& rOther);
    KratosMultiscaleROMApplication(KratosMultiscaleROMApplication const& rOther);
    
    // constitutive laws
    const RVELaw mRVELaw;
};
}
#endif
