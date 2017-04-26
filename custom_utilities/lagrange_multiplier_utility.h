//
// Created by mcaicedo on 3/29/17.
//

#ifndef KRATOSMULTIPHYSICS_LAGRANGE_MULTIPLIER_UTILITY_H
#define KRATOSMULTIPHYSICS_LAGRANGE_MULTIPLIER_UTILITY_H

#include "includes/model_part.h"
#include "multiscale_rom_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

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

/** \brief LagrangeMultiplierUtility
 * Creates blah blah blah
 */
class LagrangeMultiplierUtility
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    LagrangeMultiplierUtility(ModelPart& MainModelPart)
        : mMainModelPart(MainModelPart)
    {
        // Hello
    }

    ~LagrangeMultiplierUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void Execute()
    {
        Node<3>::Pointer pNode = mMainModelPart.pGetNode(0);
        mMainModelPart.GetProcessInfo()[LAGRANGE_MULTIPLIER_NODE] = pNode;
    }

    ///@}
    ///@name Operations
    ///@{

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mMainModelPart;

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class LagrangeMultiplierUtility
}

#endif // KRATOSMULTIPHYSICS_LAGRANGE_MULTIPLIER_UTILITY_H
