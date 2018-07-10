//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Marcelo Raschi
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

/**
 * @class LagrangeMultiplierUtility
 * @ingroup MultiscaleROMApplication
 * @brief Add the necessary DOF to represent the lagrange multipliers
 * @details
 * @author Vicente Mataix Ferrandiz
 * @author Marcelo Raschi
*/
class KRATOS_API(MULTISCALE_ROM_APPLICATION) LagrangeMultiplierUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part
     */
    LagrangeMultiplierUtility(ModelPart& rModelPart)
        : mMainModelPart(rModelPart)
    {
    }

    // Destructor
    ~LagrangeMultiplierUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void Execute()
    {
        NodeType::Pointer pNode = mMainModelPart.pGetNode(1);
        mMainModelPart.GetProcessInfo()[LAGRANGE_MULTIPLIER_NODE] = pNode;

        pNode->AddDof(LAGRANGE_MULTIPLIER_1);
        pNode->AddDof(LAGRANGE_MULTIPLIER_2);
        pNode->AddDof(LAGRANGE_MULTIPLIER_3);
        pNode->AddDof(LAGRANGE_MULTIPLIER_4);
        pNode->AddDof(LAGRANGE_MULTIPLIER_5);
        pNode->AddDof(LAGRANGE_MULTIPLIER_6);
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
