#if !defined(KRATOS_MINIMAL_KINETIC_CONDITION_2D_H_INCLUDED)
#define KRATOS_MINIMAL_KINETIC_CONDITION_2D_H_INCLUDED

#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

#include "geometries/geometry.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/serializer.h"
#include "utilities/indexed_object.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// Condition used to assign periodic boundary conditions.

/**
 * Implements a Force Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(MULTISCALE_ROM_APPLICATION) MinimalKineticCondition2D : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MinimalKineticCondition2D
    KRATOS_CLASS_POINTER_DEFINITION(MinimalKineticCondition2D);

    typedef IndexedObject IndexedObjectType;

    typedef Condition BaseType;

    typedef Node<3> NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /** @param NewId Index number of the new condition (optional)
     */
    MinimalKineticCondition2D(IndexType NewId = 0);

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    MinimalKineticCondition2D(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    MinimalKineticCondition2D(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    MinimalKineticCondition2D(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties);

    /// Copy constructor.
    MinimalKineticCondition2D(MinimalKineticCondition2D const& rOther);

    /// Destructor.
    ~MinimalKineticCondition2D() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MinimalKineticCondition2D& operator=(MinimalKineticCondition2D const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new MinimalKineticCondition2D instance
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    /// Check input to ensure that it makes sense.
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a matrix of penalty terms for the periodic variables.
    /**
     * The weight of the penalty terms is given by the member variable mWeight,
     * set using SetValueOnIntegrationPoints. The periodic variables are read
     * from
     * the value of PERIODIC_VARIABLES stored in rCurrentProcessInfo.
     * @param rLeftHandSideMatrix Local left hand side matrix (output)
     * @param rRightHandSideVector Local right hand side vector (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a matrix of penalty terms for the periodic variables.
    /**
     * @param rLeftHandSideMatrix Local left hand side matrix (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

    /// Returns RHS values for the penalized dofs.
    /**
     * @param rRightHandSideVector Local right hand side vector (output)
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo ProcessInfo instance (unused)
     */
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    /// Returns the values of the unknowns for each node
    void GetValuesVector(Vector& Values, int Step = 0) override;

    ///@}
    ///@name Conditional Data
    ///@{

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
        std::stringstream buffer;
        buffer << "MinimalKineticCondition2D #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MinimalKineticCondition2D #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        Condition::PrintData(rOStream);
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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

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

    ///@}

}; // Class MinimalKineticCondition2D

///@}

// TODO check the function of the template below,
// and if it is OK to comment it out.
// template class KRATOS_API(MULTISCALE_ROM_APPLICATION)
// KratosComponents<MinimalKineticCondition2D >;

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, MinimalKineticCondition2D& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const MinimalKineticCondition2D& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif /* KRATOS_FLUID_PERTIODIC_CONDITION_H */
