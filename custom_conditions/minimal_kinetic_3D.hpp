
#ifndef KRATOSMULTIPHYSICS_MINIMAL_KINETIC_3D_H
#define KRATOSMULTIPHYSICS_MINIMAL_KINETIC_3D_H

#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/element.h"

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
    class KRATOS_API(MULTISCALE_ROM_APPLICATION) MinimalKineticCondition3D
            : public Condition
    {
    public:

        ///@name Type Definitions
        ///@{

        /// Pointer definition of MinimalKineticCondition2D
        KRATOS_CLASS_POINTER_DEFINITION( MinimalKineticCondition3D );

        typedef IndexedObject IndexedObjectType;

        typedef Condition BaseType;

        typedef Node<3> NodeType;

        typedef Properties PropertiesType;

        typedef Geometry<NodeType> GeometryType;

        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        typedef Vector VectorType;

        typedef Matrix MatrixType;

        //typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;


        ///@}
        ///@name Life Cycle
        ///@{

        /// Empty constructor needed for serialization
        MinimalKineticCondition3D();

        /// Default constructors
        MinimalKineticCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

        MinimalKineticCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /// Copy constructor.
        MinimalKineticCondition3D(MinimalKineticCondition3D const& rOther);


        /// Destructor.
        virtual ~MinimalKineticCondition3D();


        /// Constructor.
        /** @param NewId Index number of the new condition (optional)
         */
        //MinimalKineticCondition3D(IndexType NewId = 0);

        /// Constructor using an array of nodes
        /**
         @param NewId Index of the new condition
         @param ThisNodes An array containing the nodes of the new condition
         */
        //MinimalKineticCondition3D(IndexType NewId,
        //                          const NodesArrayType& ThisNodes);

        /// Constructor using Geometry
        /**
         @param NewId Index of the new condition
         @param pGeometry Pointer to a geometry object
         */
        //MinimalKineticCondition3D(IndexType NewId,
        //                          GeometryType::Pointer pGeometry);

        /// Constructor using Properties
        /**
         @param NewId Index of the new element
         @param pGeometry Pointer to a geometry object
         @param pProperties Pointer to the element's properties
         */
        //MinimalKineticCondition3D(IndexType NewId,
        //                          GeometryType::Pointer pGeometry,
        //                          PropertiesType::Pointer pProperties);




        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.
        MinimalKineticCondition3D & operator=(MinimalKineticCondition3D const& rOther);

        ///@}
        ///@name Operations
        ///@{

        /// Create a new MinimalKineticCondition2D instance
        Condition::Pointer Create(IndexType NewId,
                                  NodesArrayType const& ThisNodes,
                                  PropertiesType::Pointer pProperties) const;

        /// Check input to ensure that it makes sense.
        int Check(const ProcessInfo& rCurrentProcessInfo);

        /// Returns a matrix of penalty terms for the periodic variables.
        /**
         * The weight of the penalty terms is given by the member variable mWeight,
         * set using SetValueOnIntegrationPoints. The periodic variables are read from
         * the value of PERIODIC_VARIABLES stored in rCurrentProcessInfo.
         * @param rLeftHandSideMatrix Local left hand side matrix (output)
         * @param rRightHandSideVector Local right hand side vector (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);

        /// Returns a matrix of the integral of Shape Functions to compute the penalty matrix in 3D cases.
        /**
         * The weight of the penalty terms is given by the member variable mWeight,
         * set using SetValueOnIntegrationPoints. The periodic variables are read from
         * the value of PERIODIC_VARIABLES stored in rCurrentProcessInfo.
         * @param rNintMatrix matrix of integral of shape functions (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateIntegralOfShapeFunctions(MatrixType& rNintMatrix,
                                                       ProcessInfo& rCurrentProcessInfo);

        /// Returns a matrix of penalty terms for the periodic variables.
        /**
         * @param rLeftHandSideMatrix Local left hand side matrix (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo);

        /// Returns RHS values for the penalized dofs.
        /**
         * @param rRightHandSideVector Local right hand side vector (output)
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo);

        /// Provides the global indices for each one of this element's local rows
        /**
         * this determines the elemental equation ID vector for all elemental
         * DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo ProcessInfo instance (unused)
         */
        virtual void GetDofList(DofsVectorType& ElementalDofList,
                                ProcessInfo& CurrentProcessInfo);

        /// Returns the values of the unknowns for each node
        virtual void GetValuesVector(Vector& Values, int Step = 0);

        /// Returns the cross product of two 3D vectors
        ////virtual void CrossProduct(array_1d<double, 3>& cross, const array_1d<double, 3>& a, const array_1d<double, 3>& b);
        //virtual void CrossProduct(Vector& cross,
        //                          const Vector& a,
        //                          const Vector& b);

        //************* GETTING METHODS

        /**
         * Returns the currently selected integration method
         * @return current integration method selected
         */
        IntegrationMethod GetIntegrationMethod() const;




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
            buffer << "MinimalKineticCondition3D #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "MinimalKineticCondition3D #" << Id();
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


        /**
         * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
         */

        struct GeneralVariables
        {
        private:

            //variables including all integration points
            const GeometryType::ShapeFunctionsGradientsType* pDN_De;
            const Matrix* pNcontainer;

        public:

            //StressMeasureType StressMeasure;

            //for axisymmetric use only
            //double  Radius;

            //general variables for large displacement use
            //double  detF;
            //double  detF0;
            double  detJ;
            //Vector  StrainVector;
            //Vector  StressVector;
            Vector  N;
            //Matrix  B;
            //Matrix  H;
            //Matrix  F;
            //Matrix  F0;
            //Matrix  DN_DX;
            //Matrix  ConstitutiveMatrix;

            //variables including all integration points
            GeometryType::JacobiansType J;
            GeometryType::JacobiansType j;
            Matrix  DeltaPosition;


            /**
             * sets the value of a specified pointer variable
             */
            //void SetShapeFunctionsGradients(const GeometryType::ShapeFunctionsGradientsType &rDN_De)
            //{
            //    pDN_De=&rDN_De;
            //};

            void SetShapeFunctions(const Matrix& rNcontainer)
            {
                pNcontainer=&rNcontainer;
            };

            /**
             * returns the value of a specified pointer variable
             */
            //const GeometryType::ShapeFunctionsGradientsType& GetShapeFunctionsGradients()
            //{
            //    return *pDN_De;
            //};

            const Matrix& GetShapeFunctions()
            {
                return *pNcontainer;
            };

            void Initialize( const unsigned int& voigt_size, const unsigned int& dimension,
                             const unsigned int& number_of_nodes )
            {

                //StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

                //doubles
                //radius
                //Radius = 0;

                //jacobians
                //detF  = 1;
                //detF0 = 1;
                detJ  = 1;

                //vectors
                //StrainVector.resize(voigt_size,false);
                //StressVector.resize(voigt_size,false);
                N.resize(number_of_nodes,false);
                //noalias(StrainVector) = ZeroVector(voigt_size);
                //noalias(StressVector) = ZeroVector(voigt_size);
                noalias(N) = ZeroVector(number_of_nodes);

                //matrices
                //B.resize(voigt_size, dimension*number_of_nodes,false);
                //H.resize(dimension, dimension,false);
                //F.resize(dimension, dimension,false);
                //F0.resize(dimension, dimension,false);
                //DN_DX.resize(number_of_nodes, dimension,false);
                //ConstitutiveMatrix.resize(voigt_size, voigt_size,false);
                DeltaPosition.resize(number_of_nodes, dimension,false);

                //noalias(B)  = ZeroMatrix(voigt_size, dimension*number_of_nodes);
                //noalias(H)  = ZeroMatrix(dimension, dimension);
                //noalias(F)  = IdentityMatrix(dimension);
                //noalias(F0) = IdentityMatrix(dimension);
                //noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);
                //noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
                noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);

                //others
                J.resize(1,false);
                j.resize(1,false);
                J[0].resize(dimension,dimension,false);
                j[0].resize(dimension,dimension,false);
                noalias(J[0]) = ZeroMatrix(dimension,dimension);
                noalias(j[0]) = ZeroMatrix(dimension,dimension);

                //pointers
                pDN_De = NULL;
                pNcontainer = NULL;
            }

        };

        /**
         * Currently selected integration methods
         */
        IntegrationMethod mThisIntegrationMethod;

        /**
         * Initialize Element General Variables
         */
        virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

        /**
         * Calculation of the Integration Weight
         */
        virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);

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
//template class KRATOS_API(MULTISCALE_ROM_APPLICATION) KratosComponents<MinimalKineticCondition2D >;

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream & operator >>(std::istream& rIStream,
                                      MinimalKineticCondition3D& rThis)
    {
        return rIStream;
    }

/// output stream function
    inline std::ostream & operator <<(std::ostream& rOStream,
                                      const MinimalKineticCondition3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << " : " << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
///@}

///@}

} // namespace Kratos.

#endif //KRATOSMULTIPHYSICS_MINIMAL_KINETIC_3D_H
