#if !defined(KRATOS_SMALL_DISPLACEMENT_ISOTROPIC_DAMAGE_3D_LAW_H_INCLUDED)
#define KRATOS_SMALL_DISPLACEMENT_ISOTROPIC_DAMAGE_3D_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/constitutive_law.h"

namespace Kratos
{
/**
* Base class of constitutive laws.
*/
class KRATOS_API(MULTISCALE_ROM_APPLICATION) LinearIsotropicDamage3DLaw
    : public ConstitutiveLaw

{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LinearIsotropicDamage3DLaw);

    // public:

    /**
    * Constructor.
    */
    LinearIsotropicDamage3DLaw();

    // Copy constructor
    LinearIsotropicDamage3DLaw(const LinearIsotropicDamage3DLaw &);

    /**
    * Destructor.
    */
    virtual ~LinearIsotropicDamage3DLaw(){};

    /**
    * Clone function (has to be implemented by any derived class)
    * @return a pointer to a new instance of this constitutive law
    * NOTE: implementation scheme:
    *      ConstitutiveLaw::Pointer p_clone(new
    * IsotropicDamageHardeningPlaneStrain2DLaw());
    *      return p_clone;
    */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
    * @return the working space dimension of the current constitutive law
    * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
    */
    SizeType WorkingSpaceDimension() override;

    /**
    * returns the size of the strain vector of the current constitutive law
    * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
    */
    SizeType GetStrainSize() override;

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    bool Has(const Variable<Matrix>& rThisVariable) override;


    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains,
    * ...)
    */
    bool Has(const Variable<array_1d<double, 3>>& rThisVariable) override;

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @param rValue output: the value of the specified variable
    */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @return the value of the specified variable
    */
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @return the value of the specified variable
    */
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;


    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @return the value of the specified variable
    */
    array_1d<double, 3>& GetValue(const Variable<array_1d<double, 3>>& rVariable,
                                          array_1d<double, 3>& rValue) override;

    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    void SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    void SetValue(const Variable<Vector>& rVariable,
                          const Vector& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    void SetValue(const Variable<Matrix>& rVariable,
                          const Matrix& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;


    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    void SetValue(const Variable<array_1d<double, 3>>& rVariable,
                          const array_1d<double, 3>& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * Is called to check whether the provided material parameters in the
    * Properties
    * match the requirements of current constitutive model.
    * @param rMaterialProperties the current Properties to be validated against.
    * @return true, if parameters are correct; false, if parameters are
    * insufficient / faulty
    * NOTE: this has to be implemented by each constitutive model. Returns false
    * in base class since
    * no valid implementation is contained here.
    */
    bool ValidateInput(const Properties& rMaterialProperties) override;

    /**
    * returns the expected strain measure of this constitutive law (by default
    * linear strains)
    * @return the expected strain measure
    */
    StrainMeasure GetStrainMeasure() override;

    /**
    * returns the stress measure of this constitutive law (by default 1st
    * Piola-Kirchhoff stress in voigt notation)
    * @return the expected stress measure
    */
    StressMeasure GetStressMeasure() override;

    /**
    * returns whether this constitutive model is formulated in incremental
    * strains/stresses
    * NOTE: by default, all constitutive models should be formulated in total
    * strains
    * @return true, if formulated in incremental strains/stresses, false
    * otherwise
    */
    bool IsIncremental() override;

    /**
    * This is to be called at the very beginning of the calculation
    * (e.g. from InitializeElement) in order to initialize all relevant
    * attributes of the constitutive law
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    */
    void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues) override;

    /**
    * to be called at the beginning of each solution step
    * (e.g. from Element::InitializeSolutionStep)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    * @param the current ProcessInfo instance
    */
    void InitializeSolutionStep(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const Vector& rShapeFunctionsValues,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * to be called at the end of each solution step
    * (e.g. from Element::FinalizeSolutionStep)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    * @param the current ProcessInfo instance
    */
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * to be called at the beginning of each step iteration
    * (e.g. from Element::InitializeNonLinearIteration)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    * @param the current ProcessInfo instance
    */
    void InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues,
                                              const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * to be called at the end of each step iteration
    * (e.g. from Element::FinalizeNonLinearIteration)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    * @param the current ProcessInfo instance
    */
    void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                            const GeometryType& rElementGeometry,
                                            const Vector& rShapeFunctionsValues,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * Computes the material response in terms of 1st Piola-Kirchhoff stresses
    * and constitutive tensor
    * @see Parameters
    */
    void CalculateMaterialResponsePK1(Parameters& rValues) override;

    /**
    * Computes the material response in terms of 2nd Piola-Kirchhoff stresses
    * and constitutive tensor
    * @see Parameters
    */
    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    /**
    * Computes the material response in terms of Kirchhoff stresses and
    * constitutive tensor
    * @see Parameters
    */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
    * Computes the material response in terms of Cauchy stresses and
    * constitutive tensor
    * @see Parameters
    */
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    /**
    * Updates the material response in terms of 1st Piola-Kirchhoff stresses
    * @see Parameters
    */
    void FinalizeMaterialResponsePK1(Parameters& rValues) override;

    /**
    * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
    * @see Parameters
    */
    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    /**
    * Updates the material response in terms of Kirchhoff stresses
    * @see Parameters
    */
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
    * Updates the material response in terms of Cauchy stresses
    * @see Parameters
    */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    /**
    * This can be used in order to reset all internal variables of the
    * constitutive law (e.g. if a model should be reset to its reference state)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current
    * integration point
    * @param the current ProcessInfo instance
    */
    void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues) override;

    /**
    * This function is designed to be called once to check compatibility with
    * element
    * @param rFeatures
    */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;

    /**
    * This function is designed to be called once to perform all the checks
    * needed
    * on the input provided. Checks can be "expensive" as the function is
    * designed
    * to catch user's errors.
    * @param rMaterialProperties
    * @param rElementGeometry
    * @param rCurrentProcessInfo
    * @return
    */
    int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) override;

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Linear Isotropic Damage 3D constitutive law\n";
    };

protected:
  double mInelasticFlag;
  double mStrainEnergy;
    ///@name Protected static Member Variables
    ///@{
    ///@}

    ///@name Protected member Variables
    ///@{
    // bool flag_C = false;
    double r;
    double r_prev;
    double tau_epsilon;
    // boost::numeric::ublas::matrix<double>& constitutiveMatrix();
    ///@}

    ///@name Protected Operators
    ///@{
    double CalculateQ(double r, const Properties& rMaterialProperties);
    virtual void CalculateConstitutiveMatrix(const Properties& props, Matrix& constitutiveMatrix);
    ///@}

    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    // bool m_initialized;
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

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    ///@}

}; /* class LinearIsotropicDamage3DLaw */
} /* namespace Kratos */
#endif /* KRATOS_SMALL_DISPLACEMENT_ISOTROPIC_DAMAGE_3D_LAW_H_INCLUDED defined \
          */
