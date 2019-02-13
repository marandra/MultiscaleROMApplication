#if !defined(KRATOS_RVE_LAW_H_INCLUDED)
#define KRATOS_RVE_LAW_H_INCLUDED

// System includes
#include <vector>
#include <unordered_map>

// Project includes
#include "includes/constitutive_law.h"

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
class KRATOS_API(MULTISCALE_ROM_APPLICATION) RVELaw : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    typedef std::unordered_map<std::size_t, Properties> PropertiesMap;

    // Counted pointer of RVELaw
    KRATOS_CLASS_POINTER_DEFINITION(RVELaw);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    RVELaw();

    /**
     * @brief Constructor used by Create()
     */
    explicit RVELaw(Kratos::Parameters Params);

    /**
     * @brief Constructor used by Clone()
     */
    RVELaw(PropertiesMap pProperties_list,
           std::vector<Matrix> B_list,
           std::vector<double> IW_list,
           std::vector<ConstitutiveLaw::Pointer> CL_list,
           std::vector<int> prop_id_list,
           double abs_tol, double rel_tol, int max_iter, int verbose);

    /**
     * @brief Copy constructor
     */
    RVELaw(const RVELaw& rOther);

    /**
     * @brief Destructor
     */
    ~RVELaw() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters) const override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::size_t GetStrainSize() override
    {
        return 6;
    };

    std::size_t WorkingSpaceDimension() override
    {
        return 3;
    };

    bool Has(const Variable<Vector>& rThisVariable) override;

    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    void SetValue(
            const Variable<Vector>& rThisVariable,
            const Vector& rValue,
            const ProcessInfo& rCurrentProcessInfo) override;

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
     * @brief To be called at the end of each solution step  (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param rCurrentProcessInfo the current ProcessInfo instance
     */
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues,
                              const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    /**
     * @brief This function provides the place to perform checks on the
     * completeness of the input.
     * @details It is designed to be called only once (or anyway, not often)
     * typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) override;

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Multiscale HPROM constitutive law";
    };

    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

protected:

private:
    PropertiesMap mProperties_map;
    std::vector<Matrix> mB_vec;
    std::vector<double> mIW_vec;
    std::vector<ConstitutiveLaw::Pointer> mCL_vec;
    std::vector<int> mPropId_vec;
    Vector mModesWeights;
    double mAbsoluteTolerance;
    double mRelativeTolerance;
    std::size_t mMaxIteration;
    int mVerbose;

    void Solve(const Matrix &A, const Vector &res, Vector &Dx);

    void Accumulate(Matrix &A, Vector &residual, const Vector &strain_macro, const ProcessInfo &process_info);

    std::string ReadFile(const std::string &filename) const;

    void CalculateIndividualMaterialResponse(Vector &stress, Matrix &constit, Vector &strain,
                                             const ProcessInfo &process_info,
                                             std::size_t ip_index);

    void GetPropertyBlock(Kratos::Parameters Materials);

    void AssignPropertyBlock(Kratos::Parameters Data);

    void TrimComponentName(std::string& rLine);

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class RVELaw
} // namespace Kratos.
#endif // KRATOS_RVE_LAW_H_INCLUDED  defined
