#if !defined(KRATOS_RVE_IDENTIFIER_LAW_H_INCLUDED)
#define KRATOS_RVE_IDENTIFIER_LAW_H_INCLUDED

#include <vector>
#include "includes/constitutive_law.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{
class KRATOS_API(MULTISCALE_ROM_APPLICATION) RVELaw : public ConstitutiveLaw
{
protected:
public:
    // Type Definitions
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    // Counted pointer of RVELaw
    KRATOS_CLASS_POINTER_DEFINITION(RVELaw);

    // default constructor, takes modelpart and parameters
    RVELaw(ModelPart::Pointer mpModelPart, Kratos::Parameters param);

    // constructor used by Clone(), takes individual data
    RVELaw(ModelPart::Pointer mpModelPart,
           std::vector<Matrix> B_list,
           std::vector<double> IW_list,
           std::vector<ConstitutiveLaw::Pointer> CL_list,
           std::vector<int> prop_id_list);

    // Clone function (has to be implemented by any derived class)
    // @return a pointer to a new instance of this constitutive law
    ConstitutiveLaw::Pointer Clone() const override;

    // Copy constructor.
    RVELaw(const RVELaw& rOther);

    // TODO define copy assignment constructor
    // Copy assignment constructor.
    // RVELaw(const RVELaw& rOther);

    // Destructor
    ~RVELaw() override;

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

    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues,
                              const ProcessInfo& rCurrentProcessInfo) override;

    bool Has(const Variable<Vector>& rThisVariable) override;

    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    std::size_t GetStrainSize() override
    {
        return 6;
    };

    std::size_t WorkingSpaceDimension() override
    {
        return 3;
    };

    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) override;

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Multiscale HPROM constitutive law";
    };

protected:
private:
    ModelPart::Pointer mpRVEModelPart;

    std::vector<Matrix> mB_vec;
    std::vector<double> mIW_vec;
    std::vector<ConstitutiveLaw::Pointer> mCL_vec;
    std::vector<int> mPropId_vec;
    Vector mModesWeights;

    void solve(const Matrix&, const Vector&, Vector&);

    void
    accumulate(Matrix &A, Vector &residual, const Vector &strain_macro);

    void calculate_individual_material_response(Vector &, Matrix &, Vector &, std::size_t);
/*
    int determinant_sign(const permutation_matrix<std::size_t>& pm)
    {
        int pm_sign = 1;
        std::size_t size = pm.size();
        for (std::size_t i = 0; i < size; ++i)
            // swap_rows would swap a pair of rows here, so we change sign
            if (i != pm(i))
                pm_sign *= -1.0;
        return pm_sign;
    };

    double determinant(Matrix m)
    {
        permutation_matrix<std::size_t> pm(m.size1());
        double det = 1.0;
        if (lu_factorize(m, pm))
        {
            det = 0.0;
        }
        else
        {
            for (auto i = 0; i < m.size1(); i++)
                det *= m(i, i); // multiply by elements on diagonal
            det = det * determinant_sign(pm);
        }
        return det;
    };
*/
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

}; // Class RVELaw
} // namespace Kratos.
#endif // KRATOS_RVE_LAW_H_INCLUDED  defined
