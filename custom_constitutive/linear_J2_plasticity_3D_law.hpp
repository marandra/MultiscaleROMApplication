#if !defined(KRATOS_SMALL_DISPLACEMENTS_ELASTO_PLASTIC_J2_3D_LAW_H_INCLUDED)
#define KRATOS_SMALL_DISPLACEMENTS_ELASTO_PLASTIC_J2_3D_LAW_H_INCLUDED

#include "includes/constitutive_law.h"

namespace Kratos
{
class KRATOS_API(MULTISCALE_ROM_APPLICATION) LinearJ2Plasticity3DLaw
    : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LinearJ2Plasticity3DLaw);
    LinearJ2Plasticity3DLaw();
    virtual ~LinearJ2Plasticity3DLaw(){};
    virtual ConstitutiveLaw::Pointer Clone() const;
    virtual SizeType WorkingSpaceDimension();
    virtual SizeType GetStrainSize();
    virtual bool Has(const Variable<double>& rThisVariable);
    virtual bool Has(const Variable<Vector>& rThisVariable);
    virtual bool Has(const Variable<Matrix>& rThisVariable);
    virtual bool Has(const Variable<array_1d<double, 2>>& rThisVariable);
    virtual bool Has(const Variable<array_1d<double, 3>>& rThisVariable);
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue);
    virtual Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue);
    virtual Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue);
    virtual array_1d<double, 2>& GetValue(const Variable<array_1d<double, 2>>& rVariable,
                                          array_1d<double, 2>& rValue);
    virtual array_1d<double, 3>& GetValue(const Variable<array_1d<double, 3>>& rVariable,
                                          array_1d<double, 3>& rValue);
    virtual void SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<Vector>& rVariable,
                          const Vector& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<Matrix>& rVariable,
                          const Matrix& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<array_1d<double, 2>>& rVariable,
                          const array_1d<double, 2>& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<array_1d<double, 3>>& rVariable,
                          const array_1d<double, 3>& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual bool ValidateInput(const Properties& rMaterialProperties);
    virtual StrainMeasure GetStrainMeasure();
    virtual StressMeasure GetStressMeasure();
    virtual bool IsIncremental();
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);
    virtual void InitializeSolutionStep(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const Vector& rShapeFunctionsValues,
                                        const ProcessInfo& rCurrentProcessInfo);
    virtual void FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo);
    virtual void InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues,
                                              const ProcessInfo& rCurrentProcessInfo);
    virtual void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                            const GeometryType& rElementGeometry,
                                            const Vector& rShapeFunctionsValues,
                                            const ProcessInfo& rCurrentProcessInfo);
    virtual void CalculateMaterialResponsePK1(Parameters& rValues);
    virtual void CalculateMaterialResponsePK2(Parameters& rValues);
    virtual void CalculateMaterialResponseKirchhoff(Parameters& rValues);
    virtual void CalculateMaterialResponseCauchy(Parameters& rValues);
    virtual void FinalizeMaterialResponsePK1(Parameters& rValues);
    virtual void FinalizeMaterialResponsePK2(Parameters& rValues);
    virtual void FinalizeMaterialResponseKirchhoff(Parameters& rValues);
    virtual void FinalizeMaterialResponseCauchy(Parameters& rValues);
    virtual void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues);
    virtual void GetLawFeatures(Features& rFeatures);
    virtual int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo);

protected:
    double mInelasticFlag;
    double mStrainEnergy;
    Vector mPlasticStrain;
    Vector mPlasticStrainOld;
    double mAccumulatedPlasticStrain;
    double mAccumulatedPlasticStrainOld;
    double yieldFunction(const double, const Properties& rMaterialProperties);
    double GetDeltaGamma(double norm_s_trial, const Properties& rMaterialProperties);
    double GetSaturationHardening(const Properties& rMaterialProperties);
    double GetPlasticPotential(const Properties& rMaterialProperties);
    virtual void CalculateTangentTensor(double dgamma,
                                        double norm_s_trial,
                                        const Vector& N_new,
                                        const Properties& props,
                                        Matrix& D);
    virtual void CalculateElasticityTensor(const Properties& props, Matrix& ElasticityTensor);

private:
    friend class Serializer;
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }
    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }
};
}
#endif
