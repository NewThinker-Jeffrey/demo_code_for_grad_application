#ifndef TrackerModel2_H
#define TrackerModel2_H

#include "Kalman/UKF.h"
#include "MathToolsCpp/Quaternion.h"

IMPORT_JK_Math
IMPORT_JK_KALMAN


struct TrackerModel_CL_UKF_Info {

    int process_status_dim ;
    int process_noise_dim ;
    int measurement_raw_dim ; // if measurement_raw_dim is 0, there is no measurement so measurement_dim and measurement_noise_dim must be 0 too.
    int measurement_dim ;
    int measurement_noise_dim ;
    int control_dim ;
    int constraint_adapt_data_size ; // (in float)
    int constraint_adapt_aux_buf_size_per_item ; // (in float)

    int data_offset ;
    int process_noise_cov_sr_ready ;
    int measurement_noise_cov_sr_ready ;
    int process_status_cov_sr_ready ;

    float kappa ;
    float alpha ;
    float beta  ;

};




class TrackerModel2 : public UKF::Model
{
public:
    struct StatusMap ;
    struct MeasureMap ;
    struct Control ;
    struct Measurement ;
    struct ProcessStatus ;
    struct TangentMeasurement ;
    struct TangentProcessStatus ;
    struct ProcessNoise ;
    struct MeasurementNoise ;
    struct LHInfo ;

    TrackerModel2();

    ~TrackerModel2();

    virtual Matrix f(const Matrix& Xk, const Matrix& Uk=Matrix(),
                     const Matrix& nonAdditive_Vk=Matrix()) ;

    virtual Matrix h(const Matrix& Xpre, const Matrix& nonAdditive_Nk=Matrix()) ;

    const ProcessStatus* getCurStatus() const ;

    void reset() ;

    bool isInitialized() const ;

    void init (const Vector3D& P, const Quaternion& Q, const Matrix& TransMat) ;

    void updateQ(bool LH_present=false, bool updateValue=false) ;

    void updateR(const Vector3D& refZ=Vector3D(0,0,0)) ;

    void setLHInfo (const LHInfo& info) ;

    LHInfo getLHInfo () ;

    void initLHInfo (TrackerModel2::LHInfo& lhInfo, Matrix& Mat) ;
    void CL_UKF_Info(TrackerModel_CL_UKF_Info* info, int offset, bool measure_valid ) ;
    void CL_UKF_Data (TrackerModel_CL_UKF_Info* info, float* data, const Matrix& control, const Matrix& measurement, bool measure_valid) ;
    int CL_UKF_Result (TrackerModel_CL_UKF_Info* info, float* data) ;
    static int CL_UKF_Data_Size (const TrackerModel_CL_UKF_Info* info) ;


public :

    bool importParam (const char* fn) ;
    bool setParam(const Float* p);

    bool exportParam(const char* fn) ;

    bool importParam (const Matrix& mat) ;
    Matrix exportParam() ;

    //// ProcessNoise
    Float IMU_GyroNoise, IMU_AccNoise ;
    Float IMU_GyroOffsetNoise,  IMU_AccOffsetNoise, GravityNoise ;
    Float LH_MatPosNoise, LH_MatOriNoise ;

    //// MeasurementNoise
    Float LH_MeasuredPosNoise,  LH_MeasuredQuatNoise ;

    //// InitDeviation
    Float Position_InitDeviation                 ;
    Float Velocity_InitDeviation                 ;
    Float Quaternion_Tangent_InitDeviation       ;
    Float IMU_GyroOffset_InitDeviation           ;
    Float IMU_AccOffset_InitDeviation            ;
    Float Gravity_InitDeviation                  ;
    Float TransMat_RotationPart_InitDeviation    ;
    Float TransMat_TranslationPart_InitDeviation ;

protected :
    virtual bool setPosterior(const Matrix& X, const Matrix& Pxx) ;
    UKF::NoiseDescrip Q_LHpresent ;
    ProcessStatus* CurStatus ;
    StatusMap* statusMap ;
    MeasureMap* measureMap ;
    bool Initialized ;
};

struct TrackerModel2::Measurement
{
    Vector3D   LH_Pos ;
    Quaternion LH_Ori ;  Vector3D LH_Ori_Tangent ;

    void fromMat(const Matrix& mat) ;
    Matrix toMat() const ;
    void fromTangentMat(const Matrix& mat) ;
    Matrix toTangentMat() const ;
};

struct TrackerModel2::MeasurementNoise
{
    void fromMat(const Matrix& mat) ;
};

struct TrackerModel2::Control
{
    Vector3D IMU_Acc ;
    Vector3D IMU_Gyro ;
    Float dt ;

    void fromMat(const Matrix& mat) ;
    Matrix toMat() const ;
};

struct TrackerModel2::ProcessStatus
{
    Vector3D P, V;
    Quaternion Q ;  Vector3D Q_Tangent ;
    Vector3D IMU_GyroOffset, IMU_AccOffset;
    Float Gravity;
    Quaternion LHQ ; Vector3D LHQ_Tangent ;
    Vector3D LHT ;
    bool LH_present ;

    void fromMat(const Matrix& mat) ;
    Matrix toMat() const ;
    void fromTangentMat(const Matrix& mat) ;
    Matrix toTangentMat() const ;

};

struct TrackerModel2::ProcessNoise
{
    Vector3D IMU_GyroNoise, IMU_AccNoise, IMU_GyroOffsetNoise, IMU_AccOffsetNoise ;
    Float GravityNoise ;
    Vector3D LHQ_Tangent_Noise ;
    Vector3D LHT_Noise ;

    void fromMat(const Matrix& mat) ;
};

struct TrackerModel2::LHInfo
{
    Matrix Value ;
    Matrix Cov ;
    Quaternion LHQ_Ref ;
    Matrix LHQ_Base ;
    Vector3D LHT ;
    Quaternion LHQ ;

    Matrix getTransMat() ;
};



struct TrackerModel2::StatusMap : public JKKalman::ManifoldMap
{
    StatusMap(TrackerModel2* model, const Matrix& refPoint=Matrix()) ;
    ~StatusMap() ;
    virtual Matrix calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const ;
    void setRefPoint(const Matrix& refPointOnManifold) ;
    Matrix project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const ;
    Matrix map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const ;
    Matrix map2Tangent_singleCol(const Matrix& pointsOnManifold) const ;

protected :
    Matrix Q_Base ;
    Quaternion Q_Ref ;
    Matrix LHQ_Base ;
    Quaternion LHQ_Ref ;
    bool LHQ_vaild ;
    TrackerModel2* model ;
    friend class TrackerModel2 ;
};

struct TrackerModel2::MeasureMap : public JKKalman::ManifoldMap
{
    MeasureMap(TrackerModel2* model, const Matrix& refPoint=Matrix()) ;
    ~MeasureMap() ;
    virtual Matrix calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const ;
    void setRefPoint(const Matrix& refPointOnManifold) ;
    Matrix project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const ;
    Matrix map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const ;
    Matrix map2Tangent_singleCol(const Matrix& pointsOnManifold) const ;

protected :
    Matrix LH_Ori_Base ;
    Quaternion LH_Ori_Ref ;
    TrackerModel2* model ;
    friend class TrackerModel2 ;
};





#endif // TrackerModel2_H
