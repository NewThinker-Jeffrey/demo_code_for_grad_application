#include <cmath>
#include "Tracker/TrackerModel2.h"

IMPORT_JK_Math
IMPORT_JK_KALMAN

///////////////////////////////////////////////////////////////////////
///
/// Static Function Area
///
///////////////////////////////////////////////////////////////////////


static Matrix Quat2Mat (const Quaternion& q)
{
    Matrix M(4,1) ;
    M(0,0) = q.w() ;
    M(1,0) = q.x() ;
    M(2,0) = q.y() ;
    M(3,0) = q.z() ;
    return M ;
}

static Quaternion Mat2Quat (const Matrix& M)
{
    Quaternion q;
    q.w() = M(0,0) ;
    q.x() = M(1,0) ;
    q.y() = M(2,0) ;
    q.z() = M(3,0) ;
    return q ;
}


static Matrix getOthoCompBase (const Matrix& _V, Float eps=0.0001)
{
    Matrix V = _V.clone() ;
    int n = V.rows(), m = V.cols() ;
    Matrix Base(n,n) ;

    int rankV = 0, rankBase = 0 ;
    int i, j;

    for (i=0; i<m; i++)
    {
        Matrix Vi = V.col(i).clone() ;
        for (j=0; j<rankBase; j++)
        {
            Float dot = Vi.transpose() * Base.col(j) ;
            Vi = Vi - dot * Base.col(j) ;
        }

        Float len = sqrt(Vi.squareSum()) ;
        if (len > eps)
        {
            Base.col(rankBase) = Vi / len ;
            rankBase ++ ;
            rankV ++ ;
            if (rankBase == n)
                break ;
        }
    }

    for (i=0; i<n; i++)
    {
        Matrix D(n,1) ;
        D.fill(0);
        D(i,0) = 1 ;

        for (j=0; j<rankBase; j++)
        {
            Float dot = D.transpose() * Base.col(j) ;
            D = D - dot * Base.col(j) ;
        }

        Float len = sqrt(D.squareSum()) ;
        if (len > eps)
        {
            Base.col(rankBase) = D / len ;
            rankBase ++ ;
            if (rankBase == n)
                break ;
        }
    }

    return Base.subMat(0,rankV, n, n-rankV).clone() ;
}

static Matrix Mat3x4_toSingleCol(const Matrix& in)
{
    Matrix out(12,1) ;
    out.subMat(0,0,3,1) = in.subMat(0,0,3,1) ;
    out.subMat(3,0,3,1) = in.subMat(0,1,3,1) ;
    out.subMat(6,0,3,1) = in.subMat(0,2,3,1) ;
    out.subMat(9,0,3,1) = in.subMat(0,3,3,1) ;
    return out ;
}

static Matrix SingleCol_toMat3x4(const Matrix& in)
{
    Matrix out(3,4) ;
    out.subMat(0,0,3,1) = in.subMat(0,0,3,1) ;
    out.subMat(0,1,3,1) = in.subMat(3,0,3,1) ;
    out.subMat(0,2,3,1) = in.subMat(6,0,3,1) ;
    out.subMat(0,3,3,1) = in.subMat(9,0,3,1) ;
    return out ;
}


///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::ProcessStatus Area
///
///////////////////////////////////////////////////////////////////////

Matrix TrackerModel2::ProcessStatus::toMat() const
{
    Float mat[24] ;
    int n ;
    if (LH_present)
    {
        n = 24 ;
    }
    else
    {
        n = 17 ;
    }

    mat[ 0]    =   P.x() ;
    mat[ 1]    =   P.y() ;
    mat[ 2]    =   P.z() ;
    mat[ 3]    =   V.x() ;
    mat[ 4]    =   V.y() ;
    mat[ 5]    =   V.z() ;
    mat[ 6]    =   Q.w() ;
    mat[ 7]    =   Q.x() ;
    mat[ 8]    =   Q.y() ;
    mat[ 9]    =   Q.z() ;
    mat[10]    =   IMU_GyroOffset.x() ;
    mat[11]    =   IMU_GyroOffset.y() ;
    mat[12]    =   IMU_GyroOffset.z() ;
    mat[13]    =   IMU_AccOffset.x()  ;
    mat[14]    =   IMU_AccOffset.y()  ;
    mat[15]    =   IMU_AccOffset.z()  ;
    mat[16]    =   Gravity ;

    if (LH_present)
    {
        mat[17] = LHQ.w() ;
        mat[18] = LHQ.x() ;
        mat[19] = LHQ.y() ;
        mat[20] = LHQ.z() ;
        mat[21] = LHT.x() ;
        mat[22] = LHT.y() ;
        mat[23] = LHT.z() ;
    }

    Matrix Mat = Matrix(n,1,mat).clone() ;
    return Mat ;
}

Matrix TrackerModel2::ProcessStatus::toTangentMat() const
{
    Float mat[22] ;
    int n ;
    if (LH_present)
    {
        n = 22 ;
    }
    else
    {
        n = 16 ;
    }

    mat[ 0]    =   P.x() ;
    mat[ 1]    =   P.y() ;
    mat[ 2]    =   P.z() ;
    mat[ 3]    =   V.x() ;
    mat[ 4]    =   V.y() ;
    mat[ 5]    =   V.z() ;
    mat[ 6]    =   Q_Tangent.x() ;
    mat[ 7]    =   Q_Tangent.y() ;
    mat[ 8]    =   Q_Tangent.z() ;
    mat[ 9]    =   IMU_GyroOffset.x() ;
    mat[10]    =   IMU_GyroOffset.y() ;
    mat[11]    =   IMU_GyroOffset.z() ;
    mat[12]    =   IMU_AccOffset.x()  ;
    mat[13]    =   IMU_AccOffset.y()  ;
    mat[14]    =   IMU_AccOffset.z()  ;
    mat[15]    =   Gravity ;

    if (LH_present)
    {
        mat[16] = LHQ_Tangent.x() ;
        mat[17] = LHQ_Tangent.y() ;
        mat[18] = LHQ_Tangent.z() ;
        mat[19] = LHT.x()         ;
        mat[20] = LHT.y()         ;
        mat[21] = LHT.z()         ;
    }

    Matrix Mat = Matrix(n,1,mat).clone() ;
    return Mat ;
}

void TrackerModel2::ProcessStatus::fromMat(const Matrix& M)
{
    Float mat[24] ;

    int n = M.rows() ;
    Matrix(n,1,mat) = M ;
    if (n != 24)
        LH_present = false ;
    else
        LH_present = true ;

                 P.x()    =  mat[ 0]  ;
                 P.y()    =  mat[ 1]  ;
                 P.z()    =  mat[ 2]  ;
                 V.x()    =  mat[ 3]  ;
                 V.y()    =  mat[ 4]  ;
                 V.z()    =  mat[ 5]  ;
                 Q.w()    =  mat[ 6]  ;
                 Q.x()    =  mat[ 7]  ;
                 Q.y()    =  mat[ 8]  ;
                 Q.z()    =  mat[ 9]  ;
    IMU_GyroOffset.x()    =  mat[10]  ;
    IMU_GyroOffset.y()    =  mat[11]  ;
    IMU_GyroOffset.z()    =  mat[12]  ;
    IMU_AccOffset.x()     =  mat[13]  ;
    IMU_AccOffset.y()     =  mat[14]  ;
    IMU_AccOffset.z()     =  mat[15]  ;
               Gravity    =  mat[16]  ;

    if (LH_present)
    {
        LHQ.w() = mat[17] ;
        LHQ.x() = mat[18] ;
        LHQ.y() = mat[19] ;
        LHQ.z() = mat[20] ;
        LHT.x() = mat[21] ;
        LHT.y() = mat[22] ;
        LHT.z() = mat[23] ;
    }
}

void TrackerModel2::ProcessStatus::fromTangentMat(const Matrix& M)
{
    Float mat[22] ;

    int n = M.rows() ;
    Matrix(n,1,mat) = M ;
    if (n != 22)
        LH_present = false ;
    else
        LH_present = true ;

                 P.x()    =  mat[ 0]  ;
                 P.y()    =  mat[ 1]  ;
                 P.z()    =  mat[ 2]  ;
                 V.x()    =  mat[ 3]  ;
                 V.y()    =  mat[ 4]  ;
                 V.z()    =  mat[ 5]  ;
         Q_Tangent.x()    =  mat[ 6]  ;
         Q_Tangent.y()    =  mat[ 7]  ;
         Q_Tangent.z()    =  mat[ 8]  ;
    IMU_GyroOffset.x()    =  mat[ 9]  ;
    IMU_GyroOffset.y()    =  mat[10]  ;
    IMU_GyroOffset.z()    =  mat[11]  ;
    IMU_AccOffset.x()     =  mat[12]  ;
    IMU_AccOffset.y()     =  mat[13]  ;
    IMU_AccOffset.z()     =  mat[14]  ;
               Gravity    =  mat[15]  ;

    if (LH_present)
    {
        LHQ_Tangent.x() = mat[16] ;
        LHQ_Tangent.y() = mat[17] ;
        LHQ_Tangent.z() = mat[18] ;
        LHT.x()         = mat[19] ;
        LHT.y()         = mat[20] ;
        LHT.z()         = mat[21] ;
    }
}

///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::ProcessNoise Area
///
///////////////////////////////////////////////////////////////////////


void TrackerModel2::ProcessNoise::fromMat(const Matrix& M)
{
    Float mat[19] ;
    bool LH_present ;
    int n = M.rows() ;
    Matrix(n,1,mat) = M ;
    if (n != 19)
        LH_present = false ;
    else
        LH_present = true ;

    IMU_GyroNoise.x()          =  mat[ 0]  ;
    IMU_GyroNoise.y()          =  mat[ 1]  ;
    IMU_GyroNoise.z()          =  mat[ 2]  ;
    IMU_AccNoise .x()          =  mat[ 3]  ;
    IMU_AccNoise .y()          =  mat[ 4]  ;
    IMU_AccNoise .z()          =  mat[ 5]  ;
    IMU_GyroOffsetNoise.x()    =  mat[ 6]  ;
    IMU_GyroOffsetNoise.y()    =  mat[ 7]  ;
    IMU_GyroOffsetNoise.z()    =  mat[ 8]  ;
    IMU_AccOffsetNoise.x()     =  mat[ 9]  ;
    IMU_AccOffsetNoise.y()     =  mat[10]  ;
    IMU_AccOffsetNoise.z()     =  mat[11]  ;
    GravityNoise               =  mat[12]  ;

    if (LH_present)
    {
        LHQ_Tangent_Noise.x() = mat [13] ;
        LHQ_Tangent_Noise.y() = mat [14] ;
        LHQ_Tangent_Noise.z() = mat [15] ;
        LHT_Noise.x() = mat[16] ;
        LHT_Noise.y() = mat[17] ;
        LHT_Noise.z() = mat[18] ;
    }
}

///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::LHInfo Area
///
///////////////////////////////////////////////////////////////////////

Matrix TrackerModel2::LHInfo::getTransMat()
{
    if (Value.isEmpty())
        return Value ;

    Matrix Mat(3,4) ;
    Mat.subMat(0,0,3,3) = LHQ.toMat() ;
    Mat.col(3) = LHT ;
    return Mat ;
}

///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::Measurement Area
///
///////////////////////////////////////////////////////////////////////

Matrix TrackerModel2::Measurement::toMat() const
{
    Float mat[7] ;
    int n=7;
    mat[ 0]    =   LH_Pos.x() ;
    mat[ 1]    =   LH_Pos.y() ;
    mat[ 2]    =   LH_Pos.z() ;
    mat[ 3]    =   LH_Ori.w() ;
    mat[ 4]    =   LH_Ori.x() ;
    mat[ 5]    =   LH_Ori.y() ;
    mat[ 6]    =   LH_Ori.z() ;
    Matrix Mat = Matrix(n,1,mat).clone() ;
    return Mat ;
}

Matrix TrackerModel2::Measurement::toTangentMat() const
{
    Float mat[6] ;
    int n=6;
    mat[ 0]    =   LH_Pos.x() ;
    mat[ 1]    =   LH_Pos.y() ;
    mat[ 2]    =   LH_Pos.z() ;
    mat[ 3]    =   LH_Ori_Tangent.x() ;
    mat[ 4]    =   LH_Ori_Tangent.y() ;
    mat[ 5]    =   LH_Ori_Tangent.z() ;
    Matrix Mat = Matrix(n,1,mat).clone() ;
    return Mat ;
}


void TrackerModel2::Measurement::fromMat(const Matrix& M)
{
    Float mat[7] ;
    Matrix(7,1,mat) = M ;
    LH_Pos.x() =   mat[ 0]    ;
    LH_Pos.y() =   mat[ 1]    ;
    LH_Pos.z() =   mat[ 2]    ;
    LH_Ori.w() =   mat[ 3]    ;
    LH_Ori.x() =   mat[ 4]    ;
    LH_Ori.y() =   mat[ 5]    ;
    LH_Ori.z() =   mat[ 6]    ;
}

void TrackerModel2::Measurement::fromTangentMat(const Matrix& M)
{
    Float mat[6] ;
    Matrix(6,1,mat) = M ;
    LH_Pos.x() =   mat[ 0]    ;
    LH_Pos.y() =   mat[ 1]    ;
    LH_Pos.z() =   mat[ 2]    ;
    LH_Ori_Tangent.x() =   mat[ 3]    ;
    LH_Ori_Tangent.y() =   mat[ 4]    ;
    LH_Ori_Tangent.z() =   mat[ 5]    ;
}


///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::MeasurementNoise Area
///
///////////////////////////////////////////////////////////////////////

void TrackerModel2::MeasurementNoise::fromMat(const Matrix& M)
{

}

///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::Control Area
///
///////////////////////////////////////////////////////////////////////

Matrix TrackerModel2::Control::toMat() const
{
    Float mat[7] ;
    int n=7;
    mat[ 0]    =   IMU_Acc.x() ;
    mat[ 1]    =   IMU_Acc.y() ;
    mat[ 2]    =   IMU_Acc.z() ;
    mat[ 3]    =   IMU_Gyro.x() ;
    mat[ 4]    =   IMU_Gyro.y() ;
    mat[ 5]    =   IMU_Gyro.z() ;
    mat[ 6]    =   dt ;
    Matrix Mat = Matrix(n,1,mat).clone() ;
    return Mat ;
}

void TrackerModel2::Control::fromMat(const Matrix& M)
{
    Float mat[7] ;
    Matrix(7,1,mat) = M ;
    IMU_Acc.x()     =  mat[ 0] ;
    IMU_Acc.y()     =  mat[ 1] ;
    IMU_Acc.z()     =  mat[ 2] ;
    IMU_Gyro.x()    =  mat[ 3] ;
    IMU_Gyro.y()    =  mat[ 4] ;
    IMU_Gyro.z()    =  mat[ 5] ;
    dt              =  mat[ 6] ;
}


///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::StatusMap Area
///
///////////////////////////////////////////////////////////////////////

TrackerModel2::StatusMap::StatusMap(TrackerModel2* model, const Matrix& refPoint)
{
    this->model = model ;
    if (! refPoint.isEmpty())
    {
        setRefPoint(refPoint) ;
    }
}

TrackerModel2::StatusMap::~StatusMap()
{

}

Matrix TrackerModel2::StatusMap::calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const
{
    int n = refPointOnManifold.cols() ;
    ProcessStatus status ;
    status.fromMat(refPointOnManifold.col(0));

    Quaternion refQ = status.Q ;
    Quaternion refLHQ = status.LHQ ;
    Float w = Weight(0,0) ;
    Matrix ret = w * refPointOnManifold.col(0) ;

    for (int i=1; i<n; i++)
    {
        w = Weight(i,0) ;
        status.fromMat(refPointOnManifold.col(i));
        Float dot = Quaternion::Dot(status.Q, refQ) ;
        if (dot < 0)
            status.Q = -status.Q ;

        if (status.LH_present)
        {
            if (!LHQ_vaild)
            {
                printf("TrackerModel2::StatusMap::calcBestRefPoint() :\n\tunexpected input with LHInfo\n") ;
            }
            else
            {
                dot = Quaternion::Dot(status.LHQ, refLHQ) ;
                if (dot < 0)
                    status.LHQ = -status.LHQ ;
            }
        }

        ret += w * status.toMat() ;
    }

    return project2Manifold_singleCol(ret) ;
}

extern bool UKF_Print_SR ;
void TrackerModel2::StatusMap::setRefPoint(const Matrix& refPointOnManifold)
{
    ProcessStatus status ;
    status.fromMat(refPointOnManifold);
    Q_Base = getOthoCompBase( Quat2Mat(status.Q), 0.0001 ) ;
    Q_Ref = status.Q ;

    if (status.LH_present)
    {
        LHQ_vaild = true ;
        LHQ_Base = getOthoCompBase(Quat2Mat(status.LHQ), 0.0001) ;
        LHQ_Ref = status.LHQ ;
    }
    else
    {
        LHQ_vaild = false ;
    }

    //if (UKF_Print_SR)
    //{
    //    printf("\n-----------\n");
    //    printf("StatusMap:\n");
    //    Q_Base.transpose().printMat() ;
    //    Matrix M(1,4) ;
    //    M(0,0) = Q_Ref.w() ; M(0,1) = Q_Ref.x() ; M(0,2) = Q_Ref.y() ; M(0,3) = Q_Ref.z() ;
    //    M.printMat() ;
    //    if (status.LH_present)
    //    {
    //        LHQ_Base.transpose().printMat() ;
    //        M(0,0) = LHQ_Ref.w() ; M(0,1) = LHQ_Ref.x() ; M(0,2) = LHQ_Ref.y() ; M(0,3) = LHQ_Ref.z() ;
    //        M.printMat() ;
    //    }
    //    printf("\n-----------\n");
    //}
}

Matrix TrackerModel2::StatusMap::project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const
{
    ProcessStatus status ;
    status.fromMat(pointsOutsideManifold);
    status.Q.normalize();
    if (status.LH_present)
    {
        status.LHQ.normalize();
    }

    return status.toMat() ;
}

Matrix TrackerModel2::StatusMap::map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const
{
    ProcessStatus status ;
    status.fromTangentMat(pointsInTangentSpace);
    Matrix mQ = Q_Base * Matrix(status.Q_Tangent) + Quat2Mat(Q_Ref) ;
    status.Q = Mat2Quat(mQ) ;
    status.Q.normalize();

    if (status.LH_present)
    {
        if (!LHQ_vaild)
        {
            printf("TrackerModel2::StatusMap::map2Manifold_singleCol() :\n\tunexpected input with LHInfo\n") ;
        }
        else
        {
            mQ = LHQ_Base * Matrix(status.LHQ_Tangent) + Quat2Mat(LHQ_Ref) ;
            status.LHQ = Mat2Quat(mQ) ;
            status.LHQ.normalize();
        }
    }
    return status.toMat() ;
}

Matrix TrackerModel2::StatusMap::map2Tangent_singleCol(const Matrix& pointsOnManifold) const
{
    ProcessStatus status ;
    status.fromMat(pointsOnManifold);
    Float dot = Quaternion::Dot(status.Q, Q_Ref) ;
    if (dot < 0)
    {
        dot = -dot ;
        status.Q = -status.Q ;
    }
    status.Q_Tangent  = Q_Base.transpose() * Quat2Mat(status.Q) / dot ;



    if (status.LH_present)
    {
        if (!LHQ_vaild)
        {
            printf("TrackerModel2::StatusMap::map2Tangent_singleCol() :\n\tunexpected input with LHInfo\n") ;
        }
        else
        {
            dot =  Quaternion::Dot(status.LHQ, LHQ_Ref) ;
            if (dot < 0)
            {
                dot = -dot ;
                status.LHQ = -status.LHQ ;
            }
            status.LHQ_Tangent = LHQ_Base.transpose() * Quat2Mat(status.LHQ) / dot ;
        }
    }

    return status.toTangentMat() ;
}



///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2::MeasureMap Area
///
///////////////////////////////////////////////////////////////////////

TrackerModel2::MeasureMap::MeasureMap(TrackerModel2* model, const Matrix& refPoint)
{
    this->model = model ;
    if (! refPoint.isEmpty())
        setRefPoint(refPoint) ;
}

TrackerModel2::MeasureMap::~MeasureMap()
{

}

Matrix TrackerModel2::MeasureMap::calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const
{
    int n = refPointOnManifold.cols() ;
    Measurement measure ;
    measure.fromMat(refPointOnManifold.col(0));

    Quaternion refQ = measure.LH_Ori ;
    Float w = Weight(0,0) ;
    Matrix ret = w * refPointOnManifold.col(0) ;

    for (int i=1; i<n; i++)
    {
        w = Weight(i,0) ;
        measure.fromMat(refPointOnManifold.col(i));
        Float dot = Quaternion::Dot(measure.LH_Ori, refQ) ;
        if (dot < 0)
            measure.LH_Ori = -measure.LH_Ori ;
        ret += w * measure.toMat() ;
    }

    return project2Manifold_singleCol(ret) ;
}


void TrackerModel2::MeasureMap::setRefPoint(const Matrix& refPointOnManifold)
{
    Measurement m ;
    m.fromMat(refPointOnManifold);
    LH_Ori_Base = getOthoCompBase(Quat2Mat(m.LH_Ori), 0.0001) ;
    LH_Ori_Ref = m.LH_Ori ;

    //if (UKF_Print_SR)
    //{
    //    printf("\n-----------\n");
    //    printf("MeasureMap:\n");
    //    LH_Ori_Base.transpose().printMat() ;
    //    Matrix M(1,4) ;
    //    M(0,0) = LH_Ori_Ref.w() ; M(0,1) = LH_Ori_Ref.x() ; M(0,2) = LH_Ori_Ref.y() ; M(0,3) = LH_Ori_Ref.z() ;
    //    M.printMat() ;
    //    printf("\n-----------\n");
    //}

}

Matrix TrackerModel2::MeasureMap::project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const
{
    Measurement m ;
    m.fromMat(pointsOutsideManifold);
    m.LH_Ori.normalize();
    return m.toMat() ;
}

Matrix TrackerModel2::MeasureMap::map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const
{
    Measurement m ;
    m.fromTangentMat(pointsInTangentSpace);
    Matrix mQ = LH_Ori_Base * Matrix(m.LH_Ori_Tangent) + Quat2Mat(LH_Ori_Ref) ;
    m.LH_Ori = Mat2Quat(mQ) ;
    m.LH_Ori.normalize();
    return m.toMat() ;
}

Matrix TrackerModel2::MeasureMap::map2Tangent_singleCol(const Matrix& pointsOnManifold) const
{
    Measurement m ;
    m.fromMat(pointsOnManifold);
    Float dot = Quaternion::Dot(m.LH_Ori, LH_Ori_Ref) ;
    if (dot < 0)
    {
        dot = -dot ;
        m.LH_Ori = -m.LH_Ori ;
    }
    m.LH_Ori_Tangent = LH_Ori_Base.transpose() * Quat2Mat(m.LH_Ori) / dot ;
    return m.toTangentMat() ;
}


///////////////////////////////////////////////////////////////////////
///
/// TrackerModel2 Area
///
///////////////////////////////////////////////////////////////////////


TrackerModel2::TrackerModel2() : UKF::Model()
{

    /*
     *  status partI: Position(3), Velocity(3), Quaternion(4), AngularSpeed(3)
     *  status partII:  IMU_GyroOffset(3), IMU_AccOffset(3), Gravity(1)
     *  status partIII: LH_Mat(12), each station has a single instance of this.
     *
     *  control :  IMU_GyroMeasurement(3), IMU_AccMeasurement(3)
     *
     *  non-additive process noise partI:  IMU_GyroNoise(3, 0.01gyro+0.001 rad/s),  IMU_AccNoise(0.01acc+0.01 m/s2)
     *  non-additive process noise partII:  noise for IMU_parame ?
     *      IMU_GyroOffsetNoise(3, 0.000001 rad/s)  IMU_AccOffsetNoise(3, 0.00001 m/s2)  GravityNoise(1, 0.0000001m/s2)
     *  non-additive process noise partIII:  noise for LH_parame ?
     *      LH_MatPosNoise(3, 0.00001m), LH_MatOriNoise(3, 0.00001 rad),
     *
     *
     *  measurement:  LH_Measured_P(3),  LH_Measured_Q(4)
     *  R_additvie is constant depending on the LightHouse (or non-constant depending on the distance from the tracker to LH).
     *  LH_MeasuredPosNoise(3, 0.02m), LH_MeasuredQuatNoise(4, 0.02 rad)
     */

    statusMap = new StatusMap(this) ;
    measureMap = new MeasureMap(this) ;
    CurStatus = new ProcessStatus ;
    UKF::Model::statusMap = statusMap ;
    UKF::Model::measureMap = measureMap ;

    //// ProcessNoise
    IMU_GyroNoise = 0.001 ;
    IMU_AccNoise = 0.01 ;
    IMU_GyroOffsetNoise = 0.000001 ;
    IMU_AccOffsetNoise = 0.00001 ;
    GravityNoise = 0.0000001 ;
    LH_MatPosNoise = 0.00001 ;
    LH_MatOriNoise = 0.00001 ;

    //// MeasurementNoise
    LH_MeasuredPosNoise = 0.1 ;
    LH_MeasuredQuatNoise = 0.02 ;

    //// InitDeviation
    Position_InitDeviation                 = 0.1 ;
    Velocity_InitDeviation                 = 1 ;
    Quaternion_Tangent_InitDeviation       = 0.05 ;
    IMU_GyroOffset_InitDeviation           = 0.05 ;
    IMU_AccOffset_InitDeviation            = 1 ;
    Gravity_InitDeviation                  = 0.05 ;
    TransMat_RotationPart_InitDeviation    = 0.05 ;
    TransMat_TranslationPart_InitDeviation = 0.05 ;

    Initialized = false ;
}

TrackerModel2::~TrackerModel2()
{
    delete statusMap ;
    delete measureMap ;
    delete CurStatus ;
}

bool TrackerModel2::importParam (const char* fn)
{
    Matrix B ;
    B.LoadMatFile(fn) ;
    return importParam(B.subMat(0,0,17,1)) ;
}

bool TrackerModel2::setParam(const Float* p)
{
     int i = 0;
    //// ProcessNoise
    IMU_GyroNoise = p[i++] ;
    IMU_AccNoise = p[i++] ;
    IMU_GyroOffsetNoise = p[i++] ;
    IMU_AccOffsetNoise = p[i++] ;
    GravityNoise = p[i++] ;
    LH_MatPosNoise = p[i++] ;
    LH_MatOriNoise = p[i++] ;

    //// MeasurementNoise
    LH_MeasuredPosNoise = p[i++] ;
    LH_MeasuredQuatNoise = p[i++] ;

    //// InitDeviation
    Position_InitDeviation                 = p[i++] ;
    Velocity_InitDeviation                 = p[i++] ;
    Quaternion_Tangent_InitDeviation       = p[i++] ;
    IMU_GyroOffset_InitDeviation           = p[i++] ;
    IMU_AccOffset_InitDeviation            = p[i++] ;
    Gravity_InitDeviation                  = p[i++] ;
    TransMat_RotationPart_InitDeviation    = p[i++] ;
    TransMat_TranslationPart_InitDeviation = p[i++] ;

    return true ;
}

bool TrackerModel2::exportParam(const char* fn)
{
    FILE* fp = fopen (fn, "wb") ;
    if (!fp)
        return false ;

    Float mat[17] ;
    int i = 0 ;

    //// ProcessNoise
    mat[i++] = IMU_GyroNoise        ;
    mat[i++] = IMU_AccNoise         ;
    mat[i++] = IMU_GyroOffsetNoise  ;
    mat[i++] = IMU_AccOffsetNoise   ;
    mat[i++] = GravityNoise         ;
    mat[i++] = LH_MatPosNoise       ;
    mat[i++] = LH_MatOriNoise       ;

    //// MeasurementNoise
    mat[i++] = LH_MeasuredPosNoise  ;
    mat[i++] = LH_MeasuredQuatNoise ;

    //// InitDeviation
    mat[i++] = Position_InitDeviation                 ;
    mat[i++] = Velocity_InitDeviation                 ;
    mat[i++] = Quaternion_Tangent_InitDeviation       ;
    mat[i++] = IMU_GyroOffset_InitDeviation           ;
    mat[i++] = IMU_AccOffset_InitDeviation            ;
    mat[i++] = Gravity_InitDeviation                  ;
    mat[i++] = TransMat_RotationPart_InitDeviation    ;
    mat[i++] = TransMat_TranslationPart_InitDeviation ;

    for (i=0; i<17; i++)
    {
        fprintf(fp, "%11.12f\n", mat[i]) ;
    }
    fclose(fp) ;

    return true ;
}

bool TrackerModel2::importParam (const Matrix& M)
{
    if (M.rows() != 17 || M.cols() != 1)
    {
        printf("TrackerModel2::importParam() : Input ParamMat has wrong dimension\n") ;
        return false ;
    }

    int i=0 ;
    Float mat[17] ;
    Matrix m(17,1,mat) ;
    m = M ;

    //// ProcessNoise
    IMU_GyroNoise = mat[i++] ;
    IMU_AccNoise = mat[i++] ;
    IMU_GyroOffsetNoise = mat[i++] ;
    IMU_AccOffsetNoise = mat[i++] ;
    GravityNoise = mat[i++] ;
    LH_MatPosNoise = mat[i++] ;
    LH_MatOriNoise = mat[i++] ;

    //// MeasurementNoise
    LH_MeasuredPosNoise = mat[i++] ;
    LH_MeasuredQuatNoise = mat[i++] ;

    //// InitDeviation
    Position_InitDeviation                 = mat[i++] ;
    Velocity_InitDeviation                 = mat[i++] ;
    Quaternion_Tangent_InitDeviation       = mat[i++] ;
    IMU_GyroOffset_InitDeviation           = mat[i++] ;
    IMU_AccOffset_InitDeviation            = mat[i++] ;
    Gravity_InitDeviation                  = mat[i++] ;
    TransMat_RotationPart_InitDeviation    = mat[i++] ;
    TransMat_TranslationPart_InitDeviation = mat[i++] ;

    return true ;
}

Matrix TrackerModel2::exportParam()
{
    Float mat[17] ;
    int i = 0 ;

    //// ProcessNoise
    mat[i++] = IMU_GyroNoise        ;
    mat[i++] = IMU_AccNoise         ;
    mat[i++] = IMU_GyroOffsetNoise  ;
    mat[i++] = IMU_AccOffsetNoise   ;
    mat[i++] = GravityNoise         ;
    mat[i++] = LH_MatPosNoise       ;
    mat[i++] = LH_MatOriNoise       ;

    //// MeasurementNoise
    mat[i++] = LH_MeasuredPosNoise  ;
    mat[i++] = LH_MeasuredQuatNoise ;

    //// InitDeviation
    mat[i++] = Position_InitDeviation                 ;
    mat[i++] = Velocity_InitDeviation                 ;
    mat[i++] = Quaternion_Tangent_InitDeviation       ;
    mat[i++] = IMU_GyroOffset_InitDeviation           ;
    mat[i++] = IMU_AccOffset_InitDeviation            ;
    mat[i++] = Gravity_InitDeviation                  ;
    mat[i++] = TransMat_RotationPart_InitDeviation    ;
    mat[i++] = TransMat_TranslationPart_InitDeviation ;

    Matrix M (17,1,mat) ;
    return M.clone() ;
}



bool TrackerModel2::setPosterior(const Matrix& X, const Matrix& Pxx)
{
    if ( Model::setPosterior(X, Pxx) )
    {
        CurStatus->fromMat(getStatus());
        return true ;
    }

    else
    {
        return false ;
    }
}

const TrackerModel2::ProcessStatus* TrackerModel2::getCurStatus() const
{
    if (isInitialized())
        return CurStatus ;
    else
        return NULL ;
}



void TrackerModel2::reset()
{
    Initialized = false ;
}

bool TrackerModel2::isInitialized() const
{
    return Initialized ;
}

void TrackerModel2::init (const Vector3D& P, const Quaternion& Q, const Matrix& TransMat)
{
    ProcessStatus status ;
    status.Gravity = 9.8 ;
    status.IMU_AccOffset = Vector3D(0,0,0) ;
    status.IMU_GyroOffset = Vector3D(0,0,0) ;

    Matrix R = TransMat.subMat(0,0,3,3) ;
    Vector3D T = TransMat.col(3) ;

    status.LHQ.fromMat(R) ;
    status.LHQ.normalize();
    status.LHQ_Tangent = Vector3D(0,0,0) ;
    status.LHT = T ;
    status.LH_present = true ;

    Vector3D mP = Vector3D(R * P) + T ;
    Quaternion mQ ;
    mQ.fromMat(R) ;
    mQ = mQ * Q ;

    status.P = mP ;
    status.Q = mQ ;
    status.Q_Tangent = Vector3D(0,0,0) ;
    status.V = Vector3D(0,0,0) ;

    statusMap->setRefPoint(status.toMat());

    Matrix Pxx(22,22) ; Pxx.fill(0) ;
    Pxx( 0, 0) = Pxx( 1, 1) = Pxx( 2, 2) = Position_InitDeviation * Position_InitDeviation ;
    Pxx( 3, 3) = Pxx( 4, 4) = Pxx( 5, 5) = Velocity_InitDeviation * Velocity_InitDeviation ;
    Pxx( 6, 6) = Pxx( 7, 7) = Pxx( 8, 8) = Quaternion_Tangent_InitDeviation * Quaternion_Tangent_InitDeviation ;
    Pxx( 9, 9) = Pxx(10,10) = Pxx(11,11) = IMU_GyroOffset_InitDeviation * IMU_GyroOffset_InitDeviation ;
    Pxx(12,12) = Pxx(13,13) = Pxx(14,14) = IMU_AccOffset_InitDeviation * IMU_AccOffset_InitDeviation ;
    Pxx(15,15) = Gravity_InitDeviation * Gravity_InitDeviation ;
    Pxx(16,16) = Pxx(17,17) = Pxx(18,18) = TransMat_RotationPart_InitDeviation * TransMat_RotationPart_InitDeviation ;
    Pxx(19,19) = Pxx(20,20) = Pxx(21,21) = TransMat_TranslationPart_InitDeviation * TransMat_TranslationPart_InitDeviation ;

    setPosterior(status.toTangentMat(), Pxx) ;


    Q_LHpresent.type = NonAdditive_Noise ;
    Q_LHpresent.nonAdditiveCov = Matrix(19,19) ;
    Q_LHpresent.SR_nonAdditiveCov = Matrix(19,19) ;
    Q_LHpresent.nonAdditiveCov.fill(0);
    Q_LHpresent.SR_nonAdditiveCov.fill(0);

    this->Q.type = NonAdditive_Noise ;
    //this->Q.nonAdditiveCov = Matrix(19,19) ;
    //this->Q.SR_nonAdditiveCov = Matrix(19,19) ;
    this->Q.nonAdditiveCov = Matrix(13,13) ;
    this->Q.SR_nonAdditiveCov = Matrix(13,13) ;
    this->Q.nonAdditiveCov.fill(0);
    this->Q.SR_nonAdditiveCov.fill(0);

    this->R.type = PureAdditive_Noise ;
    this->R.additiveCov = Matrix(6,6) ;
    this->R.additiveCov.fill(0) ;

    updateQ(false, true) ;
    updateQ(true, true) ;
    Initialized = true ;
}

void TrackerModel2::updateQ(bool LH_present, bool updateValue)
{
    if (LH_present)
    {
        this->curQ = &Q_LHpresent ;
    }
    else
    {
        this->curQ = &(this->Q) ;
    }

    if (!updateValue)
        return ;

    Matrix& Q = curQ->nonAdditiveCov ;
    Matrix& SR = curQ->SR_nonAdditiveCov ;

    Q(0,0) = Q(1,1) = Q(2,2) = IMU_GyroNoise*IMU_GyroNoise ;
    Q(3,3) = Q(4,4) = Q(5,5) = IMU_AccNoise*IMU_AccNoise ;
    Q(6,6) = Q(7,7) = Q(8,8) = IMU_GyroOffsetNoise*IMU_GyroOffsetNoise ;
    Q(9,9) = Q(10,10) = Q(11,11) = IMU_AccOffsetNoise*IMU_AccOffsetNoise ;
    Q(12,12) = GravityNoise*GravityNoise ;

    SR(0,0) = SR(1,1) = SR(2,2) = IMU_GyroNoise ;
    SR(3,3) = SR(4,4) = SR(5,5) = IMU_AccNoise ;
    SR(6,6) = SR(7,7) = SR(8,8) = IMU_GyroOffsetNoise ;
    SR(9,9) = SR(10,10) = SR(11,11) = IMU_AccOffsetNoise ;
    SR(12,12) = GravityNoise ;

    if (LH_present)
    {
        Q(13,13) = Q(14,14) = Q(15,15) = LH_MatOriNoise*LH_MatOriNoise ;
        Q(16,16) = Q(17,17) = Q(18,18) = LH_MatPosNoise*LH_MatPosNoise ;
        SR(13,13) = SR(14,14) = SR(15,15) = LH_MatOriNoise ;
        SR(16,16) = SR(17,17) = SR(18,18) = LH_MatPosNoise ;
    }
}

void TrackerModel2::updateR(const Vector3D& refZ)
{
    Matrix& R = curR->additiveCov ;
    Matrix& SR_R = curR->SR_additiveCov ;
    Vector3D X,Y,Z;
    Z = refZ.getNormalized() ;
    if (Z.x() < Z.y())
        Y = Vector3D(1,0,0) ;
    else
        Y = Vector3D(0,1,0) ;
    X = Vector3D::Cross(Y,Z).getNormalized() ;
    Y = Vector3D::Cross(Z,X).getNormalized() ;
    const Float e = 0.02 * 3.1415926 / 180 ;
    const Float l = 0.04 ; // the average dimension(of radius) of trackers.
    Float s = refZ.getLength() ;
    Float errXY, errZ ;
    errXY = e * s ;
    errZ  = LH_MeasuredPosNoise * e * s * s / l ;
    Matrix S(3,3) ;
    S.col(0) = X * errXY ;
    S.col(1) = Y * errXY ;
    S.col(2) = Z * errZ ;

    Matrix posCov = S * S.transpose() ;
    R.subMat(0,0,3,3) = posCov ;
    R(3,3) = R(4,4) = R(5,5) = LH_MeasuredQuatNoise*LH_MeasuredQuatNoise ;

    SR_R = Matrix(6,6) ;
    SR_R.fill(0);
    SR_R.subMat(0,0,3,3) = SR(posCov) ;
    SR_R(3,3) = SR_R(4,4) = SR_R(5,5) = LH_MeasuredQuatNoise ;
}



void TrackerModel2::setLHInfo (const LHInfo& info)
{
    if ( !info.Cov.isEmpty() )
    {
        Matrix X(22, 1) ;
        Matrix Pxx(22,22) ;
        X.subMat(0,0,16,1) = getPosterior().X.subMat(0,0,16,1) ;
        Pxx.fill(0);
        Pxx.subMat(0,0,16,16) = getPosterior().Pxx.subMat(0,0,16,16) ;

        X.subMat(16,0,6,1) = info.Value ;
        Pxx.subMat(16,16,6,6) = info.Cov ;

        statusMap->LHQ_Base = info.LHQ_Base ;
        statusMap->LHQ_Ref = info.LHQ_Ref ;
        statusMap->LHQ_vaild = true ;

        setPosterior(X, Pxx) ;
    }
    else
    {
        Matrix X(16, 1) ;
        Matrix Pxx(16,16) ;
        X = getPosterior().X.subMat(0,0,16,1) ;
        Pxx = getPosterior().Pxx.subMat(0,0,16,16) ;
        statusMap->LHQ_vaild = false ;
        setPosterior(X, Pxx) ;

    }
}

TrackerModel2::LHInfo TrackerModel2::getLHInfo ()
{
    LHInfo info ;
    if (getPosterior().X.rows() != 22)
    {
        info.Cov = Matrix() ;
        info.Value = Matrix() ;
    }
    else
    {
        const ProcessStatus* pStatus = getCurStatus() ;
        if (!pStatus)
        {
            info.Cov = Matrix() ;
            info.Value = Matrix() ;
        }
        else
        {
            info.Cov = getPosterior().Pxx.subMat(16,16,6,6).clone() ;
            info.Value = getPosterior().X.subMat(16,0,6,1).clone() ;
            info.LHQ = pStatus->LHQ ;
            info.LHT = pStatus->LHT ;
            info.LHQ_Ref = statusMap->LHQ_Ref ;
            info.LHQ_Base = statusMap->LHQ_Base ;
        }
    }
    return info ;
}

void TrackerModel2::initLHInfo (TrackerModel2::LHInfo& lhInfo, Matrix& Mat)
{
    lhInfo.LHQ.fromMat(Mat.subMat(0,0,3,3));
    lhInfo.LHT = Mat.col(3);
    lhInfo.Value = Matrix(6,1);
    lhInfo.Value.fill(0);
    lhInfo.Value.subMat(3,0,3,1) = Mat.col(3);
    lhInfo.Cov = Matrix(6,6);
    lhInfo.Cov.fill(0);
    lhInfo.Cov(0,0) = lhInfo.Cov(1,1) = lhInfo.Cov(2,2) = TransMat_RotationPart_InitDeviation * TransMat_RotationPart_InitDeviation ;
    lhInfo.Cov(3,3) = lhInfo.Cov(4,4) = lhInfo.Cov(5,5) = TransMat_TranslationPart_InitDeviation * TransMat_TranslationPart_InitDeviation ;
    lhInfo.LHQ_Ref = lhInfo.LHQ;
    lhInfo.LHQ_Base = getOthoCompBase(Quat2Mat(lhInfo.LHQ), 0.0001);
    //printf(">> initLHInfo:\n");
    //printf(">> LHQ_Ref:%f %f%f %f \n",lhInfo.LHQ.w(),lhInfo.LHQ.x(),lhInfo.LHQ.y(),lhInfo.LHQ.z());
}

void TrackerModel2::CL_UKF_Info(TrackerModel_CL_UKF_Info* info, int offset,
                bool measure_valid )
{
    info->alpha = 0.03 ;
    info->beta = 2 ;
    info->kappa = 0 ;

    if (measure_valid)
    {
        info->constraint_adapt_aux_buf_size_per_item = 4 ;
        info->constraint_adapt_data_size = 16 * 3 ;
        info->control_dim = 8 ;
        info->data_offset = offset ;
        info->measurement_dim = 6 ;
        info->measurement_noise_cov_sr_ready = true ;
        info->measurement_noise_dim = 6 ;
        info->measurement_raw_dim = 7 ;

        info->process_noise_cov_sr_ready = true ;
        info->process_noise_dim = 19 ;
        info->process_status_cov_sr_ready = false ;
        info->process_status_dim = 22 ;
    }

    else
    {
        info->constraint_adapt_aux_buf_size_per_item = 2 ;
        info->constraint_adapt_data_size = 16 ;
        info->control_dim = 8 ;
        info->data_offset = offset ;
        info->measurement_dim = 0 ;
        info->measurement_noise_cov_sr_ready = false ;
        info->measurement_noise_dim = 0 ;
        info->measurement_raw_dim = 0 ;

        info->process_noise_cov_sr_ready = true ;
        info->process_noise_dim = 13 ;
        info->process_status_cov_sr_ready = false ;
        info->process_status_dim = 16 ;
    }
}


int TrackerModel2::CL_UKF_Data_Size (const TrackerModel_CL_UKF_Info* info)
{
    int size = info->process_status_dim +
            info->control_dim +
            info->constraint_adapt_data_size +
            (info->process_status_dim * info->process_status_dim) +
            (info->process_noise_dim * info->process_noise_dim) +
            info->measurement_raw_dim +
            (info->measurement_dim * info->measurement_dim) ;

    return size ;
}


void TrackerModel2::CL_UKF_Data (TrackerModel_CL_UKF_Info* info, float* data,
               const Matrix& control, const Matrix& measurement_raw, bool measure_valid)
{
    Matrix(info->process_status_dim, 1, data) = Posterior.X ;
    data += info->process_status_dim ;

    Matrix(info->control_dim-1, 1, data) = control ;
    data [info->control_dim-1] = measure_valid ;
    data += info->control_dim ;

    if (measure_valid)
    {
        Matrix(3,4,data) = statusMap->Q_Base.transpose() ;
        data += 12 ;
        data[0] = statusMap->Q_Ref.w() ;
        data[1] = statusMap->Q_Ref.x() ;
        data[2] = statusMap->Q_Ref.y() ;
        data[3] = statusMap->Q_Ref.z() ;
        data += 4 ;

        Matrix(3,4,data) = statusMap->LHQ_Base.transpose() ;
        data += 12 ;
        data[0] = statusMap->LHQ_Ref.w() ;
        data[1] = statusMap->LHQ_Ref.x() ;
        data[2] = statusMap->LHQ_Ref.y() ;
        data[3] = statusMap->LHQ_Ref.z() ;
        data += 4 ;

        Matrix(3,4,data) = measureMap->LH_Ori_Base.transpose() ;
        data += 12 ;
        data[0] = measureMap->LH_Ori_Ref.w() ;
        data[1] = measureMap->LH_Ori_Ref.w() ;
        data[2] = measureMap->LH_Ori_Ref.w() ;
        data[3] = measureMap->LH_Ori_Ref.w() ;
        data += 4 ;
    }
    else
    {
        Matrix(3,4,data) = statusMap->Q_Base.transpose() ;
        data += 12 ;
        data[0] = statusMap->Q_Ref.w() ;
        data[1] = statusMap->Q_Ref.x() ;
        data[2] = statusMap->Q_Ref.y() ;
        data[3] = statusMap->Q_Ref.z() ;
        data += 4 ;
    }

    Matrix(info->process_status_dim, info->process_status_dim, data) = Posterior.Pxx ;
    data += info->process_status_dim * info->process_status_dim ;

    Matrix(info->process_noise_dim, info->process_noise_dim, data) = curQ->SR_nonAdditiveCov.transpose() ;
    data += info->process_noise_dim * info->process_noise_dim ;

    if (measure_valid)
    {
        Matrix(info->measurement_raw_dim, 1, data) = measurement_raw ;
        data += info->measurement_raw_dim ;

        Matrix(info->measurement_noise_dim, info->measurement_noise_dim, data) = curR->SR_additiveCov.transpose() ;
        data += info->measurement_noise_dim * info->measurement_noise_dim ;
    }
}

int TrackerModel2::CL_UKF_Result (TrackerModel_CL_UKF_Info* info, float* data)
{
    Matrix X, Pxx ;
    X = Matrix(info->process_status_dim, 1, data).clone() ;
    data += info->process_status_dim ;
    //printf("CL_UKF_Result(): X \n");
    //X.transpose().printMat() ;
    bool measure_valid = (int)(data [info->control_dim-1]) ;
    data += info->control_dim ;

    if (measure_valid)
    {
        statusMap->Q_Base = Matrix(3,4,data).transpose().clone() ;
        //printf("CL_UKF_Result(): statusMap->Q_Base \n");
        //statusMap->Q_Base.transpose().printMat() ;

        data += 12 ;
        statusMap->Q_Ref.w() = data[0] ;
        statusMap->Q_Ref.x() = data[1] ;
        statusMap->Q_Ref.y() = data[2] ;
        statusMap->Q_Ref.z() = data[3] ;
        data += 4 ;
        //printf("CL_UKF_Result(): statusMap->Q_Ref (%f, %f, %f, %f) \n",statusMap->Q_Ref.w(), statusMap->Q_Ref.x(), statusMap->Q_Ref.y(), statusMap->Q_Ref.z() );

        statusMap->LHQ_vaild = true ;
        statusMap->LHQ_Base = Matrix(3,4,data).transpose().clone() ;
        //printf("CL_UKF_Result(): statusMap->LHQ_Base \n");
        //statusMap->LHQ_Base.transpose().printMat() ;


        data += 12 ;
        statusMap->LHQ_Ref.w() = data[0] ;
        statusMap->LHQ_Ref.x() = data[1] ;
        statusMap->LHQ_Ref.y() = data[2] ;
        statusMap->LHQ_Ref.z() = data[3] ;
        data += 4 ;
        //printf("CL_UKF_Result(): statusMap->LHQ_Ref (%f, %f, %f, %f) \n",statusMap->LHQ_Ref.w(), statusMap->LHQ_Ref.x(), statusMap->LHQ_Ref.y(), statusMap->LHQ_Ref.z() );


        measureMap->LH_Ori_Base = Matrix(3,4,data).transpose().clone() ;
        data += 12 ;
        measureMap->LH_Ori_Ref.w() = data[0] ;
        measureMap->LH_Ori_Ref.x() = data[1] ;
        measureMap->LH_Ori_Ref.y() = data[2] ;
        measureMap->LH_Ori_Ref.z() = data[3] ;
        data += 4 ;
    }
    else
    {
        statusMap->Q_Base = Matrix(3,4,data).transpose().clone() ;
        data += 12 ;
        statusMap->Q_Ref.w() = data[0] ;
        statusMap->Q_Ref.x() = data[1] ;
        statusMap->Q_Ref.y() = data[2] ;
        statusMap->Q_Ref.z() = data[3] ;
        data += 4 ;

        statusMap->LHQ_vaild = false ;
    }

    Pxx = Matrix(info->process_status_dim, info->process_status_dim, data) ;
    data += info->process_status_dim * info->process_status_dim ;

    //printf("CL_UKF_Result(): Pxx \n");
    //Pxx.printMat() ;

    setPosterior(X, Pxx) ;
    return 0 ;
}


Matrix TrackerModel2::f(const Matrix& Xk, const Matrix& Uk, const Matrix& nonAdditive_Vk)
{
    //printf(" **** Xk : \n") ;
    //Xk.transpose().printMat() ;

    ProcessStatus X ;
    Control U ;
    ProcessNoise V ;
    X.fromMat(Xk);
    U.fromMat(Uk);
    V.fromMat(nonAdditive_Vk);


    // update additive components
    X.IMU_GyroOffset += V.IMU_GyroOffsetNoise ;
    X.IMU_AccOffset  += V.IMU_AccOffsetNoise ;
    X.Gravity        += V.GravityNoise ;

    if (X.LH_present)
    {
        Matrix Xk_tan = statusMap->map2Tangent(Xk) ;
        ProcessStatus Xtan ;
        Xtan.fromTangentMat(Xk_tan);
        Xtan.LHQ_Tangent += V.LHQ_Tangent_Noise ;
        X.LHQ = Mat2Quat(statusMap->LHQ_Base * Xtan.LHQ_Tangent) + statusMap->LHQ_Ref ;
        X.LHQ.normalize();
        X.LHT  += V.LHT_Noise ;
    }

    // update X.Q, X.P and X.V
    Vector3D  angular, acc, vel ;
    Quaternion dQ, newQ ;
    angular = (U.IMU_Gyro + V.IMU_GyroNoise) - X.IMU_GyroOffset ;
    dQ = Quaternion(1, 0.5*angular.x()*U.dt, 0.5*angular.y()*U.dt, 0.5*angular.z()*U.dt) ;
    newQ = X.Q * dQ ;
    newQ.normalize();

    acc = (U.IMU_Acc + V.IMU_AccNoise) - X.IMU_AccOffset ;
    acc = Quaternion::rotVbyQ(acc, X.Q) + Quaternion::rotVbyQ(acc, newQ) ;
    acc *= 0.5 ;
    acc.y() -= X.Gravity ;

    X.Q = newQ ;
    X.P = X.P + X.V * U.dt + 0.5*acc*U.dt*U.dt ;
    X.V = X.V + acc*U.dt ;

    return X.toMat() ;
}

Matrix TrackerModel2::h(const Matrix& Xpre, const Matrix& nonAdditive_Nk)
{
    ProcessStatus X ;
    MeasurementNoise N ;
    Measurement Y ;
    X.fromMat(Xpre);
    N.fromMat(nonAdditive_Nk);

    Y.LH_Ori = X.LHQ.conj() * X.Q ;
    Y.LH_Ori.normalize();
    Y.LH_Pos = Quaternion::rotVbyQ( (X.P - X.LHT), X.LHQ.conj() ) ;
    return Y.toMat() ;
}
