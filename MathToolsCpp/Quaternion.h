
#ifndef __MathToolsCpp__Quaternion__
#define __MathToolsCpp__Quaternion__

#include "MathToolsCpp/cppFloat.h"
#include "MathToolsCpp/Vector3D.h"
#include "MathToolsCpp/Matrix.h"

START_JK_Math_NAMESPACE
class Quaternion ;
END_JK_Math_NAMESPACE


JKMath::Quaternion operator* (const JKMath::Quaternion& lhs, const JKMath::Quaternion& rhs) ;
JKMath::Quaternion operator+ (const JKMath::Quaternion& lhs, const JKMath::Quaternion& rhs) ;
JKMath::Quaternion operator- (const JKMath::Quaternion& lhs, const JKMath::Quaternion& rhs) ;
JKMath::Quaternion operator- (const JKMath::Quaternion& rhs) ;
JKMath::Quaternion operator* (const JKMath::Quaternion& lhs, Float r) ;
JKMath::Quaternion operator* (Float r, const JKMath::Quaternion& rhs) ;
JKMath::Quaternion operator/ (const JKMath::Quaternion& lhs, Float r) ;


START_JK_Math_NAMESPACE


enum EulerOrder {
       EulerOrder_YXZ,
       EulerOrder_ZYX,
       EulerOrder_XZY,
       EulerOrder_XYZ,
       EulerOrder_YZX,
       EulerOrder_ZXY
};

class EulerAngle ;
class AxisAngle ;
class Quaternion {

public:
    Quaternion() ;
    Quaternion(Float w, Float x, Float y, Float z) ;
    Quaternion(const Float quaternion[]) ;
    Quaternion(const Quaternion& other ) ;
    Quaternion & operator=(const Quaternion& other) ;
    Quaternion & operator+=(const Quaternion& rhs) ;
    Quaternion & operator-=(const Quaternion& rhs) ;
    Quaternion & operator*=(Float r) ;
    Quaternion & operator/=(Float r) ;

    Float& w() ;
    Float& x() ;
    Float& y() ;
    Float& z() ;
    Float w() const ;
    Float x() const ;
    Float y() const ;
    Float z() const ;

    friend Quaternion (::operator*) (const Quaternion& lhs, const Quaternion& rhs) ;
    friend Quaternion (::operator+) (const Quaternion& lhs, const Quaternion& rhs) ;
    friend Quaternion (::operator-) (const Quaternion& lhs, const Quaternion& rhs) ;
    friend Quaternion (::operator-) (const Quaternion& rhs) ;
    friend Quaternion (::operator*) (const Quaternion& lhs, Float r) ;
    friend Quaternion (::operator*) (Float r, const Quaternion& rhs) ;
    friend Quaternion (::operator/) (const Quaternion& lhs, Float r) ;

    void   normalize() ;
    Quaternion conj() const ;

    Vector3D   rotVector (const Vector3D& v) const ;
    EulerAngle toEuler (EulerOrder order = EulerOrder_YXZ) const ;
    AxisAngle  toAxisAngle() const ;
    Matrix toMat () const ;
    bool fromMat (const Matrix& M) ;

public :
    static Quaternion v2q(const Vector3D& vec1, const Vector3D& vec2);
    static Vector3D rotVbyQ (const Vector3D& v, const Quaternion& q);
    static Quaternion Slerp(const Quaternion& start, const Quaternion& end, const Float endWeight);
    static Float Dot (const Quaternion& q1, const Quaternion& q2) ;

private:
    Float values[4] ;
};


class EulerAngle
{

public :
    EulerAngle() ;
    EulerAngle(Float yaw, Float pitch, Float roll) ;

public :

    Float& yaw() ;
    Float& pitch() ;
    Float& roll() ;
    Float yaw() const ;
    Float pitch() const ;
    Float roll() const ;

    EulerAngle DegreeToRad () const ;
    EulerAngle RadToDegree () const ;
    Quaternion toQuat(EulerOrder order=EulerOrder_YXZ) const;

private:
    Float Yaw, Pitch, Roll ;
};


class AxisAngle {

public :
    AxisAngle() ;
    AxisAngle ( const Vector3D& axis, Float angle ) ;

public :

    Vector3D& axis() ;
    Float& angle() ;
    const Vector3D& axis() const ;
    Float angle() const ;


    AxisAngle DegreeToRad () const ;
    AxisAngle RadToDegree () const ;
    Quaternion toQuat () const ;

private:
    Vector3D Axis ;
    Float Angle ;
};

END_JK_Math_NAMESPACE

#endif
