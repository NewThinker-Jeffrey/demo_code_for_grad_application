
#include <cmath>
#include "MathToolsCpp/Quaternion.h"
IMPORT_JK_Math

static const Float _1_ = 0.99999999 ;
static const Float _0_ = 0.00000001 ;
static const Float R2D = 180.0 / 3.1415926 ;
static const Float D2R = 3.1415926 / 180.0 ;



Quaternion::Quaternion() {
    values[0] = 1;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}
Quaternion::Quaternion(Float w, Float x, Float y, Float z) {
    values[0] = w;
    values[1] = x;
    values[2] = y;
    values[3] = z;
}
Quaternion::Quaternion(const Float quaternion[]) {
    values[0] = quaternion[0];
    values[1] = quaternion[1];
    values[2] = quaternion[2];
    values[3] = quaternion[3];
}

Quaternion::Quaternion(const Quaternion& other ) {
    values[0] = other.values[0];
    values[1] = other.values[1];
    values[2] = other.values[2];
    values[3] = other.values[3];
}


void Quaternion::normalize() {
    Float tmp = 1 / sqrt ( w()*w() + x()*x() +
                           y()*y() + z()*z() ) ;
    values[0] *= tmp ;
    values[1] *= tmp ;
    values[2] *= tmp ;
    values[3] *= tmp ;
}

Quaternion & Quaternion::operator=(const Quaternion& rhs) {
    values[0] = rhs.values[0] ;
    values[1] = rhs.values[1] ;
    values[2] = rhs.values[2] ;
    values[3] = rhs.values[3] ;
    return *this ;
}

Quaternion & Quaternion::operator+=(const Quaternion& rhs) {
    values[0] += rhs.values[0] ;
    values[1] += rhs.values[1] ;
    values[2] += rhs.values[2] ;
    values[3] += rhs.values[3] ;
    return *this ;
}

Quaternion & Quaternion::operator-=(const Quaternion& rhs) {
    values[0] -= rhs.values[0] ;
    values[1] -= rhs.values[1] ;
    values[2] -= rhs.values[2] ;
    values[3] -= rhs.values[3] ;
    return *this ;
}

Quaternion & Quaternion::operator*=(Float r) {
    values[0] *= r ;
    values[1] *= r ;
    values[2] *= r ;
    values[3] *= r ;
    return *this ;
}

Quaternion & Quaternion::operator/=(Float r) {
    values[0] /= r ;
    values[1] /= r ;
    values[2] /= r ;
    values[3] /= r ;
    return *this ;
}

Float& Quaternion::w()
{
    return values[0] ;
}

Float& Quaternion::x()
{
    return values[1] ;
}

Float& Quaternion::y()
{
    return values[2] ;
}

Float& Quaternion::z()
{
    return values[3] ;
}

Float Quaternion::w() const
{
    return values[0] ;
}

Float Quaternion::x() const
{
    return values[1] ;
}
Float Quaternion::y() const
{
    return values[2] ;
}
Float Quaternion::z() const
{
    return values[3] ;
}



Quaternion operator*(const Quaternion& lhs, const Quaternion& rhs){

    Quaternion rst;
    rst.values[0] = lhs.values[0] * rhs.values[0] - lhs.values[1] * rhs.values[1] - lhs.values[2] * rhs.values[2] - lhs.values[3] * rhs.values[3];

    rst.values[1] = lhs.values[0] * rhs.values[1] + lhs.values[1] * rhs.values[0] + lhs.values[2] * rhs.values[3] - lhs.values[3] * rhs.values[2];

    rst.values[2] = lhs.values[0] * rhs.values[2] - lhs.values[1] * rhs.values[3] + lhs.values[2] * rhs.values[0] + lhs.values[3] * rhs.values[1];

    rst.values[3] = lhs.values[0] * rhs.values[3] + lhs.values[1] * rhs.values[2] - lhs.values[2] * rhs.values[1] + lhs.values[3] * rhs.values[0];
    return rst;
}

Quaternion operator+(const Quaternion& lhs, const Quaternion& rhs)
{
    Quaternion rst;
    rst.values[0] = lhs.values[0] + rhs.values[0] ;
    rst.values[1] = lhs.values[1] + rhs.values[1] ;
    rst.values[2] = lhs.values[2] + rhs.values[2] ;
    rst.values[3] = lhs.values[3] + rhs.values[3] ;
    return rst;
}

Quaternion operator-(const Quaternion& lhs, const Quaternion& rhs)
{
    Quaternion rst;
    rst.values[0] = lhs.values[0] - rhs.values[0] ;
    rst.values[1] = lhs.values[1] - rhs.values[1] ;
    rst.values[2] = lhs.values[2] - rhs.values[2] ;
    rst.values[3] = lhs.values[3] - rhs.values[3] ;
    return rst;
}

Quaternion operator-(const Quaternion& rhs)
{
    Quaternion rst;
    rst.values[0] = - rhs.values[0] ;
    rst.values[1] = - rhs.values[1] ;
    rst.values[2] = - rhs.values[2] ;
    rst.values[3] = - rhs.values[3] ;
    return rst;
}

Quaternion operator*(const Quaternion& lhs, Float r)
{
    Quaternion rst;
    rst.values[0] = lhs.values[0] * r ;
    rst.values[1] = lhs.values[1] * r ;
    rst.values[2] = lhs.values[2] * r ;
    rst.values[3] = lhs.values[3] * r ;
    return rst;
}

Quaternion operator*(Float r, const Quaternion& rhs)
{
    Quaternion rst;
    rst.values[0] = rhs.values[0] * r ;
    rst.values[1] = rhs.values[1] * r ;
    rst.values[2] = rhs.values[2] * r ;
    rst.values[3] = rhs.values[3] * r ;
    return rst;
}

Quaternion operator/(const Quaternion& lhs, Float r)
{
    Quaternion rst;
    rst.values[0] = lhs.values[0] / r ;
    rst.values[1] = lhs.values[1] / r ;
    rst.values[2] = lhs.values[2] / r ;
    rst.values[3] = lhs.values[3] / r ;
    return rst;
}



Quaternion Quaternion::conj() const {
    return Quaternion ( w(), -x(), -y(), -z() ) ;
}

Vector3D Quaternion::rotVector (const Vector3D& v) const {
    Quaternion q = *this ;
    Quaternion V(0, v.x(), v.y(), v.z());
    Quaternion tmp = q*V*q.conj();
    return Vector3D(tmp.x(), tmp.y(), tmp.z());
}





Quaternion Quaternion::v2q(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D V1 = vec1.getNormalized() ;
    Vector3D V2 = vec2.getNormalized() ;

    Vector3D V3 = V1 + V2 ;
    if ( V3.getLength() > _0_ )
    {
        V3.normalize();
        Vector3D Cross = Vector3D::Cross(V1, V3);
        Float w = Vector3D::Dot(V1,V3);
        Quaternion q(w, Cross.x(), Cross.y(), Cross.z());
        q.normalize();
        return q;
    }
    else
    {
        Vector3D Any1 (1,0,0) ;
        Vector3D Any2 (0,1,0) ;
        Vector3D Cross1 = Vector3D::Cross(V1, Any1);
        Vector3D Cross2 = Vector3D::Cross(V1, Any2);
        Vector3D Cross = Cross1.getLength() > Cross2.getLength() ? Cross1 : Cross2 ;
        Quaternion q(0, Cross.x(), Cross.y(), Cross.z());
        q.normalize();
        return q;
    }
}

Vector3D Quaternion::rotVbyQ (const Vector3D& v, const Quaternion& q)
{
    Quaternion V(0, v.x(), v.y(), v.z());
    Quaternion tmp = q*V*q.conj();
    return Vector3D(tmp.x(), tmp.y(), tmp.z());
}

EulerAngle Quaternion::toEuler (EulerOrder order) const {
    Quaternion q = *this ;

    Float r11, r12, r21, r31, r32, r4;
    if(order == EulerOrder_YXZ)
    {
        r11 = 2 * (q.x() * q.z() + q.w() * q.y());
        r12 = q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z()*q.z();
        r21 = -2 * (q.y()*q.z() - q.w() * q.x());
        r31 = 2*(q.x() * q.y() + q.w()*q.z());
        r32 =q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z()*q.z();
        r4  = q.y() ;
    }
  /////
    else if(order == EulerOrder_ZYX)
    {
        r11 = 2 * (q.x() * q.y() + q.w() * q.z());
        r12 = q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z()*q.z();
        r21 = -2 * (q.x() * q.z() - q.w() * q.y());
        r31 = 2 * (q.y() * q.z() + q.w() * q.x());
        r32 = q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z()*q.z();
        r4  = q.z() ;
    }
    else if(order == EulerOrder_ZXY)
    {
        r11 = -2 * (q.x() * q.y() - q.w() * q.z());
        r12 = q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z()*q.z();
        r21 = 2 * (q.y() * q.z() + q.w() * q.x());
        r31 = -2 * (q.x() * q.z() - q.w() * q.y());
        r32 = q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z()*q.z();
        r4  = q.z() ;
    }
    else if(order == EulerOrder_YZX)
    {
        r11 = -2 * (q.x() * q.z() - q.w() * q.y());
        r12 = q.w()*q.w() + q.x()*q.x() - q.y()*q.y() - q.z()*q.z();
        r21 = 2 * (q.x() * q.y() + q.w() * q.z());
        r31 = -2 * (q.y() * q.z() - q.w() * q.x());
        r32 = q.w()*q.w() - q.x()*q.x() + q.y()*q.y() - q.z()*q.z();
        r4  = q.y() ;
    }
    /////////
    else if(order == EulerOrder_XYZ)
    {
        r11 = -2 * (q.y() * q.z() - q.w() * q.x());
        r12 = q.w() * q.w() - pow(q.x(), 2) - q.y() * q.y() + q.z()*q.z();
        r21 = 2 * (q.x() * q.z() + q.w() * q.y());
        r31 = -2 * (q.x() * q.y() - q.w() * q.z());
        r32 = q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z()*q.z();
        r4  = q.x() ;
    }
    else if(order == EulerOrder_XZY)
    {
        r11 = 2 * (q.y() * q.z() + q.w() * q.x());
        r12 = q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z()*q.z();
        r21 = -2 * (q.x() * q.y() - q.w() * q.z());
        r31 = 2 * (q.x() * q.z() + q.w() * q.y());
        r32 = q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z()*q.z();
        r4  = q.x() ;
    }
    else
    {
        return EulerAngle(0,0,0);
    }
    if ( r21 > 0.9999 )
    {
        r21 = r21 > 1 ? 1 : r21 ;
        EulerAngle E( 2 * atan2f(r4, q.w()), asinf(r21), 0) ;
        if(E.yaw() > 180*D2R)
            E.yaw() -= 360 * D2R;
        else if(E.yaw() < -180 * D2R)
            E.yaw() += 360 * D2R;
        return E;
    }
    else if ( r21 < -0.9999 )
    {
        r21 = r21 < (-1) ? (-1) : r21 ;
        EulerAngle E( 2 * atan2f(r4, q.w()), asinf(r21), 0) ;
        if(E.yaw() > 180*D2R)
            E.yaw() -= 360 * D2R;
        else if(E.yaw() < -180 * D2R)
            E.yaw() += 360 * D2R;
        return E;
    }
    else
    {
        return EulerAngle(atan2f(r11, r12), asinf(r21), atan2f(r31,r32));
    }
}

Matrix Quaternion::toMat() const
{
    Float w,x,y,z;
    w = values[0] ;
    x = values[1] ;
    y = values[2] ;
    z = values[3] ;

    Matrix RotMat(3, 3);
    Float sqw = w * w;
    Float sqx = x * x;
    Float sqy = y * y;
    Float sqz = z * z;

    RotMat(0, 0) = 1 - 2*(sqz + sqy);
    RotMat(1, 1) = 1 - 2*(sqz + sqx);
    RotMat(2, 2) = 1 - 2*(sqy + sqx);
    Float xy = x * y;
    Float zw = z * w;
    RotMat(0, 1) = 2.0 * (xy - zw);
    RotMat(1, 0) = 2.0 * (xy + zw);
    Float xz = x * z;
    Float yw = y * w;
    RotMat(0, 2) = 2.0 * (xz + yw);
    RotMat(2, 0) = 2.0 * (xz - yw);
    Float xw = x * w;
    Float yz = y * z;
    RotMat(1, 2) = 2.0 * (yz - xw);
    RotMat(2, 1) = 2.0 * (yz + xw);

//    RotMat(0, 0) = sqw + sqx - sqy - sqz;
//    RotMat(1, 1) = sqw - sqx + sqy - sqx;
//    RotMat(2, 2) = sqw - sqx - sqy + sqz;
//    Float xy = x * y;
//    Float zw = z * w;
//    RotMat(0, 1) = 2.0 * (xy - zw);
//    RotMat(1, 0) = 2.0 * (xy + zw);
//    Float xz = x * z;
//    Float yw = y * w;
//    RotMat(0, 2) = 2.0 * (xz + yw);
//    RotMat(2, 0) = 2.0 * (xz - yw);
//    Float xw = x * w;
//    Float yz = y * z;
//    RotMat(1, 2) = 2.0 * (yz - xw);
//    RotMat(2, 1) = 2.0 * (yz + xw);
    return RotMat;
}

bool Quaternion::fromMat(const Matrix& M)
{
    if ( M.rows() != 3 || M.cols() != 3 )
        return false ;

   Matrix A = M;
   int idx, i;
   Float U[4];
   Float max = 0;
   Float absm = 0;

   U[0] = 1 + A(0,0) + A(1,1) + A(2,2);
   U[1] = 1 + A(0,0) - A(1,1) - A(2,2);
   U[2] = 1 - A(0,0) + A(1,1) - A(2,2);
   U[3] = 1 - A(0,0) - A(1,1) + A(2,2);
   for(i=0;i<4;i++)
   {
       absm = fabs(U[i]);
       if(absm > max)
       {
           max = absm;
           idx = i;
       }
   }
   switch(idx)
   {
   case 0:
       values[0] = U[0];
       values[1] = A(2, 1) - A(1, 2);
       values[2] = A(0, 2) - A(2, 0);
       values[3] = A(1, 0) - A(0, 1);
       break;
   case 1:
       values[0] = A(2, 1) - A(1, 2);
       values[1] = U[1];
       values[2] = A(0, 1) + A(1, 0);
       values[3] = A(2, 0) + A(0, 2);
       break;
   case 2:
       values[0] = A(0, 2) - A(2, 0);
       values[1] = A(0, 1) + A(1, 0);
       values[2] = U[2];
       values[3] = A(2, 1) + A(1, 2);
       break;
   case 3:
       values[0] = A(1, 0) - A(0, 1);
       values[1] = A(0, 2) + A(2, 0);
       values[2] = A(2, 1) + A(1, 2);
       values[3] = U[3];
       break;
   default:
       break;
   }

   normalize();
   return true ;
}

AxisAngle Quaternion::toAxisAngle() const {

    Float w,x,y,z;
    Quaternion q = *this ;
    q.normalize();
    w = q.w() ;
    if ( w > _1_ )
        w = _1_ ;
    else if ( w < -_1_ )
        w = -_1_ ;

    Float SinTheta2 = 1-w*w ;
    if (SinTheta2 < 0 )
        SinTheta2 = 0 ;
    Float s = sqrt(SinTheta2);
    if ( s > _1_ )
        s = _1_ ;

    if ( s < _0_ )
    {
        s = _0_ ;
        x = y = z = 1 ;
    }
    else
    {
        x = q.x() / s ;
        y = q.y() / s ;
        z = q.z() / s ;
    }
    Float angle = 2 * asin (s) ;
    if ( w < 0 )
        angle = -angle ;
    return AxisAngle (Vector3D(x,y,z).normalize(), angle) ;
}





Float Quaternion::Dot (const Quaternion& q1, const Quaternion& q2)
{
    return  q1.w() * q2.w() +
            q1.x() * q2.x() +
            q1.y() * q2.y() +
            q1.z() * q2.z() ;
}

Quaternion Quaternion::Slerp(const Quaternion& startQ, const Quaternion& endQ, Float endWeight)
{
    Quaternion start = startQ ;
    Quaternion end = endQ ;
    Float dot = Dot (start, end) ;
    if ( dot < 0.0 )
    {
        end = - end ;
        dot = - dot ;
    }

    dot = dot > _1_ ? _1_ : dot ;
    endWeight = endWeight > 1.0 ? 1.0 : endWeight ;
    endWeight = endWeight < 0.0 ? 0.0 : endWeight ;

    Float theta = acos(dot) ;
    Float sinTheta = sin(theta) ;
    Float rStart = 1.0 - endWeight , rEnd = endWeight ;
    if ( sinTheta > _0_ )
    {
        rEnd = sin(endWeight*theta) / sinTheta ;
        rStart = sin(theta - endWeight*theta) / sinTheta ;
    }
    Quaternion rst = rStart * start + rEnd * end ;
    //rst.normalize();
    return rst ;
}


EulerAngle::EulerAngle()
{
    Yaw  =  Pitch  =  Roll  =  0 ;
}

EulerAngle::EulerAngle(Float yaw, Float pitch, Float roll)
{
    Yaw = yaw ;  Pitch = pitch ; Roll = roll ;
}

Float& EulerAngle::yaw()
{
    return Yaw ;
}

Float& EulerAngle::pitch()
{
    return Pitch ;
}

Float& EulerAngle::roll()
{
    return Roll ;
}
Float EulerAngle::yaw() const
{
    return Yaw ;
}
Float EulerAngle::pitch() const
{
    return Pitch ;
}
Float EulerAngle::roll() const
{
    return Roll ;
}


EulerAngle EulerAngle::DegreeToRad () const
{
    return EulerAngle( yaw()*D2R, pitch()*D2R, roll()*D2R ) ;
}

EulerAngle EulerAngle::RadToDegree () const
{
    return EulerAngle( yaw()*R2D, pitch()*R2D, roll()*R2D ) ;
}



Quaternion EulerAngle::toQuat(EulerOrder order) const{
    EulerAngle EA = *this ;

    Float ang[3] = {EA.yaw(), EA.pitch(), EA.roll()};
    Float cang[3] = {cosf(0.5 * ang[0]), cosf(0.5 * ang[1]), cosf(0.5 * ang[2])};
    Float sang[3] = {sinf(0.5 * ang[0]), sinf(0.5 * ang[1]), sinf(0.5 * ang[2])};
    Float w, x, y, z;
    if(order == EulerOrder_YXZ)
    {
        w = cang[0] * cang[1] * cang[2] + sang[0] * sang[1] * sang[2];
        x = cang[0] * sang[1] * cang[2] + sang[0] * cang[1] * sang[2];
        y = sang[0] * cang[1] * cang[2] - cang[0] * sang[1] * sang[2];
        z = cang[0] * cang[1] * sang[2] - sang[0] * sang[1] * cang[2];
    }
    else if(order == EulerOrder_ZYX)
    {
        w = cang[0] * cang[1] * cang[2] + sang[0] * sang[1] * sang[2];
        x = cang[0] * cang[1] * sang[2] - sang[0] * sang[1] * cang[2];
        y = cang[0] * sang[1] * cang[2] + sang[0] * cang[1] * sang[2];
        z = sang[0] * cang[1] * cang[2] - cang[0] * sang[1] * sang[2];
    }
    else if(order == EulerOrder_ZXY)
    {
        w = cang[0] * cang[1] * cang[2] - sang[0] * sang[1] * sang[2];
        x = cang[0] * sang[1] * cang[2] - sang[0] * cang[1] * sang[2];
        y = cang[0] * cang[1] * sang[2] + sang[0] * sang[1] * cang[2];
        z = cang[0] * sang[1] * sang[2] + sang[0] * cang[1] * cang[2];
    }
    else if(order == EulerOrder_YZX)
    {
        w = cang[0] * cang[1] * cang[2] - sang[0] * sang[1] * sang[2];
        x = cang[0] * cang[1] * sang[2] + sang[0] * sang[1] * cang[2];
        y = cang[0] * sang[1] * sang[2] + sang[0] * cang[1] * cang[2];
        z = cang[0] * sang[1] * cang[2] - sang[0] * cang[1] * sang[2];
    }
    else if(order == EulerOrder_XYZ)
    {
        w = cang[0] * cang[1] * cang[2] - sang[0] * sang[1] * sang[2];
        x = cang[0] * sang[1] * sang[2] + sang[0] * cang[1] * cang[2];
        y = cang[0] * sang[1] * cang[2] - sang[0] * cang[1] * sang[2];
        z = cang[0] * cang[1] * sang[2] + sang[0] * sang[1] * cang[2];
    }
    else if(order == EulerOrder_XZY)
    {
        w = cang[0] * cang[1] * cang[2] + sang[0] * sang[1] * sang[2];
        x = sang[0] * cang[1] * cang[2] - cang[0] * sang[1] * sang[2];
        y = cang[0] * cang[1] * sang[2] - sang[0] * sang[1] * cang[2];
        z = cang[0] * sang[1] * cang[2] + sang[0] * cang[1] * sang[2];
    }
    return Quaternion(w, x, y, z);
}



AxisAngle::AxisAngle()
{
    Axis = Vector3D (0,0,1) ;
    Angle = 0 ;
}

AxisAngle::AxisAngle ( const Vector3D& axis, Float angle )
{
    Axis = axis ; Angle = angle ;
}

Vector3D& AxisAngle::axis()
{
    return Axis ;
}

Float& AxisAngle::angle()
{
    return Angle ;
}

const Vector3D& AxisAngle::axis() const
{
    return Axis ;
}
Float AxisAngle::angle() const
{
    return Angle ;
}


AxisAngle AxisAngle::DegreeToRad () const
{
    return AxisAngle( Axis, Angle*D2R ) ;
}

AxisAngle AxisAngle::RadToDegree () const
{
    return AxisAngle( Axis, Angle*R2D ) ;
}

Quaternion AxisAngle::toQuat () const {
    // Float RadAng = Angle*D2R;
    Float RadAng = Angle ;
    Float qw = cos(RadAng/2);
    Float qx = Axis.x() * sin(RadAng/2);
    Float qy = Axis.y() * sin(RadAng/2);
    Float qz = Axis.z() * sin(RadAng/2);
    Quaternion Q = Quaternion(qw,qx,qy,qz);
    Q.normalize();
    return Q ;
}

