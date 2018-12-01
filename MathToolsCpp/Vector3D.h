#ifndef __MathtoolsCpp__Vector3D__
#define __MathtoolsCpp__Vector3D__
 
#include "MathToolsCpp/cppFloat.h"

START_JK_Math_NAMESPACE
class Vector3D ;
class Matrix ;
END_JK_Math_NAMESPACE


JKMath::Vector3D operator+ (const JKMath::Vector3D& lhs, const JKMath::Vector3D& rhs) ;
JKMath::Vector3D operator- (const JKMath::Vector3D& lhs) ;
JKMath::Vector3D operator- (const JKMath::Vector3D& lhs, const JKMath::Vector3D& rhs) ;
JKMath::Vector3D operator* (const JKMath::Vector3D& lhs, const Float& num) ;
JKMath::Vector3D operator* (const Float& num, const JKMath::Vector3D& rhs ) ;
JKMath::Vector3D operator/ (const JKMath::Vector3D& lhs, const Float& num) ;


START_JK_Math_NAMESPACE

class Vector3D
{
    
public:
    Vector3D() ;
    Vector3D(Float x, Float y, Float z) ;
    Vector3D(const Float V[]) ;
    Vector3D(const Vector3D& other ) ;
    Vector3D(const Matrix& mat ) ;
    operator Matrix() const ;
    Vector3D & operator=(const Matrix& mat) ;
    Vector3D & operator=(const Vector3D& other) ;
    Vector3D & operator+=(const Vector3D& rhs) ;
    Vector3D & operator-=(const Vector3D& rhs) ;
    Vector3D & operator*=(Float rhs) ;
    Vector3D & operator/=(Float rhs) ;

    Float& x() ;
    Float& y() ;
    Float& z() ;
    Float x() const ;
    Float y() const ;
    Float z() const ;


    void setValues(Float X, Float Y, Float Z) ;
    Float getLength() const ;
    Float getNorm() const ;
    const Vector3D& SetLength(Float length) ;
    const Vector3D& normalize () ;
    Vector3D getNormalized () const ;

    friend Vector3D (::operator+) (const Vector3D& lhs, const Vector3D& rhs) ;
    friend Vector3D (::operator-) (const Vector3D& lhs) ;
    friend Vector3D (::operator-) (const Vector3D& lhs, const Vector3D& rhs) ;
    friend Vector3D (::operator*) (const Vector3D& lhs, const Float& num) ;
    friend Vector3D (::operator*) (const Float& num, const Vector3D& rhs ) ;
    friend Vector3D (::operator/) (const Vector3D& lhs, const Float& num) ;

public :
    static Float Dot ( const Vector3D& vec1, const Vector3D& vec2 ) ;
    static Vector3D Cross(const Vector3D& vec1, const Vector3D& vec2) ;

private:
    Float values[3];
};

END_JK_Math_NAMESPACE
#endif
