
#include <cmath>
#include "MathToolsCpp/Vector3D.h"
#include "MathToolsCpp/Matrix.h"
IMPORT_JK_Math


Vector3D::Vector3D(){
    values[0] = values[1] = values[2] = 0 ;
}

Vector3D::Vector3D(Float x, Float y, Float z) {
    values[0] = x ;
    values[1] = y ;
    values[2] = z ;
}

Vector3D::Vector3D(const Float V[]){
    values[0] = V[0] ;
    values[1] = V[1] ;
    values[2] = V[2] ;
}

Vector3D::Vector3D(const Vector3D& other ) {
    values[0] = other.values[0] ;
    values[1] = other.values[1] ;
    values[2] = other.values[2] ;
}

Vector3D::Vector3D(const Matrix& mat ) {
    if (mat.rows()==3 && mat.cols()==1)
    {
        values[0] = mat(0,0) ;
        values[1] = mat(1,0) ;
        values[2] = mat(2,0) ;
    }
    else
    {
        values[0] = values[1] = values[2] = 0 ;
    }
}


Float& Vector3D::x()
{
    return values[0] ;
}

Float& Vector3D::y()
{
    return values[1] ;
}

Float& Vector3D::z()
{
    return values[2] ;
}

Float Vector3D::x() const
{
    return values[0] ;
}
Float Vector3D::y() const
{
    return values[1] ;
}
Float Vector3D::z() const
{
    return values[2] ;
}


void Vector3D::setValues(Float X, Float Y, Float Z) {
    values[0]=X;
    values[1]=Y;
    values[2]=Z;
}

Float Vector3D::getLength() const {
    return sqrt ( x() * x() + y() * y() + z() * z() ) ;
}

Float Vector3D::getNorm() const {
    return x() * x() + y() * y() + z() * z() ;
}

const Vector3D& Vector3D::SetLength(Float length) {
    Float L = this->getLength() ;
    if ( L > 0 ) {
        *this = *this * (length / L) ;
    }
    return *this ;
}

const Vector3D& Vector3D::normalize () {
    SetLength(1) ;
    return *this;
}

Vector3D Vector3D::getNormalized () const {
    Vector3D V = *this ;
    return V.normalize() ;
}


Vector3D operator+(const Vector3D& lhs, const Vector3D& rhs) {
    return Vector3D( lhs.values[0] + rhs.values[0], lhs.values[1] + rhs.values[1], lhs.values[2] + rhs.values[2]);
}

Vector3D operator-(const Vector3D& lhs) {
    return Vector3D(-lhs.values[0], -lhs.values[1], -lhs.values[2]);
}

Vector3D operator-(const Vector3D& lhs, const Vector3D& rhs) {
    return Vector3D( lhs.values[0] - rhs.values[0], lhs.values[1] - rhs.values[1], lhs.values[2] - rhs.values[2]);
}

Vector3D operator*(const Vector3D& lhs, const Float& num) {
     return Vector3D(lhs.values[0] * num, lhs.values[1] * num, lhs.values[2] * num);
}

Vector3D operator*(const Float& num, const Vector3D& rhs ) {
     return Vector3D(rhs.values[0] * num, rhs.values[1] * num, rhs.values[2] * num);
}

Vector3D operator/(const Vector3D& lhs, const Float& num) {
    return Vector3D(lhs.values[0] / num, lhs.values[1] / num, lhs.values[2] / num);
}

Vector3D::operator Matrix() const
{
    Matrix mat(3,1) ;
    mat(0,0) = values[0] ;
    mat(1,0) = values[1] ;
    mat(2,0) = values[2] ;
    return mat ;
}

Vector3D& Vector3D::operator=(const Matrix& mat)
{
    if (mat.rows()==3 && mat.cols()==1)
    {
        values[0] = mat(0,0) ;
        values[1] = mat(1,0) ;
        values[2] = mat(2,0) ;
    }
    return *this ;
}


Vector3D & Vector3D::operator=(const Vector3D& rhs) {
    values[0] = rhs.values[0] ;
    values[1] = rhs.values[1] ;
    values[2] = rhs.values[2] ;
    return *this ;
}

Vector3D & Vector3D::operator+=(const Vector3D& rhs) {
    *this = *this + rhs ;
    return *this ;
}

Vector3D & Vector3D::operator*=(Float rhs)
{
    *this = *this * rhs ;
    return *this ;
}

Vector3D & Vector3D::operator/=(Float rhs)
{
    *this = *this / rhs ;
    return *this ;
}

Vector3D & Vector3D::operator-=(const Vector3D& rhs) {
    *this = *this - rhs ;
    return *this ;
}



Float Vector3D::Dot ( const Vector3D& vec1, const Vector3D& vec2 ){
    return vec1.x() * vec2.x() +
            vec1.y() * vec2.y() +
            vec1.z() * vec2.z();
}

Vector3D Vector3D::Cross(const Vector3D& vec1, const Vector3D& vec2){
    return Vector3D(vec1.y() * vec2.z() - vec1.z() * vec2.y(),
                   vec1.z() * vec2.x() - vec1.x() * vec2.z(),
                   vec1.x() * vec2.y() - vec1.y() * vec2.x());
}

