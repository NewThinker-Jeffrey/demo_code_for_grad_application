#ifndef __MathToolsCpp__Matrix__
#define __MathToolsCpp__Matrix__
#include <cstddef>
#include "MathToolsCpp/cppFloat.h"

/***
 * Some features :
 *
 *      1.  Free matrix
 *
 *          "Matrix A(3,4);" creates a free matrix whose data is allocated in the
 *          heap.
 *
 *      2.  Lazy copy
 *
 *          "Matrix B=A;" creates another matrix object "B" sharing data with "A",
 *          supposing "A" is already a FREE matrix object. This operation is very
 *          efficient since no memory copy needs to be performed. The data won't
 *          be cloned until one of the matrices (sharing the same data) be modified.
 *
 *      3.  Element-wise operation
 *          It's supported to access every entry of a matrix, e.g. "A(3,4)" represents
 *          the entry at row 3 and col 4 (0-based indexes).
 *
 *      4.  Fixbinded matrix
 *
 *          It's supported to create a matrix object from an existing float array.
 *          Such a matrix is said to be "Fixbinded". Operations on the matrix will
 *          immediately affect the data in the array, and vice versa.
 *
 *          A fix-binded matrix is not a free matrix.
 *
 *      5.  Parasitic matrix (sub-matrix and transpose)
 *
 *          "Matrix B = A.subMat(startRow, startCol, rows, cols);" creates a
 *          sub-matrix "B" of "A", supposing "A" is already a matrix object.
 *          Operations on "B" will affect "A".
 *
 *          "Matrix B = A.transpose()" creates a matrix "B" that acts like the
 *          transpose of "A" and shares memory with "A". Operations on "B" will
 *          affect "A".
 *
 *          A parasitic matrix is not a free matrix.
 *
 * */

#define ENABLE_MATRIX_IO
#define DEBUG_MATRIX
//#define SAFE_MATRIX  // not supported yet .

START_JK_Math_NAMESPACE
class Matrix ;
class MatrixPrivate ;
class Vector3D ;
END_JK_Math_NAMESPACE


JKMath::Matrix operator+ (const JKMath::Matrix &A, const JKMath::Matrix &B);
JKMath::Matrix operator+ (const JKMath::Matrix &A, Float b);
JKMath::Matrix operator- (const JKMath::Matrix &A, const JKMath::Matrix &B);
JKMath::Matrix operator- (const JKMath::Matrix &A, Float b);
JKMath::Matrix operator- (const JKMath::Matrix &A);
JKMath::Matrix operator* (const JKMath::Matrix &A, const JKMath::Matrix &B);
JKMath::Matrix operator* (const JKMath::Matrix &A, const JKMath::Vector3D &B);
JKMath::Matrix operator* (const JKMath::Vector3D &A, const JKMath::Matrix &B);
JKMath::Matrix operator* (Float b, const JKMath::Matrix &A);
JKMath::Matrix operator* (const JKMath::Matrix &A, Float b);
JKMath::Matrix operator/ (const JKMath::Matrix &A, Float b);


START_JK_Math_NAMESPACE

class Matrix {

public:

    // constructors
    Matrix() ;
    Matrix(unsigned rows, unsigned cols, Float*data=NULL, int rowStep=0, int colStep=0);
    Matrix(const Matrix &A);
    Matrix &operator=(const Matrix &A);
    Matrix (double s) ;
    operator double() const ;
    Matrix &operator=(double s) ;
    ~Matrix();

    // information
    unsigned rows() const ;
    unsigned cols() const ;
    bool isFixBinded() const ;
    bool isParasitic() const ;
    bool isReadOnly() const ;
    bool isEmpty() const ;


    // some convenient interfaces .
    enum DataOption { CloneData, ShareData } ;
    inline void reset();
    void fill(Float v=0);
    Float sum() const ;
    Float squareSum() const ;


    Matrix clone() const;
    Matrix transpose(DataOption op=ShareData) const;
    Matrix subMat(unsigned i, unsigned j,
                  unsigned rows, unsigned cols,
                  DataOption op=ShareData, bool transposeFirst=false) const;
    Matrix row(int i, DataOption op=ShareData) const ;
    Matrix col(int j, DataOption op=ShareData) const ;


    Matrix transpose(DataOption op=ShareData) ;
    Matrix subMat(unsigned i, unsigned j,
                  unsigned rows, unsigned cols,
                  DataOption op=ShareData, bool transposeFirst=false) ;
    Matrix row(int i, DataOption op=ShareData) ;
    Matrix col(int j, DataOption op=ShareData) ;


    // arithmatic functions
    static Matrix& Add ( const Matrix& A, const Matrix& B, Matrix& C) ;
    static Matrix& Add ( const Matrix& A, Float b, Matrix& C) ;
    static Matrix& Sub ( const Matrix& A, const Matrix& B, Matrix& C) ;
    static Matrix& Sub ( const Matrix& A, Float b, Matrix& C) ;
    static Matrix& Mult ( const Matrix& A, const Matrix& B, Matrix& C) ;
    static Matrix& Mult ( const Matrix& A, Float b, Matrix& C) ;
    static Matrix& Div ( const Matrix& A, Float b, Matrix& C ) ;    
    static Float Dot ( const Matrix& A, const Matrix& B) ;




    // arithmatic operations .
    friend Matrix (::operator+) (const Matrix &A, const Matrix &B);
    friend Matrix (::operator+) (const Matrix &A, Float b);
    friend Matrix (::operator-) (const Matrix &A, const Matrix &B);
    friend Matrix (::operator-) (const Matrix &A, Float b);
    friend Matrix (::operator-) (const Matrix &A);
    friend Matrix (::operator*) (const Matrix &A, const Matrix &B);
    friend Matrix (::operator*) (const Matrix &A, const JKMath::Vector3D &B);
    friend Matrix (::operator*) (const JKMath::Vector3D &A, const Matrix &B);
    friend Matrix (::operator*) (Float b, const Matrix &A);
    friend Matrix (::operator*) (const Matrix &A, Float b);
    friend Matrix (::operator/) (const Matrix &A, Float b);
    Matrix& operator+=(const Matrix &B) ;
    Matrix& operator+=(Float b) ;
    Matrix& operator-=(const Matrix &B) ;
    Matrix& operator-=(Float b) ;
    Matrix& operator*=(Float b) ;
    Matrix& operator/=(Float b) ;


    Float det() const ;
    Float trace() const ;
    Matrix getInverse() const;
    bool GaussElimination ( bool backEnable=false, Matrix* follow=NULL, int* sign=NULL) ;
    static Matrix I(int n) ;
    static Matrix solve(const Matrix& A, const Matrix& Y) ;


    // for single element access
    class ElementRO ;
    class Element ;
    Element operator()(unsigned i, unsigned j) ;
    ElementRO operator()(unsigned i, unsigned j) const ;


    // IO functions .
#ifdef ENABLE_MATRIX_IO
    int printMat() const;
    int save2file(const char *filename) const;
    bool LoadMatFile(const char *filename) ;
    int save2Str (char* str, int bufLen) const ;
    bool loadMat (const char* str, int strLen=0) ;
#endif



public:
    // medium adaptor-classes for single element access .
    class ElementRO {
        friend class Matrix ;
    public:
        ElementRO (const ElementRO &other);
        operator Float() const ;

    private:
        unsigned i, j ;
        const Matrix* m;
        // a private constructor to prevent an unexpected misuse of
        // this inner-class at outer codes .
        ElementRO (const Matrix*, unsigned, unsigned) ;
    };

    class Element {
        friend class Matrix ;

    public:
        Element(const Element&other);
        Element& operator =(const Element&other) ;
        Element& operator =(const ElementRO &other) ;

        Element& operator =(Float f) ;
        operator Float() const ;

    private:
        unsigned i, j ;
        Matrix* m;
        // a private constructor to prevent an unexpected misuse of
        // this inner-class at outer codes .
        Element(Matrix*, unsigned, unsigned) ;
    };


private:
    MatrixPrivate* Priv;
    inline const MatrixPrivate* d() const ;
    inline MatrixPrivate*& d() ;
    Matrix(const MatrixPrivate*) ;
    friend class ElementRO ;
    friend class Element ;
    friend class MatrixPrivate ;
} ;


END_JK_Math_NAMESPACE


#ifdef DEBUG_MATRIX
extern int MatrixInstances     ;
extern int MatrixPrivInstances ;
extern int MatrixDataInstances ;
#include <cstdio>
static inline void MATRIX_DEBUG_INFO(const char* head=NULL)
{
    if (head)       printf("%s\n", head);
    printf("MatrixInstances    : %d\n",MatrixInstances    );
    printf("MatrixPrivInstances: %d\n",MatrixPrivInstances);
    printf("MatrixDataInstances: %d\n",MatrixDataInstances);
}
#endif



#endif // MATRIX_H
