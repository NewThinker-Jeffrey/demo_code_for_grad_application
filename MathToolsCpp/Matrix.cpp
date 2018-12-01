#include <cmath>
#include <stdint.h>
#include "float.h"
#include "MathToolsCpp/Matrix.h"
#include "MathToolsCpp/Vector3D.h"


#ifdef ENABLE_MATRIX_IO
#include <cstdio>
#include <cstdlib>
#include <cstring>
#endif


#ifdef StepVR_Algorithm_Pro
#define snprintf _snprintf
#endif



#ifdef DEBUG_MATRIX

    int MatrixInstances = 0 ;
    int MatrixPrivInstances = 0 ;
    int MatrixDataInstances = 0 ;
    #define INC_MATRIX (MatrixInstances++)
    #define DEC_MATRIX (MatrixInstances--)
    #define INC_MATRIX_PRIV (MatrixPrivInstances++)
    #define DEC_MATRIX_PRIV (MatrixPrivInstances--)
    #define INC_MATRIX_DATA (MatrixDataInstances++)
    #define DEC_MATRIX_DATA (MatrixDataInstances--)
    #define DEBUG_ERR(str) printf((str))
    #define DEBUG_INFO(str) printf((str))
    //#define DEBUG_INFO(str)

#else
    #define INC_MATRIX
    #define DEC_MATRIX
    #define INC_MATRIX_PRIV
    #define DEC_MATRIX_PRIV
    #define INC_MATRIX_DATA
    #define DEC_MATRIX_DATA
    #define DEBUG_ERR(str)
    #define DEBUG_INFO(str)
#endif


START_JK_Math_NAMESPACE
enum MatrixAttr {
    Temporary           = (0x1 << 0 ),
    Parasitic           = (0x1 << 1 ),
    Transposed          = (0x1 << 2 ),
    IndependentData     = (0x1 << 3 ),
    ReadOnly            = (0x1 << 4 )
};


class MatrixPrivate {

    friend class Matrix ;
    friend class Matrix::Element ;
    MatrixPrivate () ;
    MatrixPrivate ( unsigned rows, unsigned cols, Float*data=NULL, int rowStep=0, int colStep=0 );
    MatrixPrivate ( unsigned _rows, unsigned _cols, Matrix*_master, unsigned _i, unsigned _j, bool transposed=false) ;
    ~MatrixPrivate();

    // the methods 'share()' 'deShare()' 'deShare' actually should NOT be declared
    // as a 'const' method , they modified something in fact indeed, but
    // consider that the modified is invisible to the outer layer and they keep the
    // core data unchanged, we can trick the outer layer code without feeling guilty .
    // what's more , the "core data" here means : the dimensions(rows/cols) and the
    // elements' values of a Matrix .
    MatrixPrivate* share() const ;
    void deShare() const ;
    void clearShare() const ;
    bool isShared() const ;

    MatrixPrivate* checkWrite() ;
    MatrixPrivate& getData() ;
    const MatrixPrivate& getData() const ;


    inline Float& operator() (int i, int j) ;
    inline Float operator() (int i, int j) const;
    inline Float* data() const ;

    // attributes
    inline void setParasitic(bool) ;
    inline bool isParasitic() const ;
    inline void setIndependent(bool) ;
    inline bool isIndependent() const ;
    inline void setTransposed(bool) ;
    inline bool isTransposed() const ;
    inline void setTemporary(bool) ;
    inline bool isTemporary() const ;
    inline void setReadOnly(bool) ;
    inline bool isReadOnly() const ;
    inline void setDefaultAttr() ;


    // copyFrom(other) will copy all the elements of 'other' to 'this' , and set the
    // dimensions(rows/cols) of 'this' the same with 'other' .
    void copyFrom (const MatrixPrivate *other) ;

    union {
        // for original matrix
        struct {
            int rowStep, colStep ;
            Float* Data ;
        } ;

        // for parasitic matrix
        struct {
            unsigned i, j ;
            Matrix* master ;
            MatrixPrivate* quickData ;
        } ;
    } ;

    unsigned rows, cols;
    unsigned attr ;
    unsigned shared ;    // when shared, it suggests a readonly permmision .
} ;
END_JK_Math_NAMESPACE


IMPORT_JK_Math

MatrixPrivate::MatrixPrivate ()
{
    INC_MATRIX_PRIV;

    setDefaultAttr() ;
    setTemporary(true);
}

MatrixPrivate::MatrixPrivate ( unsigned _rows, unsigned _cols, Float*data, int _rowStep, int _colStep)
{
    INC_MATRIX_PRIV;

    setDefaultAttr() ;
    if ( data )
    {
        setIndependent(false);
        Data = data ;
        if ( _rowStep == 0 )
            _rowStep = _cols ;
        if ( _colStep == 0 )
            _colStep = 1 ;
    }

    else
    {
        setIndependent(true);
        INC_MATRIX_DATA ;
        Data = new Float[_rows*_cols] ;
        _rowStep = _cols ;
        _colStep = 1 ;
    }

    rows = _rows ;
    cols = _cols ;
    rowStep = _rowStep ;
    colStep = _colStep ;

    clearShare();
    setTemporary(false);
}


MatrixPrivate::MatrixPrivate ( unsigned _rows, unsigned _cols,
                               Matrix*_master, unsigned _i, unsigned _j, bool transposed)
{
    INC_MATRIX_PRIV;

    setDefaultAttr() ;
    rows = _rows ;
    cols = _cols ;
    i = _i ;
    j = _j ;
    master = _master ;
    quickData = new MatrixPrivate ;
    setParasitic(true);
    setTransposed(transposed);
    clearShare();
    setTemporary(false);
}


MatrixPrivate::~MatrixPrivate ()
{
    if ( ! isTemporary() && ! isParasitic() && isIndependent() && Data  )
    {
        delete[] Data ;
        DEC_MATRIX_DATA;
    }
    else if ( ! isTemporary() && isParasitic() && quickData )
    {
        delete quickData ;
    }

    DEC_MATRIX_PRIV ;
}


void MatrixPrivate::setParasitic(bool parasitic)
{
    if (parasitic)
        attr |= Parasitic ;
    else
        attr &= (~Parasitic) ;
}

bool MatrixPrivate::isParasitic() const
{
    if ( attr & Parasitic )
        return true ;
    else
        return false ;
}

void MatrixPrivate::setIndependent(bool independent)
{
    if (independent)
        attr |= IndependentData ;
    else
        attr &= (~IndependentData) ;
}


bool MatrixPrivate::isIndependent() const
{
    if (attr & IndependentData)
        return true ;
    else
        return false ;
}

void MatrixPrivate::setTransposed(bool transposed)
{
    if (transposed)
        attr |= Transposed ;
    else
        attr &= (~Transposed) ;
}


bool MatrixPrivate::isTransposed() const
{
    if (attr & Transposed)
        return true ;
    else
        return false ;
}

void MatrixPrivate::setTemporary(bool temporary)
{
    if (temporary)
        attr |= Temporary ;
    else
        attr &= (~Temporary) ;
}

bool MatrixPrivate::isTemporary() const
{
    if (attr & Temporary)
        return true ;
    else
        return false ;
}

inline void MatrixPrivate::setReadOnly(bool readOnly)
{
    if (readOnly)
        attr |= ReadOnly ;
    else
        attr &= (~ReadOnly) ;
}

inline bool MatrixPrivate::isReadOnly() const
{
    if (attr & ReadOnly)
        return true ;
    else
        return false ;
}
void MatrixPrivate::setDefaultAttr()
{
    attr = 0 ;
}



MatrixPrivate* MatrixPrivate::share() const
{
    MatrixPrivate* This = (MatrixPrivate*)this;
    This->shared ++ ;
    return This ;
}
void MatrixPrivate::deShare() const
{
    MatrixPrivate* This = (MatrixPrivate*)this;
    This->shared -- ;
}
void MatrixPrivate::clearShare() const
{
    MatrixPrivate* This = (MatrixPrivate*)this;
    This->shared = 0 ;
}
bool MatrixPrivate::isShared() const
{
    if (shared)
        return true ;
    else
        return false ;
}


MatrixPrivate* MatrixPrivate::checkWrite()
{
    if ( isTemporary() || isReadOnly() )
    {
        return this ;
    }

    else if ( isParasitic() )
    {
        master->d() = master->d()->checkWrite() ;
        return this ;
    }

    else if ( ! isIndependent() || ! isShared()  )
    {
        return this ;
    }

    else // if ( ! isTemporary() && ! isParasitic() && isIndependent() && isShared()  )
    {
        MatrixPrivate* d = new MatrixPrivate(this->rows, this->cols) ;
        d->copyFrom(this);
        this->deShare() ;
        DEBUG_INFO("MatrixPrivate::checkWrite:  MatrixPrivate split \n");
        return d ;
    }
}


const MatrixPrivate& MatrixPrivate::getData() const
{
    MatrixPrivate* This = (MatrixPrivate*)this ;
    const MatrixPrivate* d = &(This->getData()) ;
    return *d ;
}

MatrixPrivate& MatrixPrivate::getData()
{    
    if ( isTemporary() || ! isParasitic() )
    {
        return *this ;
    }

    MatrixPrivate* d=quickData ;
    d->rows = rows ;
    d->cols = cols ;

    // get the top master
    int iSrc = 0, jSrc = 0 ;
    bool transposed = false ;
    const MatrixPrivate* cur = this ;

    while ( cur->isParasitic() )
    {
        int i=iSrc, j=jSrc;
        if (cur->isTransposed())
        {
            transposed = !transposed ;
            iSrc = j + cur->j ;
            jSrc = i + cur->i ;
        }
        else
        {
            iSrc = i + cur->i ;
            jSrc = j + cur->j ;
        }
        cur = cur->master->d() ;
    }

    if ( transposed )
    {
        d->colStep = cur->rowStep ;
        d->rowStep = cur->colStep ;
    }
    else
    {
        d->colStep = cur->colStep ;
        d->rowStep = cur->rowStep ;
    }
    d->Data = cur->Data + iSrc * cur->rowStep + jSrc * cur->colStep ;
    d->setReadOnly( isReadOnly() );
    return *d ;
}


Float& MatrixPrivate::operator() (int i, int j)
{
#ifdef DEBUG_MATRIX
    if ( ! isTemporary() && isParasitic() )
    {
        DEBUG_ERR("MatrixPrivate::operator(): only temporary or original(non-Parasitic) objects can call operator()\n") ;
    }
    if ( isReadOnly() )
    {
        DEBUG_ERR("MatrixPrivate::operator(): ReadOnly Matrix should not call operator()\n") ;
    }
#endif

    return data()[rowStep*i + colStep*j] ;
}

Float MatrixPrivate::operator() (int i, int j) const
{
#ifdef DEBUG_MATRIX
    if ( ! isTemporary() && isParasitic() )
    {
        DEBUG_ERR("MatrixPrivate::operator(): only temporary or original(non-Parasitic) objects can call operator()\n") ;
    }
#endif
    return data()[rowStep*i + colStep*j] ;
}


Float* MatrixPrivate::data() const
{
    return Data ;
}


void MatrixPrivate::copyFrom (const MatrixPrivate *other)
{
    if ( NULL == other )
    {
        DEBUG_ERR("MatrixPrivate::copyFrom: couldn't copy from a NULL matrix\n");
        return ;
    }

    if ( rows != other->rows || cols != other->cols )
    {
        DEBUG_ERR("MatrixPrivate::copyFrom: couldn't copy from a matrix with different dimensions\n");
        return ;
    }
    MatrixPrivate & dThis = this->getData() ;
    const MatrixPrivate & dOther = other->getData() ;

    for ( int i=0; i<rows; i++ )
    {
        for ( int j=0; j<cols; j++ )
        {
            dThis(i,j) = dOther(i,j) ;
        }
    }
}


Matrix::ElementRO::ElementRO (const Matrix*m, unsigned i, unsigned j)
{
    this->m = m ;
    this->i = i ;
    this->j = j ;
}
Matrix::ElementRO::ElementRO (const ElementRO &other)
{
    i = other.i ;
    j = other.j ;
    m = other.m ;
}
Matrix::ElementRO::operator Float() const
{
    const MatrixPrivate & d = *(m->d()) ;
    if ( i >= m->rows() || j >= m->cols() )
        return NAN ;

    if ( d.isParasitic() )
    {
        int iSrc, jSrc ;
        if ( d.isTransposed() )
        {
            jSrc = i + d.i ;
            iSrc = j + d.j ;
        }
        else
        {
            iSrc = i + d.i ;
            jSrc = j + d.j ;
        }
        return (Float)( (*(d.master))(iSrc, jSrc) ) ;
    }

    // if not Parasitic
    return d(i,j) ;
}


Matrix::Element::Element(Matrix*m, unsigned i, unsigned j)
{
    this->m = m ;
    this->i = i ;
    this->j = j ;
}
Matrix::Element::Element(const Element&other)
{
    i = other.i ;
    j = other.j ;
    m = other.m ;
}

Matrix::Element& Matrix::Element::operator =(const Element&other)
{
    return (*this) = (Float)(other) ;
}

Matrix::Element& Matrix::Element::operator =(const ElementRO &other)
{
    return (*this) = (Float)(other) ;
}

Matrix::Element& Matrix::Element::operator =(Float f)
{        
    if ( m->isReadOnly() )
    {
        DEBUG_ERR("Matrix::Element::operator= : Should NOT assign a ReadOnly Matrix\n");
        return *this ;
    }

    if ( i<m->rows() && j<m->cols() )
    {        
        MatrixPrivate & d = *(m->d()) ;
        if ( d.isParasitic() )
        {
            int iDst, jDst ;
            if ( d.isTransposed() )
            {
                iDst = j + d.j ;
                jDst = i + d.i ;
            }
            else
            {
                iDst = i + d.i ;
                jDst = j + d.j ;
            }
            (*(d.master))(iDst, jDst) = f ;
        }

        else
        {
            m->d() = m->d()->checkWrite() ;
            MatrixPrivate & d = *(m->d()) ;
            d(i,j) = f;
        }
    }

    return *this ;
}


Matrix::Element::operator Float() const
{
    const MatrixPrivate & d = *(m->d()) ;
    if ( i >= m->rows() || j >= m->cols() )
        return NAN ;

    if ( d.isParasitic() )
    {
        int iSrc, jSrc ;
        if ( d.isTransposed() )
        {
            jSrc = i + d.i ;
            iSrc = j + d.j ;
        }
        else
        {
            iSrc = i + d.i ;
            jSrc = j + d.j ;
        }
        return (Float)( (*(d.master))(iSrc, jSrc) ) ;
    }

    // if not Parasitic
    return d(i,j) ;
}



Matrix::Matrix()
{
    INC_MATRIX ;
    d() = NULL ;
}

Matrix::Matrix(const MatrixPrivate*_d)
{
    INC_MATRIX ;
    d() = (MatrixPrivate*)_d ;
}

Matrix::Matrix(const Matrix &A)
{
    INC_MATRIX ;
    if ( A.d() )
        d() = A.d()->share() ;
    else
        d() = NULL ;
}

Matrix::Matrix(double s)
{
    INC_MATRIX ;
    d() = new MatrixPrivate(1, 1) ;
    (*this)(0,0) = s ;
}


Matrix::Matrix(unsigned rows, unsigned cols, Float*data, int rowStep, int colStep)
{
    INC_MATRIX ;
    if ( rows * cols > 0 )
        d() = new MatrixPrivate(rows, cols, data, rowStep, colStep) ;
    else
        d() = NULL ;
}

Matrix::~Matrix()
{
    reset() ;
    DEC_MATRIX ;
}


unsigned Matrix::cols() const
{
    if (d())
        return d()->cols;
    return 0 ;
}

unsigned Matrix::rows() const
{
    if (d())
        return d()->rows;
    return 0 ;
}

bool Matrix::isEmpty() const
{
    if ( rows()==0 || cols()==0 )
        return true ;
    else
        return false ;
}

bool Matrix::isFixBinded() const
{
    if ( d() && ! d()->isIndependent() )
        return true ;
    else
        return false ;
}

bool Matrix::isParasitic() const
{
    if ( d() && d()->isParasitic() )
        return true ;
    else
        return false ;
}

bool Matrix::isReadOnly() const
{
    if ( d() && d()->isReadOnly() )
        return true ;
    else
        return false ;
}






const MatrixPrivate* Matrix::d() const
{
    return Priv ;
}
MatrixPrivate*& Matrix::d()
{
    return Priv ;
}

void Matrix::reset()
{
    if (d())
    {
        if (d()->isShared())
            d()->deShare() ;
        else
            delete d() ;
        d() = NULL ;
    }
}


void Matrix::fill(Float v)
{
    if ( ! (this->d()) )
        return ;

    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::fill():  ReadOnly Matrix should NOT call this method\n");
        return ;
    }

    this->d() = this->d()->checkWrite() ;
    MatrixPrivate & d = this->d()->getData() ;

    for ( unsigned i=0; i<d.rows; i++ )
        for ( unsigned j=0; j<d.cols; j++ )
            d(i,j) = v ;
}

Float Matrix::sum() const
{
    if ( rows()*cols() > 0 )
    {
        const MatrixPrivate & d = this->d()->getData() ;
        int i, j;
        Float sum = 0 ;
        for (i=0; i<d.rows; i++)
            for (j=0; j<d.cols; j++)
                sum += d(i,j) ;
        return sum ;
    }

    return 0 ;
}

Float Matrix::squareSum() const
{
    if ( rows()*cols() > 0 )
    {
        const MatrixPrivate & d = this->d()->getData() ;
        int i, j;
        Float sum = 0 ;
        for (i=0; i<d.rows; i++)
            for (j=0; j<d.cols; j++)
                sum += d(i,j) * d(i,j) ;
        return sum ;
    }

    return 0 ;
}

Matrix Matrix::clone() const
{
    if ( ! (this->d()) )
        return Matrix();

    MatrixPrivate* d = new MatrixPrivate(rows(), cols()) ;
    d->copyFrom(this->d());
    return Matrix(d) ;
}

Matrix Matrix::transpose(DataOption op) const
{
    if ( ! (this->d()) )
        return Matrix();

    MatrixPrivate*d  = new MatrixPrivate(cols(), rows(), (Matrix*)this, 0, 0, true) ;
    d->setReadOnly(true);
    switch(op)
    {
    case CloneData:
        return Matrix(d).clone() ;

    case ShareData:
    default:
        return Matrix(d);
    }
}

Matrix Matrix::transpose(DataOption op)
{
    if ( ! (this->d()) )
        return Matrix();

    MatrixPrivate*d  = new MatrixPrivate(cols(), rows(), this, 0, 0, true) ;
    switch(op)
    {
    case CloneData:
        return Matrix(d).clone() ;

    case ShareData:
    default:
        return Matrix(d);
    }
}


Matrix Matrix::subMat(unsigned i, unsigned j,
                      unsigned rows, unsigned cols,
                      DataOption op, bool transposeFirst) const
{
    if ( ! (this->d()) || rows*cols == 0 )
        return Matrix();

    MatrixPrivate*d  = new MatrixPrivate(rows, cols, (Matrix*)this, i, j, transposeFirst) ;
    d->setReadOnly(true);
    switch(op)
    {
    case CloneData:
        return Matrix(d).clone() ;

    case ShareData:
    default:
        return Matrix(d);
    }
}

Matrix Matrix::subMat(unsigned i, unsigned j,
                      unsigned rows, unsigned cols,
                      DataOption op, bool transposeFirst)
{
    if ( ! (this->d()) || rows*cols == 0 )
        return Matrix();

    MatrixPrivate*d  = new MatrixPrivate(rows, cols, this, i, j, transposeFirst) ;
    switch(op)
    {
    case CloneData:
        return Matrix(d).clone() ;

    case ShareData:
    default:
        return Matrix(d);
    }
}

Matrix Matrix::row(int i, DataOption op) const
{
    return subMat(i, 0, 1, cols(), op) ;
}

Matrix Matrix::row(int i, DataOption op)
{
    return subMat(i, 0, 1, cols(), op) ;
}

Matrix Matrix::col(int j, DataOption op) const
{
    return subMat(0, j, rows(), 1, op) ;
}

Matrix Matrix::col(int j, DataOption op)
{
    return subMat(0, j, rows(), 1, op) ;
}

Matrix::operator double() const
{
    return (*this)(0,0) ;
}

Matrix &Matrix::operator=(double s)
{
    fill(s);
    return *this;
}

Matrix &Matrix::operator=(const Matrix &A)
{
    if ( A.d() == d() )
        return *this ;

    if ( d() && ( d()->isParasitic() || ! d()->isIndependent() ) )
    {
        if ( isReadOnly() )
        {
            DEBUG_ERR("Matrix::operator=:  Unindepent ReadOnly Matrix should NOT be assigned\n");
            return *this;
        }

        this->d() = this->d()->checkWrite() ;
        d()->copyFrom (A.d()) ;
        return *this ;
    }

    reset() ;
    if ( A.d() )
        d() = A.d()->share() ;

    return *this ;
}



Matrix::Element Matrix::operator()(unsigned i, unsigned j)
{
    return Element(this,i,j);
}

Matrix::ElementRO Matrix::operator()(unsigned i, unsigned j) const
{
    return ElementRO (this,i,j);
}


Matrix& Matrix::Add ( const Matrix& _A, const Matrix& _B, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _B.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Add: the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Add(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }


    const MatrixPrivate &A = _A.d()->getData() ;
    const MatrixPrivate &B = _B.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( A.rows != B.rows || A.cols != B.cols ||
         A.rows != C.rows || A.cols != C.cols  )
    {
        DEBUG_ERR("Matrix::Add: unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) + B(i,j) ;

    return _C ;
}

Matrix& Matrix::Sub ( const Matrix& _A, const Matrix& _B, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _B.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Sub: the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Sub(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    const MatrixPrivate &B = _B.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( A.rows != B.rows || A.cols != B.cols ||
         A.rows != C.rows || A.cols != C.cols  )
    {
        DEBUG_ERR("Matrix::Sub: unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) - B(i,j) ;

    return _C ;
}

Float Matrix::Dot ( const Matrix& _A, const Matrix& _B)
{
    if ( NULL == _A.d() || NULL == _B.d() )
    {
        return 0 ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    const MatrixPrivate &B = _B.d()->getData() ;

    if ( A.cols != B.cols || A.rows != B.rows )
    {
        DEBUG_ERR("Matrix::Dot(A,B): unmatched dimensions \n");
        return 0 ;
    }


    Float dot = 0 ;
    for ( int i=0; i<A.rows; i++ )
        for ( int j=0; j<A.cols; j++ )
            dot += A(i,j) * B(i,j) ;

    return dot ;
}

Matrix& Matrix::Mult ( const Matrix& _A, const Matrix& _B, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _B.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Mult(A,B,C): the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Mult(): the Result Matrix should NOT be ReadOnly\n");
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    const MatrixPrivate &B = _B.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( C.cols != B.cols || C.rows != A.rows ||
         A.cols != B.rows  )
    {
        DEBUG_ERR("Matrix::Mult(A,B,C): unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
        {
            C(i,j) = 0 ;
            for ( int k=0; k<A.cols; k++ )
                C(i,j) += A(i,k) * B(k,j) ;
        }

    return _C ;
}

Matrix& Matrix::Add ( const Matrix& _A, Float b, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Add(A,b,C): the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Add(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( C.cols != A.cols || C.rows != A.rows )
    {
        DEBUG_ERR("Matrix::Add(A,b,C): unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) + b ;

    return _C ;
}

Matrix& Matrix::Sub ( const Matrix&_A, Float b, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Sub(A,b,C): the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Sub(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( C.cols != A.cols || C.rows != A.rows )
    {
        DEBUG_ERR("Matrix::Sub(A,b,C): unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) - b ;

    return _C ;
}


Matrix& Matrix::Mult ( const Matrix& _A, Float b, Matrix& _C)
{
    if ( NULL == _A.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Mult(A,b,C): the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Mult(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( C.cols != A.cols || C.rows != A.rows )
    {
        DEBUG_ERR("Matrix::Mult(A,b,C): unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) * b ;

    return _C ;
}


Matrix& Matrix::Div ( const Matrix& _A, Float b, Matrix& _C )
{
    if ( NULL == _A.d() || NULL == _C.d()  )
    {
        DEBUG_ERR("Matrix::Div(A,b,C): the input and output matrix should NOT be empty \n");
        return _C ;
    }
    if ( _C.isReadOnly() )
    {
        DEBUG_ERR("Matrix::Div(): the Result Matrix should NOT be ReadOnly\n");
        return _C ;
    }

    const MatrixPrivate &A = _A.d()->getData() ;
    _C.d() = _C.d()->checkWrite() ;
    MatrixPrivate &C = _C.d()->getData() ;

    if ( C.cols != A.cols || C.rows != A.rows )
    {
        DEBUG_ERR("Matrix::Div(A,b,C): unmatched dimensions \n");
        return _C ;
    }

    for ( int i=0; i<C.rows; i++ )
        for ( int j=0; j<C.cols; j++ )
            C(i,j) = A(i,j) / b ;

    return _C ;
}


Matrix operator+(const Matrix &A, const Matrix &B)
{
    if ( A.rows() * A.cols() == 0 || B.rows() * B.cols() == 0 )
    {
        DEBUG_ERR("operator+(const Matrix &A, const Matrix &B):  Empty input Matrix!\n");
        return Matrix();
    }

    if ( A.rows() != B.rows() || A.cols() != B.cols()  )
    {
        DEBUG_ERR("operator+(const Matrix &A, const Matrix &B):  unmatched dimensions!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Add(A,B,C) ;
}

Matrix operator+(const Matrix &A, Float b)
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator+(const Matrix &A, Float b):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Add(A,b,C) ;
}

Matrix operator-(const Matrix &A, const Matrix &B)
{
    if ( A.rows() * A.cols() == 0 || B.rows() * B.cols() == 0 )
    {
        DEBUG_ERR("operator-(const Matrix &A, const Matrix &B):  Empty input Matrix!\n");
        return Matrix();
    }

    if ( A.rows() != B.rows() || A.cols() != B.cols()  )
    {
        DEBUG_ERR("operator-(const Matrix &A, const Matrix &B):  unmatched dimensions!\n");
        return Matrix();
    }


    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Sub(A,B,C) ;
}

Matrix operator-(const Matrix &A, Float b)
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator-(const Matrix &A, Float b):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Sub(A,b,C) ;
}

Matrix operator-(const Matrix &A)
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator-(const Matrix &A):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Mult(A,-1,C) ;
}

Matrix operator*(const Matrix &A, const Matrix &B)
{
    if ( A.rows() * A.cols() == 0 || B.rows() * B.cols() == 0 )
    {
        DEBUG_ERR("operator*(const Matrix &A, const Matrix &B):  Empty input Matrix!\n");
        return Matrix();
    }

    if ( A.cols() != B.rows() )
    {
        DEBUG_ERR("operator*(const Matrix &A, const Matrix &B):  unmatched dimensions!\n");
        return Matrix();
    }

    Matrix C(A.rows(),B.cols()) ;
    return Matrix::Mult(A,B,C) ;
}

Matrix operator*(const Matrix &A, const Vector3D &B)
{
    return A*Matrix(B);
}

Matrix operator*(const Vector3D &A, const Matrix &B)
{
    return Matrix(A)*B;
}

Matrix operator* (const Matrix &A, Float b)
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator*(const Matrix &A, Float b):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Mult(A,b,C) ;
}

Matrix operator* (Float b, const Matrix &A )
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator*(Float b, const Matrix &A ):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Mult(A,b,C) ;
}

Matrix operator/ (const Matrix &A, Float b)
{
    if ( A.rows() * A.cols() == 0 )
    {
        DEBUG_ERR("operator/(const Matrix &A, Float b):  Empty input Matrix!\n");
        return Matrix();
    }

    Matrix C(A.rows(),A.cols()) ;
    return Matrix::Div(A,b,C) ;
}


Matrix& Matrix::operator+=(const Matrix &B)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator+= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }

    Matrix &A = *this ;
    return Matrix::Add(A,B,A) ;
}

Matrix& Matrix::operator+=(Float b)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator+= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }
    Matrix &A = *this ;
    return Matrix::Add(A,b,A) ;
}

Matrix& Matrix::operator-=(const Matrix &B)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator-= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }
    Matrix &A = *this ;
    return Matrix::Sub(A,B,A) ;
}

Matrix& Matrix::operator-=(Float b)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator-= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }
    Matrix &A = *this ;
    return Matrix::Sub(A,b,A) ;
}

Matrix& Matrix::operator*=(Float b)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator*= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }
    Matrix &A = *this ;
    return Matrix::Mult(A,b,A) ;
}

Matrix& Matrix::operator/=(Float b)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::operator/= : ReadOnly Matrix should NOT use this operator\n");
        return *this ;
    }
    Matrix &A = *this ;
    return Matrix::Div(A,b,A) ;
}


///////////////////////////////////////////////////////
////
///    MATRIX_IO Part

#ifdef ENABLE_MATRIX_IO

int Matrix::printMat() const
{
    int maxSize = rows()*cols()*20 + 1 ;
    char* str = new char[maxSize] ;
    int ret = save2Str(str, maxSize) ;
    printf(str) ;
    delete[] str;
    return ret ;
}

int Matrix::save2Str (char* str, int bufLen) const
{
    if ( rows()*cols()==0 )
    {
        if (bufLen>0 || str)
            *str = '\0' ;
        return 0 ;
    }

    int totalBufLen = bufLen ;
    const MatrixPrivate & d = this->d()->getData() ;
    int strLen = 0 ;
    for (int i=0; i < d.rows; i++)
    {
        for (int j = 0; j < d.cols; j++)
        {
            strLen = snprintf(str, bufLen, "%5.6f\t", d(i,j)) ;
            if ( strLen < bufLen )
            {
                str += strLen ;
                bufLen -= strLen ;
            }
            else
            {
                return -1 ;
            }
        }

        strLen = snprintf(str, bufLen, "\n") ;
        if ( strLen < bufLen )
        {
            str += strLen ;
            bufLen -= strLen ;
        }
        else
        {
            return -1 ;
        }
    }

    return totalBufLen - bufLen ;
}


static bool isNumber(char c)
{
    return  c>='0' && c<='9' ;
}

static bool str2float(const char* str, int strLen, Float* f)
{
    bool numCount = 0 ;
    int dotCount = 0 ;
    int eCount = 0 ;

    if (strLen == 0)
        return false ;

    if (str[0]!='+' && str[0]!='-' && !isNumber(str[0]))
        return false ;
    if (isNumber(str[0]))
        numCount ++ ;

    for (int i=1; i<strLen; i++)
    {
        if ( isNumber(str[i]) )
        {
            numCount++ ;
            continue ;
        }

        // 'e' only allowed to appear once (after at least one number)
        if (str[i] == 'e')
        {
            if (numCount>0 && eCount<1)
            {
                eCount ++ ;
                continue ;
            }
            else
            {
                return false ;
            }
        }

        // inner the number, the '+' and '-' only allowed to appear next to (after) 'e' .
        if (str[i]=='+' || str[i]=='-')
        {
            if (str[i-1]=='e' && i+1<strLen && isNumber(str[i+1]))
                continue ;
            else
                return false ;
        }

        // dot only allowed to appear once (after at least one number and before 'e')
        if (str[i]=='.')
        {
            if (numCount>0 && eCount<1 && dotCount<1)
            {
                dotCount ++ ;
                continue ;
            }
            else
                return false ;
        }
    }

    *f = atof(str) ;
    return true ;
}

bool Matrix::loadMat (const char* _str, int strLen)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::LoadMatFile() : ReadOnly Matrix should NOT use this method\n");
        return false ;
    }

    if (strLen==0)
        strLen = strlen(_str) ;

    if (strLen<=0)
        return false ;

    int bufLen = strLen + 1 ;
    char* sub = new char[bufLen] ;
    Float* num = new Float[bufLen] ;
    int numCount = 0 ;
    bool newRowStarted = false ;
    int subLen = 0 ;
    int curRowSize = 0 ;
    int preRowSize = -1 ;
    int curRows = 0 ;
    int i,j ;
    for (i=0,j=0; i<strLen; i++)
    {
        if ( newRowStarted )
        {
            bool separator = _str[i]=='\r' || _str[i]=='\t' || _str[i]==','
                    || _str[i]=='\n' || _str[i]==' ' || _str[i]<0  ;

            if ( !separator )
            {
                sub[subLen] = _str[i] ;
                subLen ++ ;
            }

            if (separator || i==strLen-1)
            {
                if (subLen > 0)
                {
                    sub[subLen] = '\0' ;
                    bool fValid = str2float(sub, subLen, &(num[numCount])) ;
                    if (fValid)
                    {
                        numCount ++ ;
                        curRowSize ++ ;
                    }
                    else
                    {
                        delete[] sub ;
                        delete[] num ;
                        return false ;
                    }
                }

                subLen = 0 ;
                if ( _str[i]=='\n' )
                {
                    if (curRows > 1)
                    {
                        if (curRowSize != preRowSize)
                        {
                            delete[] sub ;
                            delete[] num ;
                            return false ;
                        }
                    }
                    preRowSize = curRowSize ;
                    newRowStarted = false;
                }
            }

        }
        else
        {
            bool separator = _str[i]=='\r' || _str[i]=='\t' || _str[i]==','
                    || _str[i]=='\n' || _str[i]==' ' || _str[i]<0  ;

            if ( !separator )
            {
                newRowStarted = true ;
                curRows ++ ;
                curRowSize = 0 ;
                subLen = 0 ;
                sub[subLen] = _str[i] ;
                subLen ++ ;
            }
        }
    }


    *this = Matrix(curRows, curRowSize, num).clone() ;
    delete[] sub ;
    delete[] num ;
    return true ;
}


int Matrix::save2file(const char *filename) const
{
    FILE *fp;

    if (NULL==(fp=fopen(filename, "wb"))) {
        puts("Open data file failed.  Exit.");
        return false ;
    }

    int maxSize = rows()*cols()*20 + 1 ;
    char* str = new char[maxSize] ;
    int ret = save2Str(str, maxSize) ;
    fwrite(str, 1, strlen(str), fp) ;
    delete[] str;

    fclose(fp);
    return ret ;
}

bool Matrix::LoadMatFile(const char *filename)
{
    if ( isReadOnly() )
    {
        DEBUG_ERR("Matrix::LoadMatFile() : ReadOnly Matrix should NOT use this method\n");
        return false ;
    }

    FILE *fp;
    int fileLen ;
    char *buf ;
    if (NULL==(fp=fopen(filename, "r")))
    {
        puts("Open data file failed.  Exit.");
        return  false;
    }

    fseek(fp,0,SEEK_END);
    fileLen = ftell(fp) ;
    fseek(fp,0,SEEK_SET);
    buf  =  (char*)(malloc(fileLen+1)) ;
    fread(buf,1,fileLen,fp) ;
    fclose(fp) ;
    buf[fileLen] = '\0' ;
    return loadMat(buf) ;
}


Matrix Matrix::I(int n)
{
    Matrix I(n,n) ;
    I.fill(0);
    MatrixPrivate &d = I.d()->getData() ;
    for (int i=0; i<n; i++)
        d(i,i) = 1 ;
    return I ;
}

Float Matrix::trace() const
{
    if ( isEmpty() || rows() != cols() )
        return 0 ;

    const MatrixPrivate & d = this->d()->getData() ;
    Float trace = 0 ;
    for (int i=0; i<d.rows && i<d.cols; i++)
        trace += d(i,i) ;

    return trace ;
}

Float Matrix::det() const
{
    if ( isEmpty() || rows() != cols() )
        return 0 ;

    Matrix tmp = this->clone() ;
    int sign = 1 ;
    if (!tmp.GaussElimination(true, NULL, &sign))
    {
        return 0 ;
    }

    Float det = sign ;
    const MatrixPrivate & tmpd = tmp.d()->getData() ;
    for (int i=0; i<rows(); i++)
    {
        det *= tmpd(i,i) ;
    }

    return det ;
}

Matrix Matrix::getInverse() const
{
    if ( isEmpty() || rows() != cols() )
        return Matrix() ;

    Matrix tmp = this->clone() ;
    Matrix inv = Matrix::I(rows()) ;

    if (!tmp.GaussElimination(true, &inv))
    {
        return Matrix() ;
    }

    for (int i=0; i<rows(); i++)
    {
        inv.row(i) /= tmp(i,i) ;
    }

    return inv ;
}

bool Matrix::GaussElimination ( bool backEnable, Matrix* follow, int* sign)
{
    this->d() = this->d()->checkWrite() ;
    MatrixPrivate &d = this->d()->getData() ;
    if ( follow )
        follow->d() = follow->d()->checkWrite() ;


    bool ret = true ;
    int i, j ;
    if ( sign )
        *sign = 1 ;

    for ( i=0; i<rows(); i++ )
    {
        if ( i >= cols() )     // 多余的行不再处理。比如共有8行5列，则第5行以后不再处理，因为消元法进行到第5行后，下面的3行已经全变成0了
        {
            break ;
        }

        // 寻找第 i 列的第一个不为 0 值, 放到第 i 行
        for ( j=i; j<rows(); j++ )
        {
            if ( d(j,i) >= FLT_MIN || d(j,i) <= - FLT_MIN  )
                break ;
            else
                d(j,i) = 0.0f ;
        }

        if ( j == rows() )
        {
            // 第 i 列全为0，方程组无解哦，或者行列式为 0，或者不存在逆矩阵 。
            // 跳过第 i 行
            ret = false ;
            continue ;
        }
        else if ( j > i )
        {
            // swap row i and row j

            Matrix temp = row(i).clone() ;
            row(i) = row(j) ;
            row(j) = temp ;

            if (follow)
            {
                temp = follow->row(i).clone() ;
                follow->row(i) = follow->row(j) ;
                follow->row(j) = temp ;
            }

            if (sign)
            {
                *sign  *=  -1 ;   // sign 取反
            }
        }

        // 开始消元，消掉第 i 行以下各行的第 i 列。如果 backEnable 使能，则第 i 行以上各行的第 i 列也会被消掉
        int start ;
        if ( backEnable )
            start = 0 ;
        else
            start = i+1 ;
        for ( j=start; j<rows(); j++ )
        {
            if (j==i)
                continue ;

            //  x[j] = x[i]*ratio + x[j] ;
            Float ratio = -d(j,i) / d(i,i) ;
            row(j) = row(i)*ratio + row(j) ;
            if ( follow )
                follow->row(j) = follow->row(i)*ratio + follow->row(j) ;
            d(j,i) = 0 ;
        }
    }

    return ret ;
}


Matrix Matrix::solve(const Matrix& _A, const Matrix& _Y)
{
    int n = _A.rows() ;
    if ( _A.rows() != _A.cols() || _A.rows() == 0 || _A.rows() != _Y.rows() )
    {
        return Matrix() ;
    }

    Matrix A = _A.clone() ;
    Matrix Y = _Y.clone() ;
    Matrix X (_Y.rows(), _Y.cols()) ;

    bool resultOK = A.GaussElimination ( true, &Y ) ;
    if ( resultOK )
    {
        for ( int i=0; i<n; i++ )
            X.row(i) = Y.row(i) / A(i,i) ;
        return X ;
    }
    else
    {
        return Matrix() ;
    }
}



#endif
