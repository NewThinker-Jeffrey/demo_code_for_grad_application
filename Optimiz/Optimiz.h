#ifndef OPTIMIZ_H
#define OPTIMIZ_H

#include "MathToolsCpp/Matrix.h"

#define START_JK_OPTIMIZ       namespace JK { namespace Math { namespace Optimiz {
#define END_JK_OPTIMIZ                          }}}
#define IMPORT_JK_OPTIMIZ      using namespace JK::Math::Optimiz ;
#define JKOptimiz   JK::Math::Optimiz


START_JK_OPTIMIZ

struct Function ;
struct SingleFunction ;

struct Result
{
public :
    Result() ;
    Result(const Result& other) ;
    Result& operator= (const Result& other) ;
    virtual ~Result() {}
    virtual void print() const ;

public :
    Matrix w ;
    Matrix grad ;
    Matrix hess ;
    Float cost ;
    bool valid ;
};


struct Method
{
public :
    Method(SingleFunction* c=NULL) : cost(c) {}
    virtual ~Method() {}
    virtual Result solve(const Matrix& startPoint=Matrix()) = 0 ;
public :
    SingleFunction* cost ;
} ;

struct ConstrainedMethod
{
public :
    ConstrainedMethod (SingleFunction* _cost=NULL,
                       Function* _EQ = NULL,
                       Function* _NEQ = NULL,
                       Method* _method = NULL ) ;
    virtual ~ConstrainedMethod() {}
    Result solve(const Matrix& startPoint=Matrix(),
                              const Matrix& startBeta=Matrix(),
                              const Matrix& startAlpha=Matrix(),
                              const Float startSigma = 1) ;

public :
    Method* iterMethod ;
    SingleFunction* cost ;
    Function* EQ ;
    Function* NEQ ;

public :
    Float sigmaRate ;
    Float constraint_eps ;
    Float sigmaRaiseThr ;
    int MaxIter ;
    Result result ;
} ;


struct Function
{
public :
    virtual ~Function() {}
    virtual int inputDim () const = 0 ;
    virtual int outputDim () const = 0 ;
    virtual Matrix Func(const Matrix& w) = 0 ;// input a Matrix in [nx1], output a Matrix in [mx1].
    virtual Matrix Jacobian(const Matrix& w) ;// input a Matrix in [nx1], output a Matrix in [mxn].
    virtual Matrix Hessian (const Matrix& w) ;// input a Matrix in [nx1], output a combined-Matrix in mx[nxn] .
protected :
    Matrix LazyJacobian(const Matrix& w, const Matrix& eps) ;
    Matrix LazyHessian(const Matrix& w, const Matrix& eps) ;
    bool checkInputDim (const Matrix& w) const ;
} ;

struct SingleFunction
{
public :
    virtual ~SingleFunction() {}
    virtual int inputDim () const = 0 ;
    virtual Float Func(const Matrix& w) = 0 ;// input a Matrix in [nx1], output a scalar .
    virtual Matrix Grad(const Matrix& w)    ;// input a Matrix in [nx1], output a Matrix in [1xn].
    virtual Matrix Hess (const Matrix& w)   ;// input a Matrix in [nx1], output a Matrix in [nxn].

protected :
    Matrix LazyGrad(const Matrix& w, const Matrix& eps) ;
    Matrix LazyHess(const Matrix& w, const Matrix& eps) ;
    bool checkInputDim (const Matrix& w) const ;
} ;

struct LagrangePenalty : public SingleFunction
{
public :
    LagrangePenalty(SingleFunction* _cost=NULL,
           Function* _EQ=NULL, Function* _NEQ=NULL ) :
        cost(_cost), EQ(_EQ), NEQ(_NEQ) {}

    ~LagrangePenalty() {}
    Float Func(const Matrix& w) ;
    Matrix Grad(const Matrix& w);
    int inputDim () const ;
public :
    SingleFunction* cost ;
    Function* EQ ;
    Function* NEQ ;
    Matrix alpha ;
    Matrix beta ;
    Float sigma ;
    Matrix nextAlpha ;
    Matrix nextBeta ;
    Float err ;
};





END_JK_OPTIMIZ

#endif // OPTIMIZ_H
