#include "Optimiz/Optimiz.h"
#include "cmath"

IMPORT_JK_OPTIMIZ
IMPORT_JK_Math

static const Float DEFAULT_EPS = sizeof(Float) == sizeof(double) ? 0.000001 : 0.001 ;

///////////////////////////////////////////////////////////////////////////
///
///     SingleFunction Area
///

Matrix SingleFunction::Grad(const Matrix &w)
{
    Matrix eps (w.rows(), 1) ;
    eps.fill(DEFAULT_EPS);
    return LazyGrad(w, eps) ;
}

Matrix SingleFunction::Hess(const Matrix &w)
{
    Matrix eps (w.rows(), 1) ;
    eps.fill(DEFAULT_EPS);
    return LazyHess (w, eps) ;
}

Matrix SingleFunction::LazyGrad(const Matrix &w, const Matrix &eps)
{
    int i, n=w.rows() ;
    Matrix grad (1,n) ;
    Float f = Func(w) ;
    bool minTouched = true ;
    bool maxTouched = true ;

    for (i=0; i<n; i++)
    {
        Matrix w1 = w ;
        Matrix w2 = w ;
        w1(i,0) = w(i,0) + eps(i,0) ;
        w2(i,0) = w(i,0) - eps(i,0) ;
        Float f1 = Func(w1) ;
        Float f2 = Func(w2) ;

        if ( minTouched &&
             (f1 < f || f2 < f) )
            minTouched = false ;
        if ( maxTouched &&
             (f1 > f || f2 > f) )
            maxTouched = false ;

        grad(0,i) = 0.5 * (f1 - f2) / eps(i,0) ;
    }

    if ( minTouched || maxTouched )
        grad.fill(0);

    return grad ;
}

Matrix SingleFunction::LazyHess(const Matrix &w, const Matrix &eps)
{
    int i, n=w.rows() ;
    Matrix hess (n,n) ;
    for (i=0; i<n; i++)
    {
        Matrix w1 = w ;
        Matrix w2 = w ;
        w1(i,0) = w(i,0) + eps(i,0) ;
        w2(i,0) = w(i,0) - eps(i,0) ;
        Matrix g1 = Grad(w1) ;
        Matrix g2 = Grad(w2) ;
        hess.col(i) = 0.5 * (g1 - g2).transpose() / eps(i,0) ;
    }
    return hess ;
}

bool SingleFunction::checkInputDim (const Matrix& w) const
{
    if (inputDim() == 0 || w.rows() != inputDim() || w.cols() != 1)
        return false ;
    else
        return true ;
}


///////////////////////////////////////////////////////////////////////////
///
///     Function Area
///

Matrix Function::Jacobian(const Matrix& w)
{
    Matrix eps (w.rows(), 1) ;
    eps.fill(DEFAULT_EPS);
    return LazyJacobian(w, eps) ;
}

Matrix Function::Hessian (const Matrix& w)
{
    Matrix eps (w.rows(), 1) ;
    eps.fill(DEFAULT_EPS);
    return LazyHessian(w, eps) ;
}

Matrix Function::LazyJacobian(const Matrix& w, const Matrix& eps)
{
    Matrix f = Func(w) ;
    int n=w.rows() ;
    int m=f.rows() ;
    Matrix jaco (m,n) ; // the Jacobian Matrix is in m rows and n cols .

    int i, k;
    bool* minTouched = new bool[m] ;
    bool* maxTouched = new bool[m] ;
    for (k=0; k<m; k++)
    {
        minTouched[k] = maxTouched[k] = true ;
    }

    for (i=0; i<n; i++)
    {
        Matrix w1 = w ;
        Matrix w2 = w ;
        w1(i,0) = w(i,0) + eps(i,0) ;
        w2(i,0) = w(i,0) - eps(i,0) ;
        Matrix f1 = Func(w1) ;
        Matrix f2 = Func(w2) ;

        for (k=0; k<m; k++)
        {
            if ( minTouched[k] &&
                 (f1(k,0) < f(k,0) || f2(k,0) < f(k,0)) )
                minTouched[k] = false ;
            if ( maxTouched[k] &&
                 (f1(k,0) > f(k,0) || f2(k,0) > f(k,0)) )
                maxTouched[k] = false ;
        }

        jaco.col(i) = (0.5 / eps(i,0)) * (f1 - f2) ;
    }

    for (k=0; k<m; k++)
    {
        if ( minTouched[k] || maxTouched[k] )
            jaco.row(k).fill(0);
    }

    return jaco ;
}

Matrix Function::LazyHessian(const Matrix& w, const Matrix& eps)
{
    int m = outputDim() ;
    int n = inputDim() ;
    Matrix hess (n*m, n) ;

    int i, k;
    for (i=0; i<n; i++)
    {
        Matrix w1 = w ;
        Matrix w2 = w ;
        w1(i,0) = w(i,0) + eps(i,0) ;
        w2(i,0) = w(i,0) - eps(i,0) ;
        Matrix J1 = Jacobian(w1) ;
        Matrix J2 = Jacobian(w2) ;
        Matrix dJ = (0.5 / eps(i,0)) * (J1 - J2).transpose() ;
        Matrix col (n*m,1) ;
        for (k=0; k<m; k++)
        {
            col.subMat(k*n, 0, n, 1) = dJ.col(k) ;
        }
        hess.col(i) = col ;
    }
    return hess ;
}

bool Function::checkInputDim (const Matrix& w) const
{
    if (inputDim() == 0 || w.rows() != inputDim() || w.cols() != 1)
        return false ;
    else
        return true ;
}

///////////////////////////////////////////////////////////////////////////
///
///     LagrangePenalty Area
///

int LagrangePenalty::inputDim () const
{
    if (cost)
        return cost->inputDim() ;
    else
        return 0 ;
}

Float LagrangePenalty::Func(const Matrix& w)
{
    if (!checkInputDim(w))
        return 0 ;

    Float ret = 0 ;
    err = 0 ;
    if (cost)
        ret += cost->Func(w) ;

    if (EQ && EQ->outputDim() == beta.rows() && beta.cols() == 1)
    {
        Matrix H = EQ->Func(w) ;
        ret += (beta + 0.5*sigma*H).transpose() * H ;
        err += H.transpose()*H ;
        nextBeta = beta + sigma*H ;
    }

    if (NEQ && NEQ->outputDim() == alpha.rows() && alpha.cols() == 1)
    {
        Matrix G = NEQ->Func(w) ;
        int n = G.rows() ;
        Float r = 0.5 / sigma ;
        for (int i=0; i<n; i++)
        {
            Float u = alpha(i,0) ;
            Float g = G(i,0) ;
            Float c = sigma ;
            Float t = c*g+u ;
            if ( t > 0 )
            {
                ret += r*(t*t - u*u) ;
                err += g*g ;
                nextAlpha(i,0) = alpha(i,0) + sigma*g ;
            }
            else
            {
                ret += r*(-u*u) ;
                err += (-u/c)*(-u/c) ;
                nextAlpha(i,0) = 0 ;
            }
        }
    }

    return ret ;
}

Matrix LagrangePenalty::Grad(const Matrix& w)
{
    if (!checkInputDim(w))
        return Matrix() ;

    Matrix grad (1,w.rows()) ;
    grad.fill(0);

    if (cost)
        grad += cost->Grad(w) ;

    if (EQ && EQ->outputDim() == beta.rows() && beta.cols() == 1)
    {
        Matrix jH = EQ->Jacobian(w) ;
        Matrix H = EQ->Func(w) ;
        grad += (beta + sigma*H).transpose() * jH ;
    }

    if (NEQ && NEQ->outputDim() == alpha.rows() && alpha.cols() == 1)
    {
        Matrix G = NEQ->Func(w) ;
        Matrix jG = NEQ->Jacobian(w) ;

        int n = G.rows() ;
        Float r = 0.5 / sigma ;
        for (int i=0; i<n; i++)
        {
            Float u = alpha(i,0) ;
            Float g = G(i,0) ;
            Matrix jg = jG.row(i) ;
            Float c = sigma ;
            Float t = c*g+u ;
            if ( t > 0 )
            {
                grad += (u+c*g)*jg ;
            }
        }
    }

    return grad ;
}



///////////////////////////////////////////////////////////////////////////
///
///     Result Area
///

Result::Result()
{
    valid = false ;
}

Result::Result(const Result& other)
{
    this->valid = other.valid ;
    this->cost = other.cost ;
    this->grad = other.grad ;
    this->hess = other.hess ;
    this->w = other.w ;
}

Result& Result::operator= (const Result& other)
{
    this->valid = other.valid ;
    this->cost = other.cost ;
    this->grad = other.grad ;
    this->hess = other.hess ;
    this->w = other.w ;
    return *this;
}

void Result::print() const
{
    printf("Optimiz::Result::print():\n");
    printf("w: \n");
    w.transpose().printMat() ;
    printf("grad: \n");
    grad.printMat() ;    
    printf("gradLen: %11.12lf\n", sqrt(grad.squareSum())) ;
    printf("cost: %11.12lf\n", cost) ;
}

///////////////////////////////////////////////////////////////////////////
///
///     ConstrainedMethod Area
///

ConstrainedMethod::ConstrainedMethod (
        SingleFunction* _cost, Function* _EQ,
        Function* _NEQ, Method* _method ) :
    iterMethod(_method), cost(_cost), EQ(_EQ), NEQ(_NEQ)
{
    sigmaRate = 5 ;
    constraint_eps = DEFAULT_EPS ;
    sigmaRaiseThr = 0.5 ;
    MaxIter = 1000 ;
}

Result ConstrainedMethod::solve(const Matrix& startPoint,
          const Matrix& startBeta, const Matrix& startAlpha,
               const Float startSigma )
{
    if (!iterMethod || !cost)
    {
        return Result() ;
    }

    Matrix startP = startPoint.clone() ;
    if ( startP.isEmpty() )
    {
        startP = Matrix(cost->inputDim(), 1);
        startP.fill(0);
    }

    if (!EQ && !NEQ)
    {
        iterMethod->cost = cost ;
        return iterMethod->solve(startP) ;
    }

    LagrangePenalty LP(cost, EQ, NEQ) ;
    LP.alpha = startAlpha.clone() ;
    LP.beta = startBeta.clone() ;
    LP.sigma = startSigma ;
    if ( EQ && LP.beta.isEmpty() )
    {
        LP.beta = Matrix(EQ->outputDim(),1) ;
        LP.beta.fill(1);
    }
    if ( NEQ && LP.alpha.isEmpty() )
    {
        LP.alpha = Matrix(NEQ->outputDim(),1);
        LP.alpha.fill(1);
    }
    if ( LP.sigma == 0 )
    {
        LP.sigma = 0.1 ;
    }

    iterMethod->cost = &LP ;
    Result curResult ;
    curResult.w = startP ;
    curResult.cost = LP.Func(curResult.w) ;
    Float err = -1 , preErr = -1 ;

    int iter = 0 ;
    Float constraint_eps2 = constraint_eps*constraint_eps ;
    Float oldSigma = LP.sigma ;
    bool convergingStarted = false ;
    while (iter++ < MaxIter)
    {
        curResult = iterMethod->solve(curResult.w) ;
        preErr = err ;
        err = LP.err ;

        //printf("***** ConstrainedMethod::solve() innerLoop *****\n");
        //printf("iter %d, err %11.12lf,  sigma %f\n", iter, LP.err, LP.sigma);
        //printf("beta : ");
        //LP.beta.transpose().printMat() ;
        //printf("alpha : ");
        //LP.alpha.transpose().printMat() ;

        if (err < constraint_eps2)
        {
            printf("ConstrainedMethod::solve() Converged Exactly.\n");
            //curResult.grad = LP.Grad(curResult.w) ;
            //curResult.hess = LP.Hess(curResult.w) ;
            //curResult.print();
            return curResult ;
        }

        LP.alpha = LP.nextAlpha ;
        LP.beta = LP.nextBeta ;

        if ( oldSigma != LP.sigma )
        {
            oldSigma = LP.sigma ;
            if ( convergingStarted && err / preErr > 0.8*0.8 )
            {
                printf("ConstrainedMethod::solve() Converged Too Slow.\n");
                //curResult.grad = LP.Grad(curResult.w) ;
                //curResult.hess = LP.Hess(curResult.w) ;
                //curResult.print();
                return curResult ;  // converging too slow .
            }
        }

        if ( preErr < 0 )
            continue ;

        if ( err / preErr > sigmaRaiseThr*sigmaRaiseThr )
            LP.sigma *= sigmaRate ;
        else
            convergingStarted = true ;
    }

    //curResult.grad = LP.Grad(curResult.w) ;
    //curResult.hess = LP.Hess(curResult.w) ;
    //curResult.print();
    return curResult ;
}
