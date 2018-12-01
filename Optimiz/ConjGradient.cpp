#include "Optimiz/ConjGradient.h"
#include "Optimiz/LineSearch.h"
#include <cmath>

IMPORT_JK_Math
IMPORT_JK_OPTIMIZ

ConjGradient::ConjGradient(SingleFunction* cost,
       Float _gradEps, Float _lineSearchEps, int _maxLoops) :
    Method(cost), gradEps(_gradEps),
    lineSearchEps(_lineSearchEps), maxLoops(_maxLoops)
{

}

ConjGradient::~ConjGradient()
{

}

Result ConjGradient::solve(const Matrix &_startPoint)
{
    int n = cost->inputDim() ;
    Result res ;
    Matrix startPoint = _startPoint.clone() ;
    LineSearch& lineSearch = *searcher ;
    Float probeStep = 0.001 ;
    Float gradEps2 = gradEps * gradEps ;
    res.w = startPoint ;
    res.cost = cost->Func(res.w) ;

    int loops = 0 ;
    Matrix d(n,1) ;
    Float preGradNorm = 1 ;

    printf("ConjGradient::solve(): %d\t%7.8f\t%11.12f\n", loops, 0.0, res.cost);

    while (loops++ < maxLoops)
    {
        d.fill(0);

        Matrix oldW = res.w ;
        for (int i=0; i<n; i++)
        {
            Matrix grad = cost->Grad(res.w) ;
            Float gradNorm = grad.squareSum() ;
            if (gradNorm < gradEps2)
            {
                res.grad = grad ;
                res.valid = true ;
                return res ;
            }

            if (i==0)
            {
                d = -grad.transpose() ;
            }
            else
            {
                Float numerator = grad*d ;
                Float beta = gradNorm / preGradNorm ;
                Float theta = numerator / preGradNorm ;
                d = - (1+theta)*grad.transpose() + beta*d ;
            }
            preGradNorm = gradNorm ;

            Matrix tmpw = lineSearch(cost, d, res.w, probeStep, lineSearchEps) ;
            Float tmpcost = cost->Func(tmpw) ;
            if (tmpcost <= res.cost)
            {
                Float nextProbeStep = 0.5 * sqrt((tmpw - res.w).squareSum()) ;
                if (nextProbeStep != 0)
                    probeStep = nextProbeStep ;
                res.w = tmpw ;
                res.cost = tmpcost ;
            }

            else if (tmpcost > res.cost)
            {
                //printf("ConjGradient::solve() : ascent back \n") ;
                //res.grad = cost->Grad(res.w) ;
                //res.print();
                res.valid = true ;
                return res ;
            }
        }

        //printf("ConjGradient::solve(): %d\t%7.8f\t%11.12f\n", loops, sqrt((res.w-oldW).squareSum()), res.cost);

        if ((res.w-oldW).squareSum()==0)
        {
            //printf("ConjGradient::solve() Stuck !!! \n") ;
            break ;
        }
    }

    //printf("ConjGradient::solve() not converged exactly\n") ;
    //res.grad = cost->Grad(res.w) ;
    //res.print();
    res.valid = true ;
    return res ;
}
