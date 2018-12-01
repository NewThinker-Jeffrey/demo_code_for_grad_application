

#include "Optimiz/LineSearch.h"
#include <cmath>

IMPORT_JK_Math
IMPORT_JK_OPTIMIZ

Matrix LineSearch::move (const Matrix& start, const Matrix& dir, Float stepLen)
{
    return start + dir*stepLen ;
}

Float LineSearch::distance2 (const Matrix& w1, const Matrix& w2)
{
    return (w1-w2).squareSum() ;
}

Matrix LineSearch::interpolate (const Matrix& start, const Matrix& end,
                   Float endWeight)
{
    return start*(1-endWeight) + end*endWeight ;
}

Matrix LineSearch::middleOf (const Matrix& start, const Matrix& end)
{
    return 0.5*start + 0.5*end ;
}

bool LineSearch::isStepAcceptable(const Matrix& start, const Matrix& end)
{
    return true ;
}

Matrix LineSearch::operator() (SingleFunction* func,
           const Matrix& dir, const Matrix& start,
           Float startStepLen, Float eps )
{
    Float eps2 = eps*eps ;
    Float Val[3] ;
    Matrix interval = findInterval(func, dir, start, startStepLen, Val) ;
    if (interval.isEmpty())
    {
        printf("LineSearch(): BadCase(*0*): Stuck !!! \n") ;
        return start.clone() ;
    }

    Matrix l, m, r, pl, pr ;
    Float Vl, Vm, Vr, Vpl, Vpr ;
    l = interval.col(0).clone() ;
    m = interval.col(1).clone() ;
    r = interval.col(2).clone() ;
    Vl = Val[0] ;
    Vm = Val[1] ;
    Vr = Val[2] ;

    Float localEps2 ;
    localEps2 = 0.01 * distance2(l,m) ;
    eps2 = eps2 < localEps2 ? eps2 : localEps2 ;
    localEps2 = 0.01 * distance2(r,m) ;
    eps2 = eps2 < localEps2 ? eps2 : localEps2 ;

    Float Distance2 = distance2(l,r) ;
    Float oldDistance2 = 0 ;
    while (Distance2 > eps2 && Distance2 != oldDistance2)
    {
        pl = middleOf(l,m) ;
        Vpl = func->Func(pl) ;
        if (Vpl < Vm)
        {
            r = m ; Vr = Vm ;
            m = pl; Vm = Vpl ;
            oldDistance2 = Distance2 ;
            Distance2 = distance2(l,r) ;
            continue ;
        }

        pr = middleOf(r,m) ;
        Vpr = func->Func(pr) ;
        if (Vpr < Vm)
        {
            l = m ; Vl = Vm ;
            m = pr; Vm = Vpr ;
            oldDistance2 = Distance2 ;
            Distance2 = distance2(l,r) ;
            continue ;
        }

        l = pl ; Vl = Vpl ;
        r = pr ; Vr = Vpr ;
        oldDistance2 = Distance2 ;
        Distance2 = distance2(l,r) ;
    }

    return m ;
}

Matrix LineSearch::findInterval (SingleFunction* func,
           const Matrix& dir, const Matrix& start,
           Float startStepLen, Float* pVal )
{
    int n = start.rows() ;
    Matrix interval (n, 3) ;
    Matrix l(n,1), m(n,1), r(n,1) ;
    Float vl, vm, vr ;

    Float stepLen = startStepLen ;
    Float startVal = func->Func(start) ;
    Matrix dirUnit = dir / sqrt(dir.squareSum()) ;

    l = move(start, dirUnit, stepLen) ;
    r = move(start, dirUnit, -stepLen) ;
    while ( ! isStepAcceptable(start, l) ||
            ! isStepAcceptable(start, r) )
    {
        stepLen *= 0.5 ;
        l = move(start, dirUnit, stepLen) ;
        r = move(start, dirUnit, -stepLen) ;
    }

    vl = func->Func(l) ;
    vr = func->Func(r) ;
    while ( vl >= startVal && vr >= startVal )
    {
        if (0.5*stepLen == stepLen)
        {
            vl = vm = vr = startVal ;
            if (pVal)
            {
                pVal[0] = vl; pVal[1] = vm; pVal[2] = vr;
            }
            return Matrix() ;
        }

        stepLen *= 0.5 ;
        l = move(start, dirUnit, stepLen) ;
        r = move(start, dirUnit, -stepLen) ;
        vl = func->Func(l) ;
        vr = func->Func(r) ;
    }

    if ( vl > vr )
    {
        vl = vr ;
        l = r ;
        dirUnit = -dirUnit ;
    }

    m = l; vm = vl ;
    l = start.clone() ; vl = startVal ;
    while (1)
    {
        stepLen *= 2 ;

        r = move(m, dirUnit, stepLen) ;
        if ( ! isStepAcceptable(start, r) )
        {
            printf("LineSearch::findInterval(): test point is beyond the boundary\n\t return unideal interval \n") ;
            while (! isStepAcceptable(start, r))
            {
                stepLen *= 0.5 ;
                r = move(m, dirUnit, stepLen) ;
            }
            vr = func->Func(r) ;
            break ;
        }

        vr = func->Func(r) ;
        if ( vr >= vm )
        {
            break ;
        }

        l = m ; vl = vm ;
        m = r ; vm = vr ;
    }

    interval.col(0) = l ;
    interval.col(1) = m ;
    interval.col(2) = r ;
    if (pVal)
    {
        pVal[0] = vl ;
        pVal[1] = vm ;
        pVal[2] = vr ;
    }

    return interval ;
}
