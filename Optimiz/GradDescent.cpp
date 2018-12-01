#include "Optimiz/GradDescent.h"
#include "Optimiz/LineSearch.h"
#include "cmath"
#include "float.h"

IMPORT_JK_OPTIMIZ
IMPORT_JK_Math

///////////////////////////////////////////////////////////////////////////
///
///

const Float DefaultMaxStepLen = 1 ;

GradDescent::GradDescent(SingleFunction* cost, Float eps,
        int MaxSteps, Float LineSearchEPS,
        EndCondition cond,LineSearch* searcher )
    : Method(cost)
{
    this->MaxSteps = MaxSteps ;
    this->eps = eps ;
    this->cond = cond ;
    this->LineSearchEPS = LineSearchEPS ;
    this->searcher = searcher ;
}


GradDescent::~GradDescent()
{

}


Result GradDescent::solve(const Matrix& startPos)
{
    result.valid = false ;
    if ( cost == NULL )
        return result ;

    Solution res ;
    res.w = startPos ;
    res.cost = cost->Func(res.w) ;
    res.grad = cost->Grad(res.w) ;
    res.lastDiff = FLT_MAX ;
    res.lastStepLen = DefaultMaxStepLen ;
    res.valid = false ;
    //res.print();

    int steps = 0 ;
    while (1)
    {

        //res.print();
        if ( touchEnd(res) || steps >= MaxSteps)
        {
            res.valid = true ;
            break ;
        }

        steps ++ ;        
        if ( searcher )
        {
            LineSearch& search = *searcher ;
            Matrix w = search(cost, -res.grad.transpose(), res.w, res.lastStepLen, LineSearchEPS) ;
            Float c = cost->Func(w) ;
            res.lastStepLen = sqrt(search.distance2(w,res.w)) ;
            res.lastDiff = res.cost - c ;
            res.w = w ;
            res.cost = c ;
            res.grad = cost->Grad(w) ;
            res.valid = true ;
        }
        else
        {
            res = nextStep(res) ;
        }
    }

    printf("GradDescent::MoveSteps  %d\n", steps ) ;
    //res.print();
    result.cost = res.cost ;
    result.grad = res.grad ;
    result.valid = res.valid ;
    result.w = res.w ;
    return result ;
}

bool GradDescent::touchEnd (const Solution& res)
{    
    if ( cond == StepLessThan && res.lastStepLen < eps )
        return true ;
    else if ( cond == GradLessThan && res.grad.squareSum() < eps*eps )
        return true ;
    else
        return false ;
}

GradDescent::Solution GradDescent::nextStep(const Solution& oldRes)
{
    Matrix curPos = oldRes.w ;
    Matrix curGrad = oldRes.grad ;
    Float curCost = oldRes.cost ;
    Matrix curDir = Grad2Dir(curGrad) ;
    Matrix nextPos ;
    Float  nextCost ;
    Solution res ;
    if ( VecLength(curDir) <= FLT_MIN )
    {
        // the length of curDir should be 1 in normal, unless the
        // length of curGrad vanish, which imply the terminate of iteration.
        res = oldRes ;
        res.lastDiff = 0 ;
        res.lastStepLen = 0 ;
        res.valid = true ;
        return res ;
    }

    // we start from trying a bigger stepLen than the last time,
    // to speed up the convergency if possible .
    Float len = oldRes.lastStepLen * 2.0 ;
    len = len > DefaultMaxStepLen ? DefaultMaxStepLen : len ;
    nextPos = curPos + len*curDir ;
    nextCost = cost->Func(nextPos) ;


    while ( nextCost > curCost )
    {
        if ( len < LineSearchEPS )
        {
            // sometimes when the real gradient is very small, the calculated
            // gradient we hold may has a reversed sign. e.g. the real gradient
            // is 0.0000001, but the floating point error may lead us to get a
            // negative gradient by calculation, say -0.0000001, then the grown
            // direction will be mistaken(reversed), which make that even the
            // 'len' here becomes very very small, the 'nextCost' is keeping
            // bigger than 'curCost' . in these cases, we consider the old
            // Solution as the best because the gradient there must be almost 0.
            res = oldRes ;
            res.lastStepLen = 0 ;
            res.lastDiff = 0 ;
            res.valid = true ;
            return res ;
        }

        len *= 0.5 ;
        nextPos = curPos + len*curDir ;
        nextCost = cost->Func(nextPos) ;
    }

    Matrix testPos ;
    Float testCost ;
    Float testLen = len * 0.5 ;
    testPos = curPos + testLen ;
    testCost = cost->Func(testPos) ;

    while ( testCost <= nextCost )
    {
        len = testLen ;
        nextPos = testPos ;
        nextCost = testCost ;
        testLen *= 0.5 ;
        testPos = curPos + testLen ;
        testCost = cost->Func(testPos) ;
    }

    res.lastStepLen = len ;
    res.cost = nextCost ;
    res.w = nextPos ;
    res.grad = cost->Grad(res.w) ;
    res.lastDiff = curCost - nextCost ;
    res.valid = true;
    return res ;
}



///////////////////////////////////////////////////////////////////////////
///
///

Float GradDescent::VecLength (const Matrix& singleRow)
{
    return sqrt(singleRow.squareSum()) ;
}

Matrix GradDescent::Grad2Dir(const Matrix& grad)
{
    Float len = VecLength(grad) ;
    if ( len > FLT_MIN )
    {
        Float r = -1.0/len ;
        return grad.transpose()*r ;
    }
    else
    {
        Matrix zeroDir = grad.transpose().clone() ;
        zeroDir.fill(0);
        return zeroDir ;
    }
}

///////////////////////////////////////////////////////////////////////////
///
///
GradDescent::Solution::Solution()
{
    valid = false ;
}

GradDescent::Solution::Solution(const Solution& other)
{
    valid        =  other.valid ;
    w            =  other.w            ;
    grad         =  other.grad         ;
    cost         =  other.cost         ;
    lastStepLen  =  other.lastStepLen  ;
    lastDiff     =  other.lastDiff     ;
}

GradDescent::Solution& GradDescent::Solution::operator=(const Solution& other)
{
    valid        =  other.valid ;
    w            =  other.w            ;
    grad         =  other.grad         ;
    cost         =  other.cost         ;
    lastStepLen  =  other.lastStepLen  ;
    lastDiff     =  other.lastDiff     ;
    return *this ;
}

void GradDescent::Solution::print()
{
    printf("*****  GradDescentSolution  *****\n");
    printf("new W: ");
    w.transpose().printMat() ;
    printf("new Grad: ");
    grad.printMat() ;
    printf("new GradLen: %7.8f\n", VecLength(grad));
    printf("direction: ");
    Grad2Dir(grad).transpose().printMat() ;
    printf("new cost: %7.8f\n", cost) ;
    printf("stepLen and Diff: %7.8f, %7.8f\n", lastStepLen, lastDiff) ;
}

