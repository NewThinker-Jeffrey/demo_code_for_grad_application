#ifndef GradDescent_H
#define GradDescent_H

#include "stddef.h"
#include "Optimiz/Optimiz.h"


START_JK_OPTIMIZ

class LineSearch ;
class GradDescent : public Method
{
public:
    struct Solution
    {
        Matrix w ;
        Matrix grad ;
        Float cost ;
        Float lastStepLen ;
        Float lastDiff ;        
        bool valid ;

        Solution() ;
        Solution(const Solution& other) ;
        Solution& operator=(const Solution& other) ;
        void print() ;
    };

    enum EndCondition {
        StepLessThan,
        GradLessThan
    };

public:

    GradDescent ( SingleFunction* cost=NULL, Float eps=0.0000001, int MaxSteps=1000, Float LineSearchEPS=0.0000001, EndCondition cond=StepLessThan, LineSearch* searcher=NULL) ;
    ~GradDescent () ;
    Result solve(const Matrix& startPos) ;


public :
    EndCondition cond ;
    Float eps ;
    int MaxSteps ;
    Float LineSearchEPS ;
    LineSearch* searcher ;

protected :
    bool touchEnd (const Solution& res) ;
    Solution nextStep(const Solution& oldRes) ;
    static Matrix Grad2Dir(const Matrix& grad) ;
    static Float VecLength (const Matrix& singleRow) ;
    Result result ;
};


END_JK_OPTIMIZ


#endif // GradDescent_H
