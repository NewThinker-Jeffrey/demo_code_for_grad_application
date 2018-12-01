#ifndef CONJGRADIENT_H
#define CONJGRADIENT_H

#include "Optimiz/Optimiz.h"

START_JK_OPTIMIZ

class LineSearch ;
class ConjGradient : public Method
{
public:
    ConjGradient(SingleFunction* cost=NULL,
           Float _gradEps=0.00001,
           Float lineSearchEps = 0.00001,
           int maxLoops = 1000 );
    virtual ~ConjGradient() ;
    virtual Result solve(const Matrix &startPoint) ;

    Float gradEps ;
    Float lineSearchEps ;
    int maxLoops ;
    LineSearch* searcher ;
};


END_JK_OPTIMIZ

#endif // CONJGRADIENT_H
