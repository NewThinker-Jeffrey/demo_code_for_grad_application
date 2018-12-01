#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "Optimiz/Optimiz.h"


START_JK_OPTIMIZ

class LineSearch
{
public:
    LineSearch(){}
    virtual ~LineSearch(){}
    virtual Matrix move (const Matrix& start, const Matrix& dir, Float stepLen) ;
    virtual Float distance2 (const Matrix& w1, const Matrix& w2) ;
    virtual Matrix interpolate (const Matrix& start, const Matrix& end, Float endWeight) ;
    virtual Matrix middleOf (const Matrix& start, const Matrix& end) ;
    virtual bool isStepAcceptable(const Matrix& start, const Matrix& end) ;
    virtual Matrix operator() (SingleFunction* func,
              const Matrix& dir, const Matrix& start,
              Float startStepLen=1, Float eps = 0.000001) ;


protected:
    Matrix findInterval (SingleFunction* func, const Matrix& dir, const Matrix& start, Float startStepLen, Float* val=NULL) ;
};


END_JK_OPTIMIZ

#endif // LINESEARCH_H
