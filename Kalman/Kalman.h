#ifndef KALMAN_H
#define KALMAN_H

#include "MathToolsCpp/Matrix.h"

#define START_JK_KALMAN       namespace JK { namespace Math { namespace Kalman {
#define END_JK_KALMAN                          }}}
#define IMPORT_JK_KALMAN      using namespace JK::Math::Kalman ;
#define JKKalman   JK::Math::Kalman

START_JK_KALMAN


struct PriorStatus
{
    Matrix Xpre ;
    Matrix Ypre ;
    Matrix Pxx ;
    Matrix Pxy ;
    Matrix Pyy ;

    PriorStatus() ;
    PriorStatus(const PriorStatus& other) ;
    PriorStatus& operator=(const PriorStatus& other) ;
};

struct PosteriorStatus
{
    Matrix X ;
    Matrix Pxx ;

    PosteriorStatus() ;
    PosteriorStatus(const PosteriorStatus& other) ;
    PosteriorStatus& operator=(const PosteriorStatus& other) ;
};



struct ManifoldMap
{
    virtual ~ManifoldMap() ;
    virtual void setRefPoint(const Matrix& refPointOnManifold) ;
    virtual Matrix calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const ;
    virtual Matrix project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const ;
    virtual Matrix map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const ;
    virtual Matrix map2Tangent_singleCol(const Matrix& pointsOnManifold) const ;

    Matrix project2Manifold(const Matrix& pointsOutsideManifold) const ;
    Matrix map2Manifold(const Matrix& pointsInTangentSpace) const ;
    Matrix map2Tangent(const Matrix& pointsOnManifold) const ;
};


END_JK_KALMAN



#endif // KALMAN_H
