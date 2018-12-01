
#include "Kalman/Kalman.h"

IMPORT_JK_Math
IMPORT_JK_KALMAN

///////////////////////////////////////////////////////////////////////
///
/// PriorStatus Area
///
///////////////////////////////////////////////////////////////////////

PriorStatus::PriorStatus()
{

}

PriorStatus::PriorStatus(const PriorStatus& other)
{
    Xpre = other.Xpre ;
    Ypre = other.Ypre ;
    Pxx = other.Pxx ;
    Pxy = other.Pxy ;
    Pyy = other.Pyy ;
}

PriorStatus& PriorStatus::operator=(const PriorStatus& other)
{
    Xpre = other.Xpre ;
    Ypre = other.Ypre ;
    Pxx = other.Pxx ;
    Pxy = other.Pxy ;
    Pyy = other.Pyy ;
    return *this ;
}



///////////////////////////////////////////////////////////////////////
///
/// PosteriorStatus Area
///
///////////////////////////////////////////////////////////////////////

PosteriorStatus::PosteriorStatus()
{

}

PosteriorStatus::PosteriorStatus(const PosteriorStatus& other)
{
    X = other.X ;
    Pxx = other.Pxx ;
}

PosteriorStatus& PosteriorStatus::operator=(const PosteriorStatus& other)
{
    X = other.X ;
    Pxx = other.Pxx ;
    return *this ;
}


///////////////////////////////////////////////////////////////////////
///
/// ManifoldMap Area
///
///////////////////////////////////////////////////////////////////////

ManifoldMap::~ManifoldMap()
{

}

void ManifoldMap::setRefPoint(const Matrix& refPointOnManifold)
{

}

Matrix ManifoldMap::calcBestRefPoint (const Matrix& refPointOnManifold, const Matrix& Weight) const
{
    return project2Manifold_singleCol ( refPointOnManifold * Weight ) ;
}


Matrix ManifoldMap::project2Manifold_singleCol(const Matrix& pointsOutsideManifold) const
{
    return pointsOutsideManifold;
}

Matrix ManifoldMap::map2Manifold_singleCol(const Matrix& pointsInTangentSpace) const
{
    return pointsInTangentSpace ;
}

Matrix ManifoldMap::map2Tangent_singleCol(const Matrix& pointsOnManifold) const
{
    return pointsOnManifold ;
}


Matrix ManifoldMap::project2Manifold(const Matrix& pointsOutsideManifold) const
{
    int n = pointsOutsideManifold.cols() ;
    Matrix P0 = project2Manifold_singleCol(pointsOutsideManifold.col(0)) ;
    Matrix ret(P0.rows(), n) ;
    ret.col(0) = P0 ;
    for (int i=1; i<n; i++)
    {
        ret.col(i) = project2Manifold_singleCol(pointsOutsideManifold.col(i)) ;
    }
    return ret ;
}

Matrix ManifoldMap::map2Manifold(const Matrix& pointsInTangentSpace) const
{
    int n = pointsInTangentSpace.cols() ;
    Matrix P0 = map2Manifold_singleCol(pointsInTangentSpace.col(0)) ;
    Matrix ret(P0.rows(), n) ;
    ret.col(0) = P0 ;
    for (int i=1; i<n; i++)
    {
        ret.col(i) = map2Manifold_singleCol(pointsInTangentSpace.col(i)) ;
    }
    return ret ;
}

Matrix ManifoldMap::map2Tangent(const Matrix& pointsOnManifold) const
{
    int n = pointsOnManifold.cols() ;
    Matrix P0 = map2Tangent_singleCol(pointsOnManifold.col(0)) ;
    Matrix ret(P0.rows(), n) ;
    ret.col(0) = P0 ;
    for (int i=1; i<n; i++)
    {
        ret.col(i) = map2Tangent_singleCol(pointsOnManifold.col(i)) ;
    }
    return ret ;
}




