#include "Kalman/UKF.h"
#include "math.h"

//#include "AccurateTimer.h"
#define TimerPrint if (timePrint) printf

bool UKF_Print_SR = false ;

IMPORT_JK_Math
IMPORT_JK_KALMAN


///////////////////////////////////////////////////////////////////////
///
/// UKF::Model Area
///
///////////////////////////////////////////////////////////////////////
UKF::Model::Model()
{
    curR = &R ;
    curQ = &Q ;
    curPosterior = &Posterior ;
    statusMap = NULL ;
    measureMap = NULL ;
}
UKF::Model::~Model()
{
}

Matrix UKF::Model::SR(const Matrix& A, int cols)
{
    return UKF::LDLT_SR(A, cols) ;
}

bool UKF::Model::setPosterior(const Matrix& X, const Matrix& Pxx)
{
    if (X.rows() != Pxx.rows() ||
        X.rows() != Pxx.cols() ||
        X.rows() == 0 )
    {
        printf("UKF::Model::setPosterior error: dimonsion of the input status and that of the covariance not match, \
               or the covariance matrix is not a sqaure one\n") ;
        return false ;
    }

    curPosterior->Pxx = Pxx.clone() ;
    curPosterior->X = X.clone() ;
    return true ;
}

const PosteriorStatus& UKF::Model::getPosterior() const
{
    return *curPosterior ;
}

Matrix UKF::Model::getStatus() const
{
    if (statusMap)
    {
        return statusMap->map2Manifold_singleCol(Posterior.X) ;
    }
    else
    {
        return Posterior.X ;
    }
}

bool UKF::Model::setStatusMap(ManifoldMap* Map)
{
    if (statusMap != Map)
    {
        if (statusMap)
            delete statusMap ;
        statusMap = Map ;
    }

    return true ;
}

const ManifoldMap* UKF::Model::getStatusMap() const
{
    return statusMap ;
}

bool UKF::Model::setMeasureMap(ManifoldMap* Map)
{
    if (measureMap != Map)
    {
        if (measureMap)
            delete measureMap ;
        measureMap = Map ;
    }

    return true ;
}

const ManifoldMap* UKF::Model::getMeasureMap() const
{
    return measureMap ;
}

Matrix UKF::Model::autoF(const Matrix& tangent_Xk, const Matrix& Uk,
         const Matrix& nonAdditive_Vk)
{
    if (statusMap)
    {
        //printf(" **** tangent_Xk : \n") ;
        //tangent_Xk.transpose().printMat() ;
        return f( statusMap->map2Manifold(tangent_Xk), Uk, nonAdditive_Vk ) ;
    }
    else
        return f( tangent_Xk, Uk, nonAdditive_Vk ) ;
}

Matrix UKF::Model::autoH(const Matrix& tangent_X, const Matrix& nonAdditive_Nk)
{
    if (statusMap)
        return h( statusMap->map2Manifold(tangent_X), nonAdditive_Nk ) ;
    else
        return h(tangent_X, nonAdditive_Nk) ;
}


bool UKF::Model::setR(const Matrix& R, NoiseType type, const Matrix& R2)
{
    curR->type = type ;
    curR->SR_additiveCov = Matrix() ;

    switch(curR->type)
    {
    case Heterogeneous_Noise:
        if (R.rows() != R.cols() || R2.rows() != R2.cols())
        {
            printf("UKF::Model::setR(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curR->additiveCov = R ;
        curR->nonAdditiveCov = R2 ;
        // curR->SR_additiveCov = SR(curR->additiveCov) ; //not necessary.
        curR->SR_nonAdditiveCov = SR(curR->nonAdditiveCov) ;
        break ;

    case NonAdditive_Noise:
        if (R.rows() != R.cols())
        {
            printf("UKF::Model::setR(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curR->nonAdditiveCov = R ;
        curR->SR_nonAdditiveCov = SR(curR->nonAdditiveCov) ;
        curR->additiveCov = Matrix() ;
        break ;

    case PureAdditive_Noise:
        if (R.rows() != R.cols())
        {
            printf("UKF::Model::setR(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curR->additiveCov = R ;
        curR->nonAdditiveCov = Matrix() ;
        curR->SR_nonAdditiveCov = Matrix() ;
        break ;

    default:
        break ;
    }

    return true ;
}

NoiseType UKF::Model::typeR() const
{
    return curR->type ;
}

Matrix UKF::Model::additiveR() const
{
    return curR->additiveCov ;
}

Matrix UKF::Model::nonAdditiveR() const
{
    return curR->nonAdditiveCov ;
}

Matrix UKF::Model::SR_nonAdditiveR() const
{
    return curR->SR_nonAdditiveCov ;
}

int UKF::Model::dimensionAdditiveR() const
{
    return additiveR().rows() ;
}

int UKF::Model::dimensionNonAdditiveR() const
{
    return nonAdditiveR().rows() ;
}


int UKF::Model::dimensionR() const
{
    return additiveR().rows() + nonAdditiveR().rows() ;
}


bool UKF::Model::setQ(const Matrix& Q, NoiseType type,const Matrix& Q2)
{
    curQ->type = type ;
    switch(curQ->type)
    {
    case Heterogeneous_Noise:
        if (Q.rows() != Q.cols() || Q2.rows() != Q2.cols())
        {
            printf("UKF::Model::setQ(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curQ->additiveCov = Q ;
        curQ->SR_additiveCov = SR(curQ->additiveCov) ;
        curQ->nonAdditiveCov = Q2 ;
        curQ->SR_nonAdditiveCov = SR(curQ->nonAdditiveCov) ;
        break ;

    case NonAdditive_Noise:
        if (Q.rows() != Q.cols())
        {
            printf("UKF::Model::setQ(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curQ->nonAdditiveCov = Q ;
        curQ->SR_nonAdditiveCov = SR(curQ->nonAdditiveCov) ;
        curQ->additiveCov = Matrix() ;
        curQ->SR_additiveCov = Matrix();
        break ;

    case PureAdditive_Noise:
        if (Q.rows() != Q.cols())
        {
            printf("UKF::Model::setQ(): the input covariance matrix is not a square one\n") ;
            return false ;
        }
        curQ->additiveCov = Q ;
        curQ->SR_additiveCov = SR(curQ->additiveCov) ;
        curQ->nonAdditiveCov = Matrix() ;
        curQ->SR_nonAdditiveCov = Matrix() ;
        break ;

    default:
        break ;
    }

    return true ;
}

NoiseType UKF::Model::typeQ() const
{
    return curQ->type ;
}

Matrix UKF::Model::additiveQ() const
{
    return curQ->additiveCov ;
}

Matrix UKF::Model::nonAdditiveQ() const
{
    return curQ->nonAdditiveCov ;
}

Matrix UKF::Model::SR_additiveQ() const
{
    return curQ->SR_additiveCov ;
}

Matrix UKF::Model::SR_nonAdditiveQ() const
{
    return curQ->SR_nonAdditiveCov ;
}

int UKF::Model::dimensionAdditiveQ() const
{
    return additiveQ().rows() ;
}

int UKF::Model::dimensionNonAdditiveQ() const
{
    return nonAdditiveQ().rows() ;
}

int UKF::Model::dimensionQ() const
{
    return additiveQ().rows() + nonAdditiveQ().rows() ;
}

///////////////////////////////////////////////////////////////////////
///
/// UKF Static Area
///
///////////////////////////////////////////////////////////////////////

static Float safeSqrt (Float v)
{
    if (v > 0)
        return sqrt(v) ;
    else
        return 0 ;
}

static Matrix ZeroMatrix(int rows, int cols)
{
    if (rows*cols==0)
        return Matrix() ;

    Matrix A(rows, cols) ;
    A.fill(0);
    return A ;
}

Matrix UKF::CholeskySR(const Matrix& A, int cols)
{
    // Let $A=LL^T$, where $L$ is a lower-triangle matrix and
    // $A$ is real-symetrical positive semidefinite matrix ;
    // FIRST we have $A_{11}=L_{11}^2$, which implies
    // $L_{11}=\sqrt{A_{11}}$ ;
    // NEXT consider the whole first column of $L$, it can become
    // that of $A$ by times $L_{11}, so we get all the entries of
    // the first column of $L$ : $L_{i1}=A_{i1}/L_{i1}$ ;
    // THEN if we have already get the first (k-1) column(s) of $L$,
    // given $A_{kk}=\sum_{i=1}^{k}{L_{ki}^2}$, we quickly get
    // $L_{kk}=\sqrt{A_{kk}-\sum_{i=1}^{k-1}{L_{ki}^2}}$ ;
    // Further more, given $A_{ik}=\sum_{j=1}^{k}{L_{ij}L_{kj}}$,
    // so $L_{ik}=(A_{ik}-\sum_{j=1}^{k-1}{L_{ij}L_{kj}})/L_{kk}$ .
    // Following this routine, finally we get all the entries of $L$ .

    int i, j, k;
    int n = A.rows() ;
    if (cols==0)
        cols = n ;

    Float *L = new Float[n*cols] ;
    for (i=0; i<n*cols; i++)
        L[i] = 0 ;

    for (k=0; k<cols; k++)
    {
        Float sum = 0 ;
        Float entry, entry1, entry2 ;
        for (i=0; i<k; i++)
        {
            entry = L[k*cols+i] ;
            sum += entry * entry ;
        }
        L[k*cols+k] = safeSqrt (A(k,k) - sum) ;

        for (i=k+1; i<n; i++)
        {
            sum = 0 ;
            for (j=0; j<k; j++)
            {
                entry1 = L[i*cols+j] ;
                entry2 = L[k*cols+j] ;
                sum += entry1 * entry2 ;
            }
            L[i*cols+k] = (A(i,k)-sum) / L[k*cols+k] ;
        }
    }

    Matrix sL = Matrix(n,cols,L).clone() ;
    delete[] L ;
    return sL ;
}

Matrix UKF::LDLT_SR(const Matrix& A, int cols)
{
    // Let $A=LDL^T$, where $L$ is a unit-lower-triangle matrix,
    // for any $i=1:n$ $L_{ii}=1$ holds, and $D$ is a diagonal
    // matrix with its nonzero entries all greater than 0 , and $A$
    // is real-symetrical positive semidefinite matrix ;
    // we always have $A_{ij}=\sum_{k=1}{j}{L_{ik}D_{kk}L_{jk}}$,
    // and given $L_{jj}=L_{ii}=1$ we can rewrite the equation :
    // $A_{ij}=\sum_{k=1}{j-1}{L_{ik}D_{kk}L_{jk}} + L_{ij}D_{jj}$,
    // for some special cases in which $i=j$, more direct result holds:
    // $A_{jj}=\sum_{k=1}{j-1}{D_{kk}L_{jk}^2} + D_{jj}$ .
    // If we have already known the first (j-1) columns of $L$ and $D$,
    // we can get $D_{jj}$ easily :
    // D_{jj} = $A_{jj} - \sum_{k=1}{j-1}{D_{kk}L_{jk}^2},
    // and then for each $i=(j+1):n$, we have
    // $L_{ij} = ( A_{ij} - \sum_{k=1}{j-1}{L_{ik}D_{kk}L_{jk}} ) / D_{jj}$
    // Following this routine, finally we get all the entries of $L$ and $D$.

    int i, j, k;
    int n = A.rows() ;
    if (cols==0)
        cols=n;

    Float *L = new Float[n*cols] ;
    Float *D = new Float[cols] ;
    for (i=0; i<n*cols; i++)
        L[i] = 0 ;

    for (j=0; j<cols; j++)
    {
        Float sum = 0 ;
        for (k=0; k<j; k++)
        {
            Float entry = L[j*cols+k] ;
            sum += entry * entry * D[k] ;
        }
        D[j] = A(j,j) - sum ;
        L[j*cols+j] = 1 ;

        for (i=j+1; i<n; i++)
        {
            sum = 0 ;
            for (k=0; k<j; k++)
            {
                Float entry1 = L[j*cols+k] ;
                Float entry2 = L[i*cols+k] ;
                sum += entry1 * entry2 * D[k] ;
            }
            L[i*cols+j] = ( A(i,j) - sum ) / D[j] ;
        }
    }


    // $L = L*\sqrt{D}$
    for (i=0; i<cols; i++)
    {
        Float r = 0 ;
        if ( D[i] > 0 )
            r = sqrt(D[i]) ;

        for (j=0; j<n; j++)
        {
            L[j*cols+i] *= r ;
        }
    }
    Matrix sL = Matrix(n,cols,L).clone() ;
    delete[] L ;
    delete[] D ;
    return sL ;
}



// by projecting all the original random samples into a lower-dimension plane, we get
// a new group of samples constrained in the plane. the result of 'ProjectCovX()' is
// the covariance matrix of the projected samples.
// and, the arguement 'OthoSpace' is the complement space of the plane.
Matrix UKF::ProjectCovX (const Matrix& CovX, const Matrix& OthoSpace)
{
    Matrix P = CovX.clone() ;
    Matrix Q ;
    int n = OthoSpace.cols() ;

    // First consider the cases in which 'OthoSpace' is only composed of one unit vector 'q'.
    // supposed $ E[xx^T]=P, qq^T = Q $, target :
    // $$ R = E[(x-(x\cdot q)q)(x-(x\cdot q)q)^T] = \\
    //    E[xx^T] + E[(x\cdot q)^2]qq^T - E[(x\cdot q)xq^T] - E[(x\cdot q)qx^T] =\\
    //    (P + E[(x\cdot q)^2]Q) - (E[(x\cdot q)xq^T] + E[(x\cdot q)qx^T])
    // $$
    // $$ E[(x\cdot q)^2] = E[(\sum_{i}{x_iq_i})(\sum_{j}{x_jq_j})] = \\
    //    E[\sum_{ij}{(x_iq_i)(x_jq_j)}] = E[\sum_{ij}{(x_ix_j)(q_iq_j)}] =\\
    //    \sum_{ij}{E[x_ix_j]q_iq_j} = \sum_{ij}{P_{ij}Q_{ij}} = P \cdot Q
    // $$
    // $$ E[(x\cdot q)xq^T] = E[\sum_k{(x_kq_k)xq^T}] = E[\sum_k{(x_kx)(q_kq^T)}] =\\
    //    \sum_k{(E[x_kx]Q_k^T)} = \sum_k{P_kQ_k^T} = PQ
    // $$
    // $$ E[(x\cdot q)qx^T] = E[\sum_k{(x_kq_k)qx^T}] = E[\sum_k{(q_kq)(x_kx^T)}] =\\
    //    \sum_k{(Q_k E[x_kx^T])}= \sum_k{Q_kP_k^T} = QP
    // $$
    // Therefore, the target is simplified as below :
    // $ R = P + (P\cdot Q)Q - (PQ + QP) $
    // denoted $ S=PQ $, then $ QP=S^T $ and $ P\cdot Q=tr(S) $, substitute the notions
    // in the equation above we finally have :
    // $ R = P + tr(S)Q - (S + S^T) $
    //
    // if 'OthoSpace' is composed of more than one vectors, say 'n' vectors, we then
    // split the 'OthoSpace' into 'n' separate parts that each composed of one vector,
    // and reiterate the process illustrated above 'n' times.

    for (int i=0; i<n; i++)
    {
        Matrix q = OthoSpace.col(i).clone() ;
        Q = q * q.transpose() ;
        Matrix S = P*Q ;
        P += S.trace()*Q - (S + S.transpose()) ;
    }

    return P ;  // it's easy to verified that $q^T R q=0$
}


///////////////////////////////////////////////////////////////////////
///
/// UKF Area
///
///////////////////////////////////////////////////////////////////////


UKF::UKF(Float alpha, Float beta, Float kappa)
{
    this->kappa = kappa ;
    this->alpha = alpha ;
    this->beta  = beta  ;
    this->timePrint = false ;
}

UKF::~UKF()
{
}


bool UKF::update(Model& model, const Matrix& Y, const Matrix& U)
{
    return update(model, Y, U, true) ;
}

bool UKF::predict(Model& model, const Matrix& U)
{
    return update(model, Matrix(), U, false) ;
}

bool UKF::checkDimension(Model& model)
{
    const PosteriorStatus& posterior = model.getPosterior() ;
    if (posterior.X.rows() != posterior.Pxx.rows() ||
        posterior.X.rows() != posterior.Pxx.cols() ||
        posterior.X.rows() == 0 )
    {
        printf("UKF::checkDimension() error: posterior status not initialized, \
               or dimonsion of the posterior status and that of its covariance not match, \
               or the covariance matrix is not a sqaure one\n") ;
        return false ;
    }

    int nX = posterior.X.rows() ;
    switch (model.typeQ())
    {
    case Heterogeneous_Noise:
    case PureAdditive_Noise:
        if ( nX != model.dimensionAdditiveQ() )
        {
            printf("UKF::checkDimension() error: dimension of the posterior\
                   status and the additive process-noise not match \n");
            return false ;
        }
        break ;

    case NonAdditive_Noise:
    default:
        break ;
    }

    return true ;
}


bool UKF::update(Model& model, const Matrix& Ymeasure, const Matrix& U, bool withMeasurement)
{
    if (!checkDimension(model))
        return false ;

    const PosteriorStatus& posterior = model.getPosterior() ;
    int i ;
    int nX = posterior.X.rows() ;
    int nQ = model.dimensionQ() ;
    int nQ_additve = model.dimensionAdditiveQ() ;
    int nQ_nonAdditive = model.dimensionNonAdditiveQ() ;
    int nR = 0 ;
    if (withMeasurement)
        nR = model.dimensionR() ;

    int L = nX + nQ + nR ;
    Float lambda = alpha*alpha*(L+kappa)-L ;
    Float _Wi = 1 / (2*(L+lambda)) ;
    Float _W0m = lambda / (L+lambda) ;
    Float _W0c = _W0m + 1 + beta - (alpha*alpha) ;
    Float r = safeSqrt(L+lambda) ;
    Float Wi, W0m, W0c ;

    Wi = _Wi ;
    W0m = _W0m + 2*Wi*(nQ_additve + nR);
    W0c = _W0c + 2*Wi*(nQ_additve + nR);

    Matrix ZeroV = ZeroMatrix(nQ_nonAdditive,1) ;

    //AccurateTimer timer ;
    //TimerPrint("******[ %d ]*********\n", withMeasurement) ;
    //timer.begin();
    Matrix fX0 = model.autoF(posterior.X, U, ZeroV) ;
    //timer.end();
    //TimerPrint("timerAutoF: %f\n", timer.s()) ;


    Matrix SigmaX (fX0.rows(), 2*(nX+nQ_nonAdditive)+1) ;
    Matrix WeightX(2*(nX+nQ_nonAdditive)+1, 1) ;
    WeightX.fill(Wi);
    WeightX(0,0) = W0m ;
    SigmaX.col(0) = fX0 ;

    // calc SigmaX ;
    //timer.begin();
    Matrix SR_Pxx = model.SR(posterior.Pxx) ;
    //if (UKF_Print_SR)
    //{
    //    printf("SR_Pxx:\n");
    //    SR_Pxx.transpose().printMat() ;
    //    printf("SR_nonAdditiveQ:\n");
    //    model.SR_nonAdditiveQ().transpose().printMat() ;
    //    if (withMeasurement)
    //    {
    //        printf("SR_additiveR:\n");
    //        model.curR->SR_additiveCov.transpose().printMat() ;
    //        printf("SR_additiveR:\n");
    //        model.curR->additiveCov.transpose().printMat() ;
    //    }
    //}

    for (i=0; i<nX; i++)
    {
        Matrix fX1 = model.autoF(posterior.X+r*SR_Pxx.col(i), U, ZeroV) ;
        Matrix fX2 = model.autoF(posterior.X-r*SR_Pxx.col(i), U, ZeroV) ;
        SigmaX.col(2*i+1) = fX1;
        SigmaX.col(2*i+2) = fX2;
        //if (UKF_Print_SR)
        //{
        //    printf("+s%d: \n", i);
        //    (posterior.X+r*SR_Pxx.col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    model.statusMap->map2Manifold(posterior.X+r*SR_Pxx.col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    fX1.transpose().printMat() ;
        //    printf("-s%d: \n", i);
        //    (posterior.X-r*SR_Pxx.col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    model.statusMap->map2Manifold(posterior.X-r*SR_Pxx.col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    fX2.transpose().printMat() ;
        //}

    }

    for (i=0; i<nQ_nonAdditive; i++)
    {
        Matrix fX1 = model.autoF(posterior.X, U,  r*model.SR_nonAdditiveQ().col(i)) ;
        Matrix fX2 = model.autoF(posterior.X, U, -r*model.SR_nonAdditiveQ().col(i)) ;
        SigmaX.col(2*nX+2*i+1) = fX1;
        SigmaX.col(2*nX+2*i+2) = fX2;
        //if (UKF_Print_SR)
        //{
        //    printf("+n%d: \n", i);
        //    (r*model.SR_nonAdditiveQ().col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    fX1.transpose().printMat() ;
        //    printf("-n%d: \n", i);
        //    (-r*model.SR_nonAdditiveQ().col(i)).transpose().printMat() ;
        //    printf("----\n");
        //    fX2.transpose().printMat() ;
        //}
    }

    //if (UKF_Print_SR)
    //{
    //    printf("L: %d, _Wi: %f,   _W0m: %f,  _W0c: %f, r: %f\n", L, _Wi, _W0m, _W0c, r) ;
    //    printf("posterior.X:\n");
    //    posterior.X.transpose().printMat() ;
    //    printf("Control:\n");
    //    U.transpose().printMat() ;
    //}

    //timer.end();
    //TimerPrint("calc SigmaX: %f\n", timer.s()) ;


    //timer.begin();
    Matrix SigmaX_Tangent = SigmaX ;
    if (model.statusMap)
    {
        Matrix Xref  = model.statusMap->calcBestRefPoint(SigmaX, WeightX);
        model.statusMap->setRefPoint( Xref ) ;
        SigmaX_Tangent = model.statusMap->map2Tangent(SigmaX) ;
    }
    //timer.end();
    //TimerPrint("calc SigmaX_Tangent: %f\n", timer.s()) ;


    // calc Xmean and Xcov ;
    //timer.begin();
    Matrix Xmean = SigmaX_Tangent * WeightX ;
    Matrix Xcov = W0c * (SigmaX_Tangent.col(0) - Xmean) *
            (SigmaX_Tangent.col(0) - Xmean).transpose() ;
    for (i=1; i<SigmaX_Tangent.cols(); i++)
    {
        Matrix X_Tangent = SigmaX_Tangent.col(i) ;
        Xcov += Wi * (X_Tangent-Xmean) * (X_Tangent-Xmean).transpose() ;
    }    
    if (nQ_additve > 0)
    {
        Xcov += model.additiveQ() ;
    }
    //timer.end();
    //TimerPrint("calc Xmean and Xcov: %f\n", timer.s()) ;

    //if (UKF_Print_SR)
    //{
    //    printf("----\n");
    //    printf("SigmaX_Tangent: \n");
    //    SigmaX_Tangent.printMat() ;
    //    printf("----\n");
    //    printf("Xcov: \n");
    //    Xcov.printMat() ;
    //    printf("----\n");
    //    printf("Xmean: \n");
    //    Xmean.transpose().printMat() ;
    //    printf("----\n");
    //}


    if (!withMeasurement)
    {
        //if (UKF_Print_SR)
        //{
        //    UKF_Print_SR = false ;
        //}
        return model.setPosterior(Xmean, Xcov) ;
    }


    // measurement env .
    int nR_additve = model.dimensionAdditiveR() ;
    int nR_nonAdditive = model.dimensionNonAdditiveR() ;

    Wi = _Wi ;
    W0m = _W0m + 2*Wi*(nR_additve);
    W0c = _W0c + 2*Wi*(nR_additve);

    Matrix ZeroN = ZeroMatrix(nR_nonAdditive,1) ;

    //timer.begin();
    Matrix Y0 = model.h(fX0, ZeroN) ;
    //timer.end();
    //TimerPrint("model.h(): %f\n", timer.s()) ;


    Matrix SigmaY(Y0.rows(), 2*(nX+nQ+nR_nonAdditive)+1) ;
    Matrix WeightY(2*(nX+nQ+nR_nonAdditive)+1, 1);
    WeightY.fill(Wi);
    WeightY(0,0) = W0m ;
    SigmaY.col(0) = Y0 ;

    // calc SigmaY ;
    //timer.begin();
    for (i=1; i<SigmaX.cols(); i++)
    {
        Matrix fX = SigmaX.col(i).clone() ;
        Matrix Y = model.h(fX, ZeroN) ;
        SigmaY.col(i) = Y ;
    }
    for (i=0; i<model.dimensionAdditiveQ(); i++)
    {
        Matrix fX_Tangnt, Y;

        fX_Tangnt = SigmaX_Tangent.col(0) + r * model.SR_additiveQ().col(i) ;
        Y = model.autoH(fX_Tangnt, ZeroN) ;
        SigmaY.col(SigmaX.cols()+2*i) = Y ;

        fX_Tangnt = SigmaX_Tangent.col(0) - r * model.SR_additiveQ().col(i) ;
        Y = model.autoH(fX_Tangnt, ZeroN) ;
        SigmaY.col(SigmaX.cols()+2*i+1) = Y ;
    }
    for (i=0; i<model.dimensionNonAdditiveR(); i++)
    {
        Matrix Y, Nk ;
        Nk = r * model.SR_nonAdditiveR().col(i) ;

        Y = model.h(fX0, Nk) ;
        SigmaY.col(SigmaX.cols()+2*model.dimensionAdditiveQ()+2*i) = Y ;
        Y = model.h(fX0, -Nk) ;
        SigmaY.col(SigmaX.cols()+2*model.dimensionAdditiveQ()+2*i+1) = Y ;
    }
    //timer.end();
    //TimerPrint("calc SigmaY: %f\n", timer.s()) ;


    //timer.begin();
    Matrix SigmaY_Tangent = SigmaY ;
    if (model.measureMap)
    {
        Matrix Yref  = model.measureMap->calcBestRefPoint(SigmaY, WeightY);
        //Matrix Yref  = Ymeasure ;
        model.measureMap->setRefPoint( Yref ) ;
        SigmaY_Tangent = model.measureMap->map2Tangent(SigmaY) ;
    }
    //timer.end();
    //TimerPrint("calc SigmaY_Tangent: %f\n", timer.s()) ;


    // calc Ymean, Ycov and Pxy ;
    //timer.begin();
    Matrix Ymean = SigmaY_Tangent * WeightY ;
    Matrix Ycov = W0c * (SigmaY_Tangent.col(0) - Ymean) *
            (SigmaY_Tangent.col(0) - Ymean).transpose() ;
    Matrix Pxy  = W0c * (SigmaX_Tangent.col(0)- Xmean) *
            (SigmaY_Tangent.col(0) - Ymean).transpose() ;

    for (i=1; i<SigmaY_Tangent.cols(); i++)
    {
        Matrix Y = SigmaY_Tangent.col(i) ;
        Matrix X ;

        int j = i - SigmaX_Tangent.cols() ;
        if (j<0)
        {
            X = SigmaX_Tangent.col(i).clone() ;
        }
        else if (j<2*model.dimensionAdditiveQ())
        {
            if ((j%2) == 0)
            {
                X = SigmaX_Tangent.col(0) + r * model.SR_additiveQ().col(j/2) ;
            }
            else
            {
                X = SigmaX_Tangent.col(0) - r * model.SR_additiveQ().col(j/2) ;
            }
        }
        else
        {
            X = SigmaX_Tangent.col(0) ;
        }

        Ycov  += Wi * (Y - Ymean) * (Y - Ymean).transpose() ;
        Pxy   += Wi * (X - Xmean) * (Y - Ymean).transpose() ;
    }

    if (nR_additve > 0)
    {
        Ycov += model.additiveR() ;
    }

    //if (UKF_Print_SR)
    //{
    //    printf("----\n");
    //    printf("SigmaY_Tangent: \n");
    //    SigmaY_Tangent.printMat() ;
    //    printf("----\n");
    //    printf("Pxy: \n");
    //    Pxy.printMat() ;
    //    printf("----\n");
    //    printf("Ycov: \n");
    //    Ycov.printMat() ;
    //    printf("----\n");
    //    printf("Ymean: \n");
    //    Ymean.transpose().printMat() ;
    //    printf("----\n");
    //    printf("Ymeasure: \n");
    //    model.measureMap->map2Tangent(Ymeasure).transpose().printMat() ;
    //    Ymeasure.transpose().printMat() ;
    //}


    //timer.end();
    //TimerPrint("calc Ymean, Ycov and Pxy: %f\n", timer.s()) ;

    // clac Kalman Gain Matrix and update status .
    //timer.begin();
    Matrix Ycov_Inv = Ycov.getInverse() ;
    //timer.end();
    //TimerPrint("Ycov.getInverse(): %f\n", timer.s()) ;

    //if (Ycov_Inv.isEmpty())
    //{
    //    printf("----------\n") ;
    //    Ycov.printMat() ;
    //}

    //timer.begin();
    Matrix K = Pxy * Ycov_Inv ;
    Matrix postX = Xmean + K * ( model.measureMap->map2Tangent(Ymeasure) - Ymean) ;
    Matrix postPxx = Xcov - K * Ycov * K.transpose() ;
    //timer.end();
    //TimerPrint("Kalman update: %f\n", timer.s()) ;


    //if (UKF_Print_SR)
    //{
    //    printf("----\n");
    //    printf("Ycov_Inv: \n");
    //    Ycov_Inv.printMat() ;
    //    printf("----\n");
    //    printf("Kalman: \n");
    //    K.printMat() ;
    //    printf("----\n");
    //    printf("postPxx: \n");
    //    postPxx.printMat() ;
    //    printf("----\n");
    //    printf("postX: \n");
    //    postX.transpose().printMat() ;
    //    printf("----\n");
    //    UKF_Print_SR = false ;
    //}

    return model.setPosterior(postX, postPxx) ;
}
