package goblas

import 

// Dget03 computes the residual for a general matrix times its inverse:
//    norm( I - ainv*A) / ( N * norm(a) * norm(ainv) * eps),
// where eps is the machine epsilon.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dget03( n, a, lda, ainv, ldainv, work, ldwork, rwork,
//                          rcond, resid)
//
//       .. Scalar Arguments ..
//       inTEGER            lda, ldainv, ldwork, N
//       DOUBLE PRECISION   rcond, resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), ainv( ldainv, *), rwork(*),
//      $                   work( ldwork, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dget03 computes the residual for a general matrix times its inverse:
//    norm( I - ainv*A) / ( N * norm(a) * norm(ainv) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original N x N matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.  lda >= max(1,N).
// \endverbatim
//
// \param[in] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (ldainv,N)
//          The inverse of the matrix A.
// \endverbatim
//
// \param[in] ldainv
// \verbatim
//          ldainv is inTEGER
//          The leading dimension of the array ainv.  ldainv >= max(1,N).
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (ldwork,N)
// \endverbatim
//
// \param[in] ldwork
// \verbatim
//          ldwork is inTEGER
//          The leading dimension of the array work.  ldwork >= max(1,N).
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The reciprocal of the condition number of a, computed as
//          ( 1/norm(a)) / norm(ainv).
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(I - ainv*A) / ( N * norm(a) * norm(ainv) * eps)
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG Ltd.
//
// \date December 2016
//
// \ingroup double_lin
//
//  =====================================================================
func Dget03(n *int, a *[][]float64, lda *int, ainv *[][]float64, ldainv *int, work *[][]float64, ldwork *int, rwork *[]float64, rcond *float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	eps := new(float64)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Quick exit if N = 0.
	//
	if (*(n)) <= 0 {
		(*(rcond)) = (*one)
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0 or ainvnm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), (n), (a), (lda), (rwork)))
	(*ainvnm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), (n), (ainv), (ldainv), (rwork)))
	if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
		(*(rcond)) = (*zero)
		(*(resid)) = (*one) / (*eps)
		return
	}
	(*(rcond)) = ((*one) / (*anorm)) / (*ainvnm)
	//
	//     Compute I - A * ainv
	//
	Dgemm(func() *[]byte {y :=[]byte("No transpose"); return &y}(), func() *[]byte {y :=[]byte("No transpose"); return &y}(), (n), (n), (n), -(*one), (ainv), (ldainv), (a), (lda), zero, (work), (ldwork))
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		(*(work))[(*i)-1][(*i)-1] = (*one) + (*(work))[(*i)-1][(*i)-1]
		//Label10:
	}
	//
	//     Compute norm(I - ainv*A) / (n * norm(a) * norm(ainv) * eps)
	//
	(*(resid)) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), (n), (work), (ldwork), (rwork)))
	//
	(*(resid)) = (((*(resid)) * (*(rcond))) / (*eps)) / DBLE((*(n)))
	//
	return
	//
	//     End of Dget03
	//
}
