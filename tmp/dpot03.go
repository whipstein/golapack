package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dpot03 computes the residual for a symmetric matrix times its
// inverse:
//    norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps),
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
//       SUBROUTinE Dpot03( uplo, n, a, lda, ainv, ldainv, work, ldwork,
//                          rwork, rcond, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            lda, ldainv, ldwork, N
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
// Dpot03 computes the residual for a symmetric matrix times its
// inverse:
//    norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the upper or lower triangular part of the
//          symmetric matrix A is stored:
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original symmetric matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,N)
// \endverbatim
//
// \param[in,out] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (ldainv,N)
//          On entry, the inverse of the matrix a, stored as a symmetric
//          matrix in the same format as A.
//          In this version, ainv is expanded into a full matrix and
//          multiplied by a, so the opposing triangle of ainv will be
//          changed; i.e., if the upper triangular part of ainv is
//          stored, the lower triangular part will be used as work space.
// \endverbatim
//
// \param[in] ldainv
// \verbatim
//          ldainv is intEGER
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
//          ldwork is intEGER
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
//          norm(I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)
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
func Dpot03(uplo *byte, n *int, a *[][]float64, lda *int, ainv *[][]float64, ldainv *int, work *[][]float64, ldwork *int, rwork *[]float64, rcond *float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
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
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (lda), (rwork)))
	(*ainvnm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (ainv), (ldainv), (rwork)))
	if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
		(*(rcond)) = (*zero)
		(*(resid)) = (*one) / (*eps)
		return
	}
	(*(rcond)) = ((*one) / (*anorm)) / (*ainvnm)
	//
	//     Expand ainv into a full matrix and call Dsymm to multiply
	//     ainv on the left by A.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
				(*(ainv))[(*j)-1][(*i)-1] = (*(ainv))[(*i)-1][(*j)-1]
				//Label10:
			}
			//Label20:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
				(*(ainv))[(*j)-1][(*i)-1] = (*(ainv))[(*i)-1][(*j)-1]
				//Label30:
			}
			//Label40:
		}
	}
	Dsymm(func() *[]byte {y := []byte("Left"); return &y }(), (uplo), (n), (n), -(*one), (a), (lda), (ainv), (ldainv), zero, (work), (ldwork))
	//
	//     Add the identity matrix to work .
	//
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		(*(work))[(*i)-1][(*i)-1] = (*(work))[(*i)-1][(*i)-1] + (*one)
		//Label50:
	}
	//
	//     Compute norm(I - A*ainv) / (n * norm(a) * norm(ainv) * eps)
	//
	(*(resid)) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (n), (work), (ldwork), (rwork)))
	//
	(*(resid)) = (((*(resid)) * (*(rcond))) / (*eps)) / DBLE((*(n)))
	//
	return
	//
	//     End of Dpot03
	//
}
