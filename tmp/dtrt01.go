package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtrt01 computes the residual for a triangular matrix A times its
// inverse:
//    resid = norm( A*ainv - I) / ( N * norm(a) * norm(ainv) * eps),
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
//       SUBROUTinE Dtrt01( uplo, diag, n, a, lda, ainv, ldainv, rcond,
//                          work, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, uplo
//       intEGER            lda, ldainv, N
//       DOUBLE PRECISION   rcond, resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), ainv( ldainv, *), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtrt01 computes the residual for a triangular matrix A times its
// inverse:
//    resid = norm( A*ainv - I) / ( N * norm(a) * norm(ainv) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the matrix A is upper or lower triangular.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//          Specifies whether or not the matrix A is unit triangular.
//          = 'N':  Non-unit triangular
//          = 'U':  Unit triangular
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The triangular matrix A.  If uplo = 'U', the leading n by n
//          upper triangular part of the array A contains the upper
//          triangular matrix, and the strictly lower triangular part of
//          A is not referenced.  If uplo = 'L', the leading n by n lower
//          triangular part of the array A contains the lower triangular
//          matrix, and the strictly upper triangular part of A is not
//          referenced.  If diag = 'U', the diagonal elements of A are
//          also not referenced and are assumed to be 1.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,N).
// \endverbatim
//
// \param[in,out] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (ldainv,N)
//          On entry, the (triangular) inverse of the matrix a, in the
//          same storage format as A.
//          On exit, the contents of ainv are destroyed.
// \endverbatim
//
// \param[in] ldainv
// \verbatim
//          ldainv is intEGER
//          The leading dimension of the array ainv.  ldainv >= max(1,N).
// \endverbatim
//
// \param[out] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The reciprocal condition number of a, computed as
//          1/(norm(a) * norm(ainv)).
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(A*ainv - I) / ( N * norm(a) * norm(ainv) * eps)
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
func Dtrt01(uplo *byte, diag *byte, n *int, a *[][]float64, lda *int, ainv *[][]float64, ldainv *int, rcond *float64, work *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
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
	//     Quick exit if N = 0
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
	(*anorm) = (*Dlantr(func() *byte {y := byte('1'); return &y }(), (uplo), (diag), (n), (n), (a), (lda), (work)))
	(*ainvnm) = (*Dlantr(func() *byte {y := byte('1'); return &y }(), (uplo), (diag), (n), (n), (ainv), (ldainv), (work)))
	if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
		(*(rcond)) = (*zero)
		(*(resid)) = (*one) / (*eps)
		return
	}
	(*(rcond)) = ((*one) / (*anorm)) / (*ainvnm)
	//
	//     Set the diagonal of ainv to 1 if ainv has unit diagonal.
	//
	if blas.Lsame((diag), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			(*(ainv))[(*j)-(1)][(*j)-(1)] = (*one)
			//Label10:
		}
	}
	//
	//     Compute A * ainv, overwriting ainv.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dtrmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("No transpose"); return &y }(), (diag), j, (a), (lda), &((*(ainv))[0][(*j)-(1)]), func() *int {y := 1; return &y }())
			//Label20:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dtrmv(func() *[]byte {y := []byte("Lower"); return &y }(), func() *[]byte {y := []byte("No transpose"); return &y }(), (diag), (*(n))-(*j)+1, &((*(a))[(*j)-(1)][(*j)-(1)]), (lda), &((*(ainv))[(*j)-(1)][(*j)-(1)]), func() *int {y := 1; return &y }())
			//Label30:
		}
	}
	//
	//     Subtract 1 from each diagonal element to form A*ainv - I.
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		(*(ainv))[(*j)-(1)][(*j)-(1)] = (*(ainv))[(*j)-(1)][(*j)-(1)] - (*one)
		//Label40:
	}
	//
	//     Compute norm(A*ainv - I) / (n * norm(a) * norm(ainv) * eps)
	//
	(*(resid)) = (*Dlantr(func() *byte {y := byte('1'); return &y }(), (uplo), func() *[]byte {y := []byte("Non-unit"); return &y }(), (n), (n), (ainv), (ldainv), (work)))
	//
	(*(resid)) = (((*(resid)) * (*(rcond))) / DBLE((*(n)))) / (*eps)
	//
	return
	//
	//     End of Dtrt01
	//
}
