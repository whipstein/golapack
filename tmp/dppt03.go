package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dppt03 computes the residual for a symmetric packed matrix times its
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
//       SUBROUTinE Dppt03( uplo, n, a, ainv, work, ldwork, rwork, rcond,
//                          resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            ldwork, N
//       DOUBLE PRECISION   rcond, resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a(*), ainv(*), rwork(*),
//      $                   work( ldwork, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dppt03 computes the residual for a symmetric packed matrix times its
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
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The original symmetric matrix a, stored as a packed
//          triangular matrix.
// \endverbatim
//
// \param[in] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The (symmetric) inverse of the matrix a, stored as a packed
//          triangular matrix.
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
func Dppt03(uplo *byte, n *int, a *[]float64, ainv *[]float64, work *[][]float64, ldwork *int, rwork *[]float64, rcond *float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
	jj := new(int)
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
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
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
	(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (rwork)))
	(*ainvnm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (ainv), (rwork)))
	if (*anorm) <= (*zero) || (*ainvnm) == (*zero) {
		(*(rcond)) = (*zero)
		(*(resid)) = (*one) / (*eps)
		return
	}
	(*(rcond)) = ((*one) / (*anorm)) / (*ainvnm)
	//
	//     uplo = 'U':
	//     Copy the leading N-1 x N-1 submatrix of ainv to work(1:N,2:N) and
	//     expand it to a full matrix, then multiply by A one column at a
	//     time, moving the result one column to the left.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		//
		//        Copy ainv
		//
		(*jj) = 1
		for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
			Dcopy(J, &((*(ainv))[(*jj)-1]), func() *int {y := 1; return &y }(), &((*(work))[0][(*j)+0]), func() *int {y := 1; return &y }())
			Dcopy((*j)-1, &((*(ainv))[(*jj)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*j)-1][1]), (ldwork))
			(*jj) = (*jj) + (*j)
			//Label10:
		}
		(*jj) = (((*(n))-1)*(*(n)))/2 + 1
		Dcopy((*(n))-1, &((*(ainv))[(*jj)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*(n))-1][1]), (ldwork))
		//
		//        Multiply by A
		//
		for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
			Dspmv(func() *[]byte {y := []byte("Upper"); return &y }(), (n), -(*one), (a), &((*(work))[0][(*j)+0]), func() *int {y := 1; return &y }(), zero, &((*(work))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label20:
		}
		Dspmv(func() *[]byte {y := []byte("Upper"); return &y }(), (n), -(*one), (a), &((*(ainv))[(*jj)-1]), func() *int {y := 1; return &y }(), zero, &((*(work))[0][(*(n))-1]), func() *int {y := 1; return &y }())
		//
		//     uplo = 'L':
		//     Copy the trailing N-1 x N-1 submatrix of ainv to work(1:N,1:N-1)
		//     and multiply by a, moving each column to the right.
		//
	} else {
		//
		//        Copy ainv
		//
		Dcopy((*(n))-1, &((*(ainv))[1]), func() *int {y := 1; return &y }(), &((*(work))[0][0]), (ldwork))
		(*jj) = (*(n)) + 1
		for (*j) = 2; (*j) <= (*(n)); (*j)++ {
			Dcopy((*(n))-(*j)+1, &((*(ainv))[(*jj)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*j)-1][(*j)-0]), func() *int {y := 1; return &y }())
			Dcopy((*(n))-(*j), &((*(ainv))[(*jj)+0]), func() *int {y := 1; return &y }(), &((*(work))[(*j)-1][(*j)-1]), (ldwork))
			(*jj) = (*jj) + (*(n)) - (*j) + 1
			//Label30:
		}
		//
		//        Multiply by A
		//
		for (*j) = (*(n)); (*j) <= 2; (*j) += -1 {
			Dspmv(func() *[]byte {y := []byte("Lower"); return &y }(), (n), -(*one), (a), &((*(work))[0][(*j)-0]), func() *int {y := 1; return &y }(), zero, &((*(work))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label40:
		}
		Dspmv(func() *[]byte {y := []byte("Lower"); return &y }(), (n), -(*one), (a), &((*(ainv))[0]), func() *int {y := 1; return &y }(), zero, &((*(work))[0][0]), func() *int {y := 1; return &y }())
		//
	}
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
	//     End of Dppt03
	//
}
