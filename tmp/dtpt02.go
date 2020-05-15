package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtpt02 computes the residual for the computed solution to a
// triangular system of linear equations  A*x = b  or  A'*x = b  when
// the triangular matrix A is stored in packed format.  Here A' is the
// transpose of A and x and b are N by nrhs matrices.  The test ratio is
// the maximum over the number of right hand sides of
//    norm(b - op(a)*x) / ( norm(op(a)) * norm(x) * eps),
// where op(a) denotes A or A' and eps is the machine epsilon.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dtpt02( uplo, trans, diag, n, nrhs, AP, x, ldx, b, ldb,
//                          work, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       inTEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   AP(*), B( ldb, *), work(*), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtpt02 computes the residual for the computed solution to a
// triangular system of linear equations  A*x = b  or  A'*x = b  when
// the triangular matrix A is stored in packed format.  Here A' is the
// transpose of A and x and b are N by nrhs matrices.  The test ratio is
// the maximum over the number of right hand sides of
//    norm(b - op(a)*x) / ( norm(op(a)) * norm(x) * eps),
// where op(a) denotes A or A' and eps is the machine epsilon.
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
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the operation applied to A.
//          = 'N':  A *x = b  (No transpose)
//          = 'T':  A'*x = b  (Transpose)
//          = 'C':  A'*x = b  (Conjugate transpose = Transpose)
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
//          N is inTEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The number of right hand sides, i.e., the number of columns
//          of the matrices X and B.  nrhs >= 0.
// \endverbatim
//
// \param[in] AP
// \verbatim
//          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The upper or lower triangular matrix a, packed columnwise in
//          a linear array.  The j-th column of A is stored in the array
//          AP as follows:
//          if uplo = 'U', AP((j-1)*j/2 + i) = a(i,j) for 1<=i<=j;
//          if uplo = 'L',
//             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = a(i,j) for j<=i<=n.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The computed solution vectors for the system of linear
//          equations.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is inTEGER
//          The leading dimension of the array X.  ldx >= max(1,N).
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          The right hand side vectors for the system of linear
//          equations.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B.  ldb >= max(1,N).
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
//          The maximum over the number of right hand sides of
//          norm(op(a)*x - b) / ( norm(op(a)) * norm(x) * eps).
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
func Dtpt02(uplo *byte, trans *byte, diag *byte, n *int, nrhs *int, ap *[]float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, work *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	j := new(int)
	anorm := new(float64)
	bnorm := new(float64)
	eps := new(float64)
	xnorm := new(float64)
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
	//     Quick exit if N = 0 or nrhs = 0
	//
	if (*(n)) <= 0 || (*(nrhs)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Compute the 1-norm of A or A'.
	//
	if blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		(*anorm) = (*Dlantp(func() *byte {y := byte('1'); return &y }(), (uplo), (diag), (n), (ap), (work)))
	} else {
		(*anorm) = (*Dlantp(func() *byte {y := byte('I'); return &y }(), (uplo), (diag), (n), (ap), (work)))
	}
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(op(a)*x - b) / ( norm(op(a)) * norm(x) * eps).
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		Dcopy((n), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		Dtpmv((uplo), (trans), (diag), (n), (ap), (work), func() *int {y := 1; return &y }())
		Daxpy((n), -(*one), &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		(*bnorm) = (*Dasum((n), (work), func() *int {y := 1; return &y }()))
		(*xnorm) = (*Dasum((n), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }()))
		if (*xnorm) <= (*zero) {
			(*(resid)) = (*one) / (*eps)
		} else {
			(*(resid)) = (MAX((*(resid)), (((*bnorm)/(*anorm))/(*xnorm))/(*eps)))
		}
		//Label10:
	}
	//
	return
	//
	//     End of Dtpt02
	//
}
