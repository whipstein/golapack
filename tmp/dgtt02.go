package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dgtt02 computes the residual for the solution to a tridiagonal
// system of equations:
//    resid = norm(B - op(a)*X) / (norm(a) * norm(x) * eps),
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
//       SUBROUTinE Dgtt02( trans, n, nrhs, DL, d, du, x, ldx, b, ldb,
//                          resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       intEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   B( ldb, *), d(*), DL(*), du(*),
//      $                   X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dgtt02 computes the residual for the solution to a tridiagonal
// system of equations:
//    resid = norm(B - op(a)*X) / (norm(a) * norm(x) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER
//          Specifies the form of the residual.
//          = 'N':  B - A * X  (No transpose)
//          = 'T':  B - A'* X  (Transpose)
//          = 'C':  B - A'* X  (Conjugate transpose = Transpose)
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGTER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of columns
//          of the matrices B and X.  nrhs >= 0.
// \endverbatim
//
// \param[in] DL
// \verbatim
//          DL is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) sub-diagonal elements of A.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (n)
//          The diagonal elements of A.
// \endverbatim
//
// \param[in] du
// \verbatim
//          du is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) super-diagonal elements of A.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The computed solution vectors X.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is intEGER
//          The leading dimension of the array X.  ldx >= max(1,N).
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          On entry, the right hand side vectors for the system of
//          linear equations.
//          On exit, B is overwritten with the difference B - op(a)*X.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B.  ldb >= max(1,N).
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(B - op(a)*X) / (norm(a) * norm(x) * eps)
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
func Dgtt02(trans *byte, n *int, nrhs *int, dl *[]float64, d *[]float64, du *[]float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, resid *float64) {
	one := new(float64)
	zero := new(float64)
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
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
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
	(*(resid)) = (*zero)
	if (*(n)) <= 0 || (*(nrhs)) == 0 {
		return
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(B - op(a)*X) / ( norm(a) * norm(x) * eps).
	//
	if blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		(*anorm) = (*Dlangt(func() *byte {y := byte('1'); return &y }(), (n), (dl), (d), (du)))
	} else {
		(*anorm) = (*Dlangt(func() *byte {y := byte('I'); return &y }(), (n), (dl), (d), (du)))
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
	//     Compute B - op(a)*X.
	//
	Dlagtm((trans), (n), (nrhs), -(*one), (dl), (d), (du), (x), (ldx), one, (b), (ldb))
	//
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		(*bnorm) = (*Dasum((n), &((*(b))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
		(*xnorm) = (*Dasum((n), &((*(x))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
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
	//     End of Dgtt02
	//
}
