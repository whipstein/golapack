package goblas

import 

// Dptt02 computes the residual for the solution to a symmetric
// tridiagonal system of equations:
//    resid = norm(B - A*X) / (norm(a) * norm(x) * eps),
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
//       SUBROUTinE Dptt02( n, nrhs, d, e, x, ldx, b, ldb, resid)
//
//       .. Scalar Arguments ..
//       intEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   B( ldb, *), d(*), E(*), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dptt02 computes the residual for the solution to a symmetric
// tridiagonal system of equations:
//    resid = norm(B - A*X) / (norm(a) * norm(x) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGTER
//          The order of the matrix A.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of columns
//          of the matrices B and X.  nrhs >= 0.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (n)
//          The n diagonal elements of the tridiagonal matrix A.
// \endverbatim
//
// \param[in] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) subdiagonal elements of the tridiagonal matrix A.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The n by nrhs matrix of solution vectors X.
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
//          On entry, the n by nrhs matrix of right hand side vectors B.
//          On exit, B is overwritten with the difference B - A*X.
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
//          norm(B - A*X) / (norm(a) * norm(x) * eps)
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
func Dptt02(n *int, nrhs *int, d *[]float64, e *[]float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, resid *float64) {
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
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Quick return if possible
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Compute the 1-norm of the tridiagonal matrix A.
	//
	(*anorm) = (*dlanst(func() *byte {y := byte('1'); return &y }(), (n), (d), (e)))
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute B - A*X.
	//
	Dlaptm((n), (nrhs), -(*one), (d), (e), (x), (ldx), one, (b), (ldb))
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(B - A*X) / ( norm(a) * norm(x) * eps).
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		(*bnorm) = (*Dasum((n), &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }()))
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
	//     End of Dptt02
	//
}
