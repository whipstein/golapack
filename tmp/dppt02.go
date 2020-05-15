package goblas

import 

// Dppt02 computes the residual in the solution of a symmetric system
// of linear equations  A*x = b  when packed storage is used for the
// coefficient matrix.  The ratio computed is
//
//    resid = norm(B - A*X) / ( norm(a) * norm(x) * eps),
//
// where eps is the machine precision.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dppt02( uplo, n, nrhs, a, x, ldx, b, ldb, rwork,
//                          resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a(*), B( ldb, *), rwork(*), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dppt02 computes the residual in the solution of a symmetric system
// of linear equations  A*x = b  when packed storage is used for the
// coefficient matrix.  The ratio computed is
//
//    resid = norm(B - A*X) / ( norm(a) * norm(x) * eps),
//
// where eps is the machine precision.
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
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of b, the matrix of right hand sides.
//          nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The original symmetric matrix a, stored as a packed
//          triangular matrix.
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
//          ldx is intEGER
//          The leading dimension of the array X.   ldx >= max(1,N).
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          On entry, the right hand side vectors for the system of
//          linear equations.
//          On exit, B is overwritten with the difference B - A*X.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B.  ldb >= max(1,N).
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          The maximum over the number of right hand sides of
//          norm(B - A*X) / ( norm(a) * norm(x) * eps).
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
func Dppt02(uplo *byte, n *int, nrhs *int, a *[]float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, rwork *[]float64, resid *float64) {
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
	//     Quick exit if N = 0 or nrhs = 0.
	//
	if (*(n)) <= 0 || (*(nrhs)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (rwork)))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute  B - A*X  for the matrix of right hand sides B.
	//
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		Dspmv((uplo), (n), -(*one), (a), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), one, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
		//Label10:
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm( B - A*X) / ( norm(a) * norm(x) * eps) .
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
		//Label20:
	}
	//
	return
	//
	//     End of Dppt02
	//
}
