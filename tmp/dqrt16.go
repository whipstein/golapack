package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dqrt16 computes the residual for a solution of a system of linear
// equations  A*x = b  or  A'*x = b:
//    resid = norm(B - A*X) / ( max(m,n) * norm(a) * norm(x) * eps),
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
//       SUBROUTinE Dqrt16( trans, m, n, nrhs, a, lda, x, ldx, b, ldb,
//                          rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       intEGER            lda, ldb, ldx, m, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), rwork(*),
//      $                   X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt16 computes the residual for a solution of a system of linear
// equations  A*x = b  or  A'*x = b:
//    resid = norm(B - A*X) / ( max(m,n) * norm(a) * norm(x) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the form of the system of equations:
//          = 'N':  A *x = b
//          = 'T':  A'*x = b, where A' is the transpose of A
//          = 'C':  A'*x = b, where A' is the transpose of A
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix A.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix A.  N >= 0.
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
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original M x N matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,M).
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
//          The leading dimension of the array X.  If trans = 'N',
//          ldx >= max(1,N); if trans = 'T' or 'C', ldx >= max(1,M).
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
//          The leading dimension of the array B.  IF trans = 'N',
//          ldb >= max(1,M); if trans = 'T' or 'C', ldb >= max(1,N).
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (m)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          The maximum over the number of right hand sides of
//          norm(B - A*X) / ( max(m,n) * norm(a) * norm(x) * eps).
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
func Dqrt16(trans *byte, m *int, n *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, b *[][]float64, ldb *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	j := new(int)
	n1 := new(int)
	n2 := new(int)
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
	//     Quick exit if M = 0 or N = 0 or nrhs = 0
	//
	if (*(m)) <= 0 || (*(n)) <= 0 || (*(nrhs)) == 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	if (*blas.Lsame((trans), func() *byte {y := byte('T'); return &y }())) || (*blas.Lsame((trans), func() *byte {y := byte('C'); return &y }())) {
		(*anorm) = (*Dlange(func() *byte {y := byte('I'); return &y }(), (m), (n), (a), (lda), (rwork)))
		(*n1) = (*(n))
		(*n2) = (*(m))
	} else {
		(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (a), (lda), (rwork)))
		(*n1) = (*(m))
		(*n2) = (*(n))
	}
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	//
	//     Compute  B - A*X  (or  B - A'*X) and store in B.
	//
	Dgemm((trans), func() *[]byte {y := []byte("No transpose"); return &y }(), n1, (nrhs), n2, -(*one), (a), (lda), (x), (ldx), one, (b), (ldb))
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(B - A*X) / ( max(m,n) * norm(a) * norm(x) * eps) .
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		(*bnorm) = (*Dasum(n1, &((*(b))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
		(*xnorm) = (*Dasum(n2, &((*(x))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
		if (*anorm) == (*zero) && (*bnorm) == (*zero) {
			(*(resid)) = (*zero)
		} else if (*anorm) <= (*zero) || (*xnorm) <= (*zero) {
			(*(resid)) = (*one) / (*eps)
		} else {
			(*(resid)) = (MAX((*(resid)), (((*bnorm)/(*anorm))/(*xnorm))/(MAX((*(m)), (*(n)))*(*eps))))
		}
		//Label10:
	}
	//
	return
	//
	//     End of Dqrt16
	//
}
