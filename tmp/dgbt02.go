package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dgbt02 computes the residual for a solution of a banded system of
// equations  A*x = b  or  A'*x = b:
//    resid = norm( B - A*X) / ( norm(a) * norm(x) * eps).
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
//       SUBROUTinE Dgbt02( trans, m, n, kl, ku, nrhs, a, lda, x, ldx, b,
//                          ldb, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       intEGER            kl, ku, lda, ldb, ldx, m, n, nrhs
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dgbt02 computes the residual for a solution of a banded system of
// equations  A*x = b  or  A'*x = b:
//    resid = norm( B - A*X) / ( norm(a) * norm(x) * eps).
// where eps is the machine precision.
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
// \param[in] kl
// \verbatim
//          kl is intEGER
//          The number of subdiagonals within the band of A.  kl >= 0.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//          The number of superdiagonals within the band of A.  ku >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of B.  nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original matrix A in band storage, stored in rows 1 to
//          kl+ku+1.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,kl+ku+1).
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
func Dgbt02(trans *byte, m *int, n *int, kl *int, ku *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, b *[][]float64, ldb *int, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i1 := new(int)
	i2 := new(int)
	j := new(int)
	kd := new(int)
	n1 := new(int)
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
	//     Quick return if N = 0 pr nrhs = 0
	//
	if (*(m)) <= 0 || (*(n)) <= 0 || (*(nrhs)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*kd) = (*(ku)) + 1
	(*anorm) = (*zero)
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		(*i1) = (MAX((*kd)+1-(*j), 1))
		(*i2) = (Min((*kd)+(*(m))-(*j), (*(kl))+(*kd)))
		(*anorm) = (MAX((*anorm), Dasum((*i2)-(*i1)+1, &((*(a))[(*i1)-1][(*j)-1]), func() *int {y := 1; return &y }())))
		//Label10:
	}
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	if (*blas.Lsame((trans), func() *byte {y := byte('T'); return &y }())) || (*blas.Lsame((trans), func() *byte {y := byte('C'); return &y }())) {
		(*n1) = (*(n))
	} else {
		(*n1) = (*(m))
	}
	//
	//     Compute  B - A*X (or  B - A'*X)
	//
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		Dgbmv((trans), (m), (n), (kl), (ku), -(*one), (a), (lda), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), one, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
		//Label20:
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(B - A*X) / ( norm(a) * norm(x) * eps).
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		(*bnorm) = (*Dasum(n1, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }()))
		(*xnorm) = (*Dasum(n1, &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }()))
		if (*xnorm) <= (*zero) {
			(*(resid)) = (*one) / (*eps)
		} else {
			(*(resid)) = (MAX((*(resid)), (((*bnorm)/(*anorm))/(*xnorm))/(*eps)))
		}
		//Label30:
	}
	//
	return
	//
	//     End of Dgbt02
	//
}
