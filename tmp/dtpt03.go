package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtpt03 computes the residual for the solution to a scaled triangular
// system of equations A*x = s*b  or  A'*x = s*b  when the triangular
// matrix A is stored in packed format.  Here A' is the transpose of a,
// s is a scalar, and x and b are N by nrhs matrices.  The test ratio is
// the maximum over the number of right hand sides of
//    norm(s*b - op(a)*x) / ( norm(op(a)) * norm(x) * eps),
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
//       SUBROUTinE Dtpt03( uplo, trans, diag, n, nrhs, AP, scale, cnorm,
//                          tscal, x, ldx, b, ldb, work, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       inTEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid, scale, tscal
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   AP(*), B( ldb, *), cnorm(*), work(*),
//      $                   X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtpt03 computes the residual for the solution to a scaled triangular
// system of equations A*x = s*b  or  A'*x = s*b  when the triangular
// matrix A is stored in packed format.  Here A' is the transpose of a,
// s is a scalar, and x and b are N by nrhs matrices.  The test ratio is
// the maximum over the number of right hand sides of
//    norm(s*b - op(a)*x) / ( norm(op(a)) * norm(x) * eps),
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
//          = 'N':  A *x = s*b  (No transpose)
//          = 'T':  A'*x = s*b  (Transpose)
//          = 'C':  A'*x = s*b  (Conjugate transpose = Transpose)
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
// \param[in] scale
// \verbatim
//          scale is DOUBLE PRECISION
//          The scaling factor s used in solving the triangular system.
// \endverbatim
//
// \param[in] cnorm
// \verbatim
//          cnorm is DOUBLE PRECISION array, dimension (n)
//          The 1-norms of the columns of a, not counting the diagonal.
// \endverbatim
//
// \param[in] tscal
// \verbatim
//          tscal is DOUBLE PRECISION
//          The scaling factor used in computing the 1-norms in cnorm.
//          cnorm actually contains the column norms of tscal*A.
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
//          norm(op(a)*x - s*b) / ( norm(op(a)) * norm(x) * eps).
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
func Dtpt03(uplo *byte, trans *byte, diag *byte, n *int, nrhs *int, ap *[]float64, scale *float64, cnorm *[]float64, tscal *float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, work *[]float64, resid *float64) {
	one := new(float64)
	zero := new(float64)
	ix := new(int)
	j := new(int)
	jj := new(int)
	bignum := new(float64)
	eps := new(float64)
	err := new(float64)
	smlnum := new(float64)
	tnorm := new(float64)
	xnorm := new(float64)
	xscal := new(float64)
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
	//     Quick exit if N = 0.
	//
	if (*(n)) <= 0 || (*(nrhs)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*smlnum) = (*Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()))
	(*bignum) = (*one) / (*smlnum)
	Dlabad(smlnum, bignum)
	//
	//     Compute the norm of the triangular matrix A using the column
	//     norms already computed by Dlatps.
	//
	(*tnorm) = (*zero)
	if blas.Lsame((diag), func() *byte {y := byte('N'); return &y }()) {
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			(*jj) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*tnorm) = (MAX((*tnorm), (*(tscal))*ABS(((*(ap))[(*jj)-1]))+(*(cnorm))[(*j)-1]))
				(*jj) = (*jj) + (*j) + 1
				//Label10:
			}
		} else {
			(*jj) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*tnorm) = (MAX((*tnorm), (*(tscal))*ABS(((*(ap))[(*jj)-1]))+(*(cnorm))[(*j)-1]))
				(*jj) = (*jj) + (*(n)) - (*j) + 1
				//Label20:
			}
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			(*tnorm) = (MAX((*tnorm), (*(tscal))+(*(cnorm))[(*j)-1]))
			//Label30:
		}
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(op(a)*x - s*b) / ( norm(op(a)) * norm(x) * eps).
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		Dcopy((n), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		(*ix) = (*Idamax((n), (work), func() *int {y := 1; return &y }()))
		(*xnorm) = (MAX((*one), ABS(((*(x))[(*ix)-1][(*j)-1]))))
		(*xscal) = ((*one) / (*xnorm)) / DBLE((*(n)))
		Dscal((n), xscal, (work), func() *int {y := 1; return &y }())
		Dtpmv((uplo), (trans), (diag), (n), (ap), (work), func() *int {y := 1; return &y }())
		Daxpy((n), -(*(scale))*(*xscal), &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		(*ix) = (*Idamax((n), (work), func() *int {y := 1; return &y }()))
		(*err) = (*(tscal)) * ABS(((*(work))[(*ix)-1]))
		(*ix) = (*Idamax((n), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }()))
		(*xnorm) = (ABS(((*(x))[(*ix)-1][(*j)-1])))
		if (*err)*(*smlnum) <= (*xnorm) {
			if (*xnorm) > (*zero) {
				(*err) = (*err) / (*xnorm)
			}
		} else {
			if (*err) > (*zero) {
				(*err) = (*one) / (*eps)
			}
		}
		if (*err)*(*smlnum) <= (*tnorm) {
			if (*tnorm) > (*zero) {
				(*err) = (*err) / (*tnorm)
			}
		} else {
			if (*err) > (*zero) {
				(*err) = (*one) / (*eps)
			}
		}
		(*(resid)) = (MAX((*(resid)), (*err)))
		//Label40:
	}
	//
	return
	//
	//     End of Dtpt03
	//
}
