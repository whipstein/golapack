package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtbt03 computes the residual for the solution to a scaled triangular
// system of equations  A*x = s*b  or  A'*x = s*b  when A is a
// triangular band matrix. Here A' is the transpose of a, s is a scalar,
// and x and b are N by nrhs matrices.  The test ratio is the maximum
// over the number of right hand sides of
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
//       SUBROUTinE Dtbt03( uplo, trans, diag, n, kd, nrhs, AB, ldab,
//                          scale, cnorm, tscal, x, ldx, b, ldb, work,
//                          resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       intEGER            kd, ldab, ldb, ldx, n, nrhs
//       DOUBLE PRECISION   resid, scale, tscal
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   AB( ldab, *), B( ldb, *), cnorm(*),
//      $                   work(*), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtbt03 computes the residual for the solution to a scaled triangular
// system of equations  A*x = s*b  or  A'*x = s*b  when A is a
// triangular band matrix. Here A' is the transpose of a, s is a scalar,
// and x and b are N by nrhs matrices.  The test ratio is the maximum
// over the number of right hand sides of
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
//          N is intEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kd
// \verbatim
//          kd is intEGER
//          The number of superdiagonals or subdiagonals of the
//          triangular band matrix A.  kd >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of columns
//          of the matrices X and B.  nrhs >= 0.
// \endverbatim
//
// \param[in] AB
// \verbatim
//          AB is DOUBLE PRECISION array, dimension (ldab,N)
//          The upper or lower triangular band matrix a, stored in the
//          first kd+1 rows of the array. The j-th column of A is stored
//          in the j-th column of the array AB as follows:
//          if uplo = 'U', AB(kd+1+i-j,j) = a(i,j) for max(1,j-kd)<=i<=j;
//          if uplo = 'L', AB(1+i-j,j)    = a(i,j) for j<=i<=min(n,j+kd).
// \endverbatim
//
// \param[in] ldab
// \verbatim
//          ldab is intEGER
//          The leading dimension of the array AB.  ldab >= kd+1.
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
//          cnorm actually contains the column norms of tscal*a.
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
//          ldb is intEGER
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
func Dtbt03(uplo *byte, trans *byte, diag *byte, n *int, kd *int, nrhs *int, ab *[][]float64, ldab *int, scale *float64, cnorm *[]float64, tscal *float64, x *[][]float64, ldx *int, b *[][]float64, ldb *int, work *[]float64, resid *float64) {
	one := new(float64)
	zero := new(float64)
	ix := new(int)
	j := new(int)
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
	//     Quick exit if N = 0
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
	//     norms already computed by Dlatbs.
	//
	(*tnorm) = (*zero)
	if blas.Lsame((diag), func() *byte {y := byte('N'); return &y }()) {
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*tnorm) = (MAX((*tnorm), (*(tscal))*ABS(((*(ab))[(*(kd))+0][(*j)-(1)]))+(*(cnorm))[(*j)-(1)]))
				//Label10:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*tnorm) = (MAX((*tnorm), (*(tscal))*ABS(((*(ab))[0][(*j)-(1)]))+(*(cnorm))[(*j)-(1)]))
				//Label20:
			}
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			(*tnorm) = (MAX((*tnorm), (*(tscal))+(*(cnorm))[(*j)-(1)]))
			//Label30:
		}
	}
	//
	//     Compute the maximum over the number of right hand sides of
	//        norm(op(a)*x - s*b) / ( norm(op(a)) * norm(x) * eps).
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		Dcopy((n), &((*(x))[0][(*j)-(1)]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		(*ix) = (*Idamax((n), (work), func() *int {y := 1; return &y }()))
		(*xnorm) = (MAX((*one), ABS(((*(x))[(*ix)-(1)][(*j)-(1)]))))
		(*xscal) = ((*one) / (*xnorm)) / DBLE((*(kd))+1)
		Dscal((n), xscal, (work), func() *int {y := 1; return &y }())
		dtbmv((uplo), (trans), (diag), (n), (kd), (ab), (ldab), (work), func() *int {y := 1; return &y }())
		Daxpy((n), -(*(scale))*(*xscal), &((*(b))[0][(*j)-(1)]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		(*ix) = (*Idamax((n), (work), func() *int {y := 1; return &y }()))
		(*err) = (*(tscal)) * ABS(((*(work))[(*ix)-(1)]))
		(*ix) = (*Idamax((n), &((*(x))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
		(*xnorm) = (ABS(((*(x))[(*ix)-(1)][(*j)-(1)])))
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
	//     End of Dtbt03
	//
}
