package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dsyt01 re_constructs a symmetric indefinite matrix A from its
// block L*D*L' or U*D*U' factorization and computes the residual
//    norm( C - A) / ( N * norm(a) * eps),
// where C is the re_constructed matrix and eps is the machine epsilon.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dsyt01_Aa( uplo, n, a, lda, afac, ldafac, ipiv, c, ldc,
//                             rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       inTEGER            lda, ldafac, ldc, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       inTEGER            ipiv(*)
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *), c( ldc, *),
//      $                   rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dsyt01 re_constructs a symmetric indefinite matrix A from its
// block L*D*L' or U*D*U' factorization and computes the residual
//    norm( C - A) / ( N * norm(a) * eps),
// where C is the re_constructed matrix and eps is the machine epsilon.
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
//          N is inTEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original symmetric matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.  lda >= max(1,N)
// \endverbatim
//
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factored form of the matrix A.  afac contains the block
//          diagonal matrix D and the multipliers used to obtain the
//          factor L or U from the block L*D*L' or U*D*U' factorization
//          as computed by Dsytrf.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is inTEGER
//          The leading dimension of the array afac.  ldafac >= max(1,N).
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is inTEGER array, dimension (n)
//          The pivot indices from Dsytrf.
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION array, dimension (ldc,N)
// \endverbatim
//
// \param[in] ldc
// \verbatim
//          ldc is inTEGER
//          The leading dimension of the array C.  ldc >= max(1,N).
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
//          If uplo = 'L', norm(L*D*L' - A) / ( N * norm(a) * eps)
//          If uplo = 'U', norm(U*D*U' - A) / ( N * norm(a) * eps)
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
//  @precisions fortran d -> z c
//
// \ingroup double_lin
//
//  =====================================================================
func Dsyt01_Aa(uplo *byte, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, ipiv *[]int, c *[][]float64, ldc *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
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
	//     Quick exit if N = 0.
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Determine eps and the norm of A.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (lda), (rwork)))
	//
	//     Initialize C to the tridiagonal matrix T.
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, zero, (c), (ldc))
	Dlacpy(func() *byte {y := byte('F'); return &y }(), func() *int {y := 1; return &y }(), (n), &((*(afac))[0][0]), (*(ldafac))+1, &((*(c))[0][0]), (*(ldc))+1)
	if (*(n)) > 1 {
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			Dlacpy(func() *byte {y := byte('F'); return &y }(), func() *int {y := 1; return &y }(), (*(n))-1, &((*(afac))[0][1]), (*(ldafac))+1, &((*(c))[0][1]), (*(ldc))+1)
			Dlacpy(func() *byte {y := byte('F'); return &y }(), func() *int {y := 1; return &y }(), (*(n))-1, &((*(afac))[0][1]), (*(ldafac))+1, &((*(c))[1][0]), (*(ldc))+1)
		} else {
			Dlacpy(func() *byte {y := byte('F'); return &y }(), func() *int {y := 1; return &y }(), (*(n))-1, &((*(afac))[1][0]), (*(ldafac))+1, &((*(c))[0][1]), (*(ldc))+1)
			Dlacpy(func() *byte {y := byte('F'); return &y }(), func() *int {y := 1; return &y }(), (*(n))-1, &((*(afac))[1][0]), (*(ldafac))+1, &((*(c))[1][0]), (*(ldc))+1)
		}
		//
		//        Call Dtrmm to form the product U' * D (or L * D).
		//
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			Dtrmm(func() *[]byte {y := []byte("Left"); return &y }(), (uplo), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (*(n))-1, (n), one, &((*(afac))[0][1]), (ldafac), &((*(c))[1][0]), (ldc))
		} else {
			Dtrmm(func() *[]byte {y := []byte("Left"); return &y }(), (uplo), func() *[]byte {y := []byte("No transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (*(n))-1, (n), one, &((*(afac))[1][0]), (ldafac), &((*(c))[1][0]), (ldc))
		}
		//
		//        Call Dtrmm again to multiply by U (or L).
		//
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			Dtrmm(func() *[]byte {y := []byte("Right"); return &y }(), (uplo), func() *[]byte {y := []byte("No transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (n), (*(n))-1, one, &((*(afac))[0][1]), (ldafac), &((*(c))[0][1]), (ldc))
		} else {
			Dtrmm(func() *[]byte {y := []byte("Right"); return &y }(), (uplo), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (n), (*(n))-1, one, &((*(afac))[1][0]), (ldafac), &((*(c))[0][1]), (ldc))
		}
	}
	//
	//     Apply symmetric pivots
	//
	for (*j) = (*(n)); (*j) <= 1; (*j) += -1 {
		(*i) = (*(ipiv))[(*j)-1]
		if (*i) != (*j) {
			Dswap((n), &((*(c))[(*j)-1][0]), (ldc), &((*(c))[(*i)-1][0]), (ldc))
		}
	}
	for (*j) = (*(n)); (*j) <= 1; (*j) += -1 {
		(*i) = (*(ipiv))[(*j)-1]
		if (*i) != (*j) {
			Dswap((n), &((*(c))[0][(*j)-1]), func() *int {y := 1; return &y }(), &((*(c))[0][(*i)-1]), func() *int {y := 1; return &y }())
		}
	}
	//
	//
	//     Compute the difference  C - A .
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j); (*i)++ {
				(*(c))[(*i)-1][(*j)-1] = (*(c))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
			}
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
				(*(c))[(*i)-1][(*j)-1] = (*(c))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
			}
		}
	}
	//
	//     Compute norm( C - A) / ( N * norm(a) * eps)
	//
	(*(resid)) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (c), (ldc), (rwork)))
	//
	if (*anorm) <= (*zero) {
		if (*(resid)) != (*zero) {
			(*(resid)) = (*one) / (*eps)
		}
	} else {
		(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	}
	//
	return
	//
	//     End of Dsyt01
	//
}
