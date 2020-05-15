package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dsyt01 re_constructs a symmetric indefinite matrix A from its
// block L*D*l' or U*D*U' factorization and computes the residual
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
//       SUBROUTinE Dsyt01( uplo, n, a, lda, afac, ldafac, ipiv, c, ldc,
//                          rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            lda, ldafac, ldc, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       intEGER            ipiv(*)
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
// block L*D*l' or U*D*U' factorization and computes the residual
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
//          N is intEGER
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
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,N)
// \endverbatim
//
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factored form of the matrix A.  afac contains the block
//          diagonal matrix D and the multipliers used to obtain the
//          factor L or U from the block L*D*l' or U*D*U' factorization
//          as computed by Dsytrf.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is intEGER
//          The leading dimension of the array afac.  ldafac >= max(1,N).
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is intEGER array, dimension (n)
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
//          ldc is intEGER
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
//          If uplo = 'L', norm(L*D*l' - A) / ( N * norm(a) * eps)
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
// \date November 2013
//
// \ingroup double_lin
//
//  =====================================================================
func Dsyt01(uplo *byte, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, ipiv *[]int, c *[][]float64, ldc *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	anorm := new(float64)
	eps := new(float64)
	//
	//  -- lapACK test routine (version 3.5.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2013
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
	//     Initialize C to the identity matrix.
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, one, (c), (ldc))
	//
	//     Call dlavsy to form the product d * U' (or d * L').
	//
	dlavsy((uplo), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (n), (n), (afac), (ldafac), (ipiv), (c), (ldc), info)
	//
	//     Call dlavsy again to multiply by U (or L).
	//
	dlavsy((uplo), func() *[]byte {y := []byte("No transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (n), (n), (afac), (ldafac), (ipiv), (c), (ldc), info)
	//
	//     Compute the difference  C - A .
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j); (*i)++ {
				(*(c))[(*i)-(1)][(*j)-(1)] = (*(c))[(*i)-(1)][(*j)-(1)] - (*(a))[(*i)-(1)][(*j)-(1)]
				//Label10:
			}
			//Label20:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
				(*(c))[(*i)-(1)][(*j)-(1)] = (*(c))[(*i)-(1)][(*j)-(1)] - (*(a))[(*i)-(1)][(*j)-(1)]
				//Label30:
			}
			//Label40:
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
