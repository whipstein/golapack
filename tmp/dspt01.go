package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dspt01 re_constructs a symmetric indefinite packed matrix A from its
// block L*D*l' or U*D*U' factorization and computes the residual
//      norm( C - A) / ( N * norm(a) * eps),
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
//       SUBROUTinE Dspt01( uplo, n, a, afac, ipiv, c, ldc, rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            ldc, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       intEGER            ipiv(*)
//       DOUBLE PRECISION   a(*), afac(*), c( ldc, *), rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dspt01 re_constructs a symmetric indefinite packed matrix A from its
// block L*D*l' or U*D*U' factorization and computes the residual
//      norm( C - A) / ( N * norm(a) * eps),
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
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The original symmetric matrix a, stored as a packed
//          triangular matrix.
// \endverbatim
//
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The factored form of the matrix a, stored as a packed
//          triangular matrix.  afac contains the block diagonal matrix D
//          and the multipliers used to obtain the factor L or U from the
//          block L*D*l' or U*D*U' factorization as computed by Dsptrf.
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is intEGER array, dimension (n)
//          The pivot indices from Dsptrf.
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
// \date December 2016
//
// \ingroup double_lin
//
//  =====================================================================
func Dspt01(uplo *byte, n *int, a *[]float64, afac *[]float64, ipiv *[]int, c *[][]float64, ldc *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	jc := new(int)
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
	(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (rwork)))
	//
	//     Initialize C to the identity matrix.
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, one, (c), (ldc))
	//
	//     Call dlavsp to form the product d * U' (or d * L').
	//
	dlavsp((uplo), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (n), (n), (afac), (ipiv), (c), (ldc), info)
	//
	//     Call dlavsp again to multiply by U ( or L).
	//
	dlavsp((uplo), func() *[]byte {y := []byte("No transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (n), (n), (afac), (ipiv), (c), (ldc), info)
	//
	//     Compute the difference  C - A .
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		(*jc) = 0
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j); (*i)++ {
				(*(c))[(*i)-1][(*j)-1] = (*(c))[(*i)-1][(*j)-1] - (*(a))[(*jc)+(*i)-1]
				//Label10:
			}
			(*jc) = (*jc) + (*j)
			//Label20:
		}
	} else {
		(*jc) = 1
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
				(*(c))[(*i)-1][(*j)-1] = (*(c))[(*i)-1][(*j)-1] - (*(a))[(*jc)+(*i)-(*j)-1]
				//Label30:
			}
			(*jc) = (*jc) + (*(n)) - (*j) + 1
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
	//     End of Dspt01
	//
}
