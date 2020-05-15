package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dsyt01_3 re_constructs a symmetric indefinite matrix A from its
// block L*D*L' or U*D*U' factorization computed by DsytrfRk
// (or Dsytrf_Bk) and computes the residual
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
//       SUBROUTinE Dsyt01_3( uplo, n, a, lda, afac, ldafac, e, ipiv, c,
//                            ldc, rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       inTEGER            lda, ldafac, ldc, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       inTEGER            ipiv(*)
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *), c( ldc, *),
//      $                   E(*), rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dsyt01_3 re_constructs a symmetric indefinite matrix A from its
// block L*D*L' or U*D*U' factorization computed by DsytrfRk
// (or Dsytrf_Bk) and computes the residual
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
//          diagonal of the block diagonal matrix D and factors U or L
//          as computed by DsytrfRk and Dsytrf_BK:
//            a) OnlY diagonal elements of the symmetric block diagonal
//               matrix D on the diagonal of a, i.e. d(k,k) = a(k,k);
//               (superdiagonal (or subdiagonal) elements of D
//                should be provided on entry in array E), and
//            b) If uplo = 'U': factor U in the superdiagonal part of A.
//               If uplo = 'L': factor L in the subdiagonal part of A.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is inTEGER
//          The leading dimension of the array afac.
//          ldafac >= max(1,N).
// \endverbatim
//
// \param[in] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (n)
//          On entry, contains the superdiagonal (or subdiagonal)
//          elements of the symmetric block diagonal matrix D
//          with 1-by-1 or 2-by-2 diagonal blocks, where
//          If uplo = 'U': E(i) = d(i-1,i),i=2:N, E1 not referenced;
//          If uplo = 'L': E(i) = d(i+1,i),i=1:N-1, E(n) not referenced.
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is inTEGER array, dimension (n)
//          The pivot indices from DsytrfRk (or Dsytrf_Bk).
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
// \date June 2017
//
// \ingroup double_lin
//
//  =====================================================================
func Dsyt01_3(uplo *byte, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, e *[]float64, ipiv *[]int, c *[][]float64, ldc *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	anorm := new(float64)
	eps := new(float64)
	//
	//  -- lapACK test routine (version 3.7.1) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2017
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
	//     a) Revert to multiplyers of L
	//
	DsyconvfRook((uplo), func() *byte {y := byte('R'); return &y }(), (n), (afac), (ldafac), (e), (ipiv), info)
	//
	//     1) Determine eps and the norm of A.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (lda), (rwork)))
	//
	//     2) Initialize C to the identity matrix.
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, one, (c), (ldc))
	//
	//     3) Call dlavsyRook to form the product d * U' (or d * L').
	//
	dlavsyRook((uplo), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (n), (n), (afac), (ldafac), (ipiv), (c), (ldc), info)
	//
	//     4) Call dlavsyRook again to multiply by U (or L).
	//
	dlavsyRook((uplo), func() *[]byte {y := []byte("No transpose"); return &y }(), func() *[]byte {y := []byte("Unit"); return &y }(), (n), (n), (afac), (ldafac), (ipiv), (c), (ldc), info)
	//
	//     5) Compute the difference  C - A.
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
	//     6) Compute norm( C - A) / ( N * norm(a) * eps)
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
	//     b) Convert to factor of L (or U)
	//
	DsyconvfRook((uplo), func() *byte {y := byte('C'); return &y }(), (n), (afac), (ldafac), (e), (ipiv), info)
	//
	return
	//
	//     End of Dsyt01_3
	//
}
