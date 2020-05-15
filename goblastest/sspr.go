package main

// Sspr ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
//       INTEGER INCX,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       REAL AP(*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSPR    performs the symmetric rank 1 operation
//
//    A := alpha*x*x**T + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
// n by n symmetric matrix, supplied in packed form.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the matrix A is supplied in the packed
//           array AP as follows:
//
//              UPLO = 'U' or 'u'   The upper triangular part of A is
//                                  supplied in AP.
//
//              UPLO = 'L' or 'l'   The lower triangular part of A is
//                                  supplied in AP.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the n
//           element vector x.
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//           On entry, INCX specifies the increment for the elements of
//           X. INCX must not be zero.
// \endverbatim
//
// \param[in,out] AP
// \verbatim
//          AP is REAL array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on. On exit, the array
//           AP is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on. On exit, the array
//           AP is overwritten by the lower triangular part of the
//           updated matrix.
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
// \ingroup single_blas_level2
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  Level 2 Blas routine.
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Sspr(uplo byte, n int, alpha float32, x *[]float32, incx int, ap *[]float32) {
	var temp float32
	var info, ix, jx, k, kk, kx int
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(uplo, 'U') && !Lsame(uplo, 'L') {
		info = 1
	} else if n < 0 {
		info = 2
	} else if incx == 0 {
		info = 5
	}
	if info != 0 {
		Xerbla("Sspr", info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) || (alpha == 0) {
		return
	}
	//
	//     Set the start point in X if the increment is not unity.
	//
	if incx <= 0 {
		kx = 1 - (n-1)*incx
	} else if incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of the array AP
	//     are accessed sequentially with one pass through AP.
	//
	kk = 1
	if Lsame(uplo, 'U') {
		//
		//        Form  A  when upper triangle is stored in AP.
		//
		if incx == 1 {
			for j := 1; j <= n; j++ {
				if (*x)[j-(1)] != 0 {
					temp = alpha * (*x)[j-(1)]
					k = kk
					for i := 1; i <= j; i++ {
						(*ap)[k-(1)] += (*x)[i-(1)] * temp
						k++
					}
				}
				kk += j
			}
		} else {
			jx = kx
			for j := 1; j <= n; j++ {
				if (*x)[jx-(1)] != 0 {
					temp = alpha * (*x)[jx-(1)]
					ix = kx
					for k = kk; k <= kk+j-1; k++ {
						(*ap)[k-(1)] += (*x)[ix-(1)] * temp
						ix += incx
					}
				}
				jx += incx
				kk += j
			}
		}
	} else {
		//
		//        Form  A  when lower triangle is stored in AP.
		//
		if incx == 1 {
			for j := 1; j <= n; j++ {
				if (*x)[j-(1)] != 0 {
					temp = alpha * (*x)[j-(1)]
					k = kk
					for i := j; i <= n; i++ {
						(*ap)[k-(1)] += (*x)[i-(1)] * temp
						k++
					}
				}
				kk += n - j + 1
			}
		} else {
			jx = kx
			for j := 1; j <= n; j++ {
				if (*x)[jx-(1)] != 0 {
					temp = alpha * (*x)[jx-(1)]
					ix = jx
					for k = kk; k <= kk+n-j; k++ {
						(*ap)[k-(1)] += (*x)[ix-(1)] * temp
						ix += incx
					}
				}
				jx += incx
				kk += n - j + 1
			}
		}
	}

	return
	//
	//     End of Sspr  .
	//
}
