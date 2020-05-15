package main

// Ssyr ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
//       INTEGER INCX,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSYR   performs the symmetric rank 1 operation
//
//    A := alpha*x*x**T + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
// n by n symmetric matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the array A is to be referenced as
//           follows:
//
//              UPLO = 'U' or 'u'   Only the upper triangular part of A
//                                  is to be referenced.
//
//              UPLO = 'L' or 'l'   Only the lower triangular part of A
//                                  is to be referenced.
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
// \param[in,out] A
// \verbatim
//          A is REAL array, dimension ( LDA, N )
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular part of the symmetric matrix and the strictly
//           lower triangular part of A is not referenced. On exit, the
//           upper triangular part of the array A is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular part of the symmetric matrix and the strictly
//           upper triangular part of A is not referenced. On exit, the
//           lower triangular part of the array A is overwritten by the
//           lower triangular part of the updated matrix.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           max( 1, n ).
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
func Ssyr(uplo byte, n int, alpha float32, x *[]float32, incx int, a *[][]float32, lda int) {
	var temp float32
	var info, ix, jx, kx int
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
	} else if lda < max(1, n) {
		info = 7
	}
	if info != 0 {
		Xerbla("Ssyr", info)
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
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through the triangular part
	//     of A.
	//
	if Lsame(uplo, 'U') {
		//
		//        Form  A  when A is stored in upper triangle.
		//
		if incx == 1 {
			for j := 1; j <= n; j++ {
				if (*x)[j-1] != 0 {
					temp = alpha * (*x)[j-1]
					for i := 1; i <= j; i++ {
						(*a)[i-1][j-1] += (*x)[i-1] * temp
					}
				}
			}
		} else {
			jx = kx
			for j := 1; j <= n; j++ {
				if (*x)[jx-1] != 0 {
					temp = alpha * (*x)[jx-1]
					ix = kx
					for i := 1; i <= j; i++ {
						(*a)[i-1][j-1] += (*x)[ix-1] * temp
						ix += incx
					}
				}
				jx += incx
			}
		}
	} else {
		//
		//        Form  A  when A is stored in lower triangle.
		//
		if incx == 1 {
			for j := 1; j <= n; j++ {
				if (*x)[j-1] != 0 {
					temp = alpha * (*x)[j-1]
					for i := j; i <= n; i++ {
						(*a)[i-1][j-1] += (*x)[i-1] * temp
					}
				}
			}
		} else {
			jx = kx
			for j := 1; j <= n; j++ {
				if (*x)[jx-1] != 0 {
					temp = alpha * (*x)[jx-1]
					ix = jx
					for i := j; i <= n; i++ {
						(*a)[i-1][j-1] += (*x)[ix-1] * temp
						ix += incx
					}
				}
				jx += incx
			}
		}
	}

	return
	//
	//     End of Ssyr  .
	//
}
