package main

// Ssbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER INCX,INCY,K,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,*),X(*),Y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSBMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n symmetric band matrix, with k super-diagonals.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the band matrix A is being supplied as
//           follows:
//
//              UPLO = 'U' or 'u'   The upper triangular part of A is
//                                  being supplied.
//
//              UPLO = 'L' or 'l'   The lower triangular part of A is
//                                  being supplied.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry, K specifies the number of super-diagonals of the
//           matrix A. K must satisfy  0 .le. K.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, N )
//           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
//           by n part of the array A must contain the upper triangular
//           band part of the symmetric matrix, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1 ) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array A is not referenced.
//           The following program segment will transfer the upper
//           triangular part of a symmetric band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, J = 1, N
//                    M = K + 1 - J
//                    DO 10, I = MAX( 1, J - K ), J
//                       A( M + I, J ) = matrix( I, J )
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
//           by n part of the array A must contain the lower triangular
//           band part of the symmetric matrix, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array A is not referenced.
//           The following program segment will transfer the lower
//           triangular part of a symmetric band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, J = 1, N
//                    M = 1 - J
//                    DO 10, I = J, MIN( N, J + K )
//                       A( M + I, J ) = matrix( I, J )
//              10    CONTINUE
//              20 CONTINUE
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           ( k + 1 ).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the
//           vector x.
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//           On entry, INCX specifies the increment for the elements of
//           X. INCX must not be zero.
// \endverbatim
//
// \param[in] BETA
// \verbatim
//          BETA is REAL
//           On entry, BETA specifies the scalar beta.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCY ) ).
//           Before entry, the incremented array Y must contain the
//           vector y. On exit, Y is overwritten by the updated vector y.
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//           On entry, INCY specifies the increment for the elements of
//           Y. INCY must not be zero.
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
//  The vector and matrix arguments are not referenced when N = 0, or M = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Ssbmv(uplo byte, n int, k int, alpha float32, a *[][]float32, lda int, x *[]float32, incx int, beta float32, y *[]float32, incy int) {
	var temp1, temp2 float32
	var info, ix, iy, jx, jy, kplus1, kx, ky, l int
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
	} else if k < 0 {
		info = 3
	} else if lda < (k + 1) {
		info = 6
	} else if incx == 0 {
		info = 8
	} else if incy == 0 {
		info = 11
	}
	if info != 0 {
		Xerbla("Ssbmv", info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) || ((alpha == 0) && (beta == 1)) {
		return
	}
	//
	//     Set up the start points in  X  and  Y.
	//
	if incx > 0 {
		kx = 1
	} else {
		kx = 1 - (n-1)*incx
	}
	if incy > 0 {
		ky = 1
	} else {
		ky = 1 - (n-1)*incy
	}
	//
	//     Start the operations. In this version the elements of the array A
	//     are accessed sequentially with one pass through A.
	//
	//     First form  y := beta*y.
	//
	if beta != 1 {
		if incy == 1 {
			if beta == 0 {
				for i := 0; i < n; i++ {
					(*y)[i] = 0
				}
			} else {
				for i := 0; i < n; i++ {
					(*y)[i] *= beta
				}
			}
		} else {
			iy = ky
			if beta == 0 {
				for i := 0; i < n; i++ {
					(*y)[iy-1] = 0
					iy += incy
				}
			} else {
				for i := 0; i < n; i++ {
					(*y)[iy-1] *= beta
					iy += incy
				}
			}
		}
	}
	if alpha == 0 {
		return
	}
	if Lsame(uplo, 'U') {
		//
		//        Form  y  when upper triangle of A is stored.
		//
		kplus1 = k + 1
		if (incx == 1) && (incy == 1) {
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[j]
				temp2 = 0
				l = kplus1 - j
				for i := max(0, j-k); i < j-1; i++ {
					(*y)[i] += temp1 * (*a)[l+i][j]
					temp2 += (*a)[l+i][j] * (*x)[i]
				}
				(*y)[j] += temp1*(*a)[kplus1-1][j] + alpha*temp2
			}
		} else {
			jx = kx
			jy = ky
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[jx-1]
				temp2 = 0
				ix = kx
				iy = ky
				l = kplus1 - j
				for i := max(0, j-k); i < j-1; i++ {
					(*y)[iy-1] += temp1 * (*a)[l+i][j]
					temp2 += (*a)[l+i][j] * (*x)[ix-1]
					ix += incx
					iy += incy
				}
				(*y)[jy-1] += temp1*(*a)[kplus1-1][j] + alpha*temp2
				jx += incx
				jy += incy
				if j > k {
					kx += incx
					ky += incy
				}
			}
		}
	} else {
		//
		//        Form  y  when lower triangle of A is stored.
		//
		if (incx == 1) && (incy == 1) {
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[j]
				temp2 = 0
				(*y)[j] += temp1 * (*a)[0][j]
				l = 1 - j
				for i := j; i < min(n, j+k); i++ {
					(*y)[i] += temp1 * (*a)[l+i][j]
					temp2 += (*a)[l+i][j] * (*x)[i]
				}
				(*y)[j] += alpha * temp2
			}
		} else {
			jx = kx
			jy = ky
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[jx-1]
				temp2 = 0
				(*y)[jy-1] += temp1 * (*a)[0][j]
				l = 1 - j
				ix = jx
				iy = jy
				for i := j; i < min(n, j+k); i++ {
					ix += incx
					iy += incy
					(*y)[iy-1] += temp1 * (*a)[l+i][j]
					temp2 += (*a)[l+i][j] * (*x)[ix-1]
				}
				(*y)[jy-1] += alpha * temp2
				jx += incx
				jy += incy
			}
		}
	}
	//
	return
	//
	//     End of SSBMV .
	//
}
