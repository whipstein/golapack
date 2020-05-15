package main

// Sspmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER INCX,INCY,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       REAL AP(*),X(*),Y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSPMV  performs the matrix-vector operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n symmetric matrix, supplied in packed form.
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
// \param[in] AP
// \verbatim
//          AP is REAL array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on.
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
// \param[in] BETA
// \verbatim
//          BETA is REAL
//           On entry, BETA specifies the scalar beta. When BETA is
//           supplied as zero then Y need not be set on input.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCY ) ).
//           Before entry, the incremented array Y must contain the n
//           element vector y. On exit, Y is overwritten by the updated
//           vector y.
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
func Sspmv(uplo byte, n int, alpha float32, ap *[]float32, x *[]float32, incx int, beta float32, y *[]float32, incy int) {
	var temp1, temp2 float32
	var info, ix, iy, jx, jy, k, kk, kx, ky int
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
		info = 6
	} else if incy == 0 {
		info = 9
	}
	if info != 0 {
		Xerbla("Sspmv", info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) || (alpha == 0 && beta == 1) {
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
	//     Start the operations. In this version the elements of the array AP
	//     are accessed sequentially with one pass through AP.
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
	kk = 1
	if Lsame(uplo, 'U') {
		//
		//        Form  y  when AP contains the upper triangle.
		//
		if (incx == 1) && (incy == 1) {
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[j]
				temp2 = 0
				k = kk
				for i := 0; i < j-1; i++ {
					(*y)[i] += temp1 * (*ap)[k-1]
					temp2 = temp2 + (*ap)[k-1]*(*x)[i]
					k++
				}
				(*y)[j] += temp1*(*ap)[kk+j-1] + alpha*temp2
				kk += j
			}
		} else {
			jx = kx
			jy = ky
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[jx-1]
				temp2 = 0
				ix = kx
				iy = ky
				for k := kk - 1; k < kk+j-2; k++ {
					(*y)[iy-1] += temp1 * (*ap)[k]
					temp2 += (*ap)[k] * (*x)[ix-1]
					ix += incx
					iy += incy
				}
				(*y)[jy-1] += temp1*(*ap)[kk+j-1] + alpha*temp2
				jx += incx
				jy += incy
				kk += j
			}
		}
	} else {
		//
		//        Form  y  when AP contains the lower triangle.
		//
		if (incx == 1) && (incy == 1) {
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[j]
				temp2 = 0
				(*y)[j] += temp1 * (*ap)[kk-1]
				k = kk + 1
				for i := j; i < n; i++ {
					(*y)[i] += temp1 * (*ap)[k-1]
					temp2 += (*ap)[k-1] * (*x)[i]
					k++
				}
				(*y)[j] += alpha * temp2
				kk += (n - j + 2)
			}
		} else {
			jx = kx
			jy = ky
			for j := 0; j < n; j++ {
				temp1 = alpha * (*x)[jx-1]
				temp2 = 0
				(*y)[jy-1] += temp1 * (*ap)[kk-1]
				ix = jx
				iy = jy
				for k := kk; k < kk+n-j; k++ {
					ix += incx
					iy += incy
					(*y)[iy-1] += temp1 * (*ap)[k]
					temp2 += (*ap)[k-1] * (*x)[ix-1]
				}
				(*y)[jy-1] += alpha * temp2
				jx += incx
				jy += incy
				kk += (n - j + 1)
			}
		}
	}
	//
	return
	//
	//     End of Sspmv .
	//
}
