package main

// Sspr2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
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
// SSPR2  performs the symmetric rank 2 operation
//
//    A := alpha*x*y**T + alpha*y*x**T + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an
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
// \param[in] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCY ) ).
//           Before entry, the incremented array Y must contain the n
//           element vector y.
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//           On entry, INCY specifies the increment for the elements of
//           Y. INCY must not be zero.
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
func Sspr2(uplo byte, n int, alpha float32, x *[]float32, incx int, y *[]float32, incy int, ap *[]float32) {
	var temp1, temp2 float32
	var info, ix, iy, jx, jy, kk, kx, ky int
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	info = 0
	if !Lsame(uplo, 'U') && !Lsame(uplo, 'L') {
		info = 1
	} else if n < 0 {
		info = 2
	} else if incx == 0 {
		info = 5
	} else if incy == 0 {
		info = 7
	}
	if info != 0 {
		Xerbla("Sspr2", info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) || (alpha == 0) {
		return
	}
	//
	//     Set up the start points in X and Y if the increments are not both
	//     unity.
	//
	if (incx != 1) || (incy != 1) {
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
		jx = kx
		jy = ky
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
		if (incx == 1) && (incy == 1) {
			for j := 1; j <= n; j++ {
				if ((*x)[j-1] != 0) || ((*y)[j-1] != 0) {
					temp1 = alpha * (*y)[j-1]
					temp2 = alpha * (*x)[j-1]
					k := kk
					for i := 1; i <= j; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
				}
				kk += j
			}
		} else {
			for j := 1; j <= n; j++ {
				if ((*x)[jx-1] != 0) || ((*y)[jy-1] != 0) {
					temp1 = alpha * (*y)[jy-1]
					temp2 = alpha * (*x)[jx-1]
					ix = kx
					iy = ky
					for k := kk; k <= kk+j-1; k++ {
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += incx
						iy += incy
					}
				}
				jx += incx
				jy += incy
				kk += j
				//Label40:
			}
		}
	} else {
		//
		//        Form  A  when lower triangle is stored in AP.
		//
		if (incx == 1) && (incy == 1) {
			for j := 1; j <= n; j++ {
				if ((*x)[j-1] != 0) || ((*y)[j-1] != 0) {
					temp1 = alpha * (*y)[j-1]
					temp2 = alpha * (*x)[j-1]
					k := kk
					for i := j; i <= n; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
				}
				kk += n - j + 1
			}
		} else {
			for j := 1; j <= n; j++ {
				if ((*x)[jx-1] != 0) || ((*y)[jy-1] != 0) {
					temp1 = alpha * (*y)[jy-1]
					temp2 = alpha * (*x)[jx-1]
					ix = jx
					iy = jy
					for k := kk; k <= kk+n-j; k++ {
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += incx
						iy += incy
					}
				}
				jx += incx
				jy += incy
				kk += n - j + 1
			}
		}
	}

	return
	//
	//     End of Sspr2 .
	//
}