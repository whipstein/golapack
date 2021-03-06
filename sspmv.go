package golapack

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
func Sspmv(uplo *byte, n *int, alpha *float32, ap *[]float32, apoff *int, x *[]float32, xoff, incx *int, beta *float32, y *[]float32, yoff, incy *int) {
	var temp1, temp2 float32
	var i, info, ix, iy, j, jx, jy, k, kk, kx, ky int
	var one float32 = 1.0
	var zero float32 = 0.0

	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
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
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 1
	} else if (*n) < 0 {
		info = 2
	} else if (*incx) == 0 {
		info = 6
	} else if (*incy) == 0 {
		info = 9
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("SSPMV "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || (((*alpha) == zero) && ((*beta) == one)) {
		return
	}
	//
	//     Set up the start points in  X  and  Y.
	//
	if (*incx) > 0 {
		kx = 1
	} else {
		kx = 1 - ((*n)-1)*(*incx)
	}
	if (*incy) > 0 {
		ky = 1
	} else {
		ky = 1 - ((*n)-1)*(*incy)
	}
	//
	//     Start the operations. In this version the elements of the array AP
	//     are accessed sequentially with one pass through AP.
	//
	//     First form  y := beta*y.
	//
	if (*beta) != one {
		if (*incy) == 1 {
			if (*beta) == zero {
				for i = 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = zero
					//Label10:
				}
			} else {
				for i = 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = (*beta) * (*y)[i-1+(*yoff)]
					//Label20:
				}
			}
		} else {
			iy = ky
			if (*beta) == zero {
				for i = 1; i <= (*n); i++ {
					(*y)[iy-1+(*yoff)] = zero
					iy = iy + (*incy)
					//Label30:
				}
			} else {
				for i = 1; i <= (*n); i++ {
					(*y)[iy-1+(*yoff)] = (*beta) * (*y)[iy-1+(*yoff)]
					iy = iy + (*incy)
					//Label40:
				}
			}
		}
	}
	if (*alpha) == zero {
		return
	}
	kk = 1
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  y  when AP contains the upper triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[j-1+(*xoff)]
				temp2 = zero
				k = kk
				for i = 1; i <= j-1; i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp1*(*ap)[k-1+(*apoff)]
					temp2 = temp2 + (*ap)[k-1+(*apoff)]*(*x)[i-1+(*xoff)]
					k++
					//Label50:
				}
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + temp1*(*ap)[kk+j-1-1+(*apoff)] + (*alpha)*temp2
				kk = kk + j
				//Label60:
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[jx-1+(*xoff)]
				temp2 = zero
				ix = kx
				iy = ky
				for k = kk; k <= kk+j-2; k++ {
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp1*(*ap)[k-1+(*apoff)]
					temp2 = temp2 + (*ap)[k-1+(*apoff)]*(*x)[ix-1+(*xoff)]
					ix = ix + (*incx)
					iy = iy + (*incy)
					//Label70:
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + temp1*(*ap)[kk+j-1-1+(*apoff)] + (*alpha)*temp2
				jx = jx + (*incx)
				jy = jy + (*incy)
				kk = kk + j
				//Label80:
			}
		}
	} else {
		//
		//        Form  y  when AP contains the lower triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[j-1+(*xoff)]
				temp2 = zero
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + temp1*(*ap)[kk-1+(*apoff)]
				k = kk + 1
				for i = j + 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp1*(*ap)[k-1+(*apoff)]
					temp2 = temp2 + (*ap)[k-1+(*apoff)]*(*x)[i-1+(*xoff)]
					k++
					//Label90:
				}
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + (*alpha)*temp2
				kk = kk + ((*n) - j + 1)
				//Label100:
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[jx-1+(*xoff)]
				temp2 = zero
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + temp1*(*ap)[kk-1+(*apoff)]
				ix = jx
				iy = jy
				for k = kk + 1; k <= kk+(*n)-j; k++ {
					ix = ix + (*incx)
					iy = iy + (*incy)
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp1*(*ap)[k-1+(*apoff)]
					temp2 = temp2 + (*ap)[k-1+(*apoff)]*(*x)[ix-1+(*xoff)]
					//Label110:
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp2
				jx = jx + (*incx)
				jy = jy + (*incy)
				kk = kk + ((*n) - j + 1)
				//Label120:
			}
		}
	}
	//
	//
	//     End of SSPMV .
	//
}
