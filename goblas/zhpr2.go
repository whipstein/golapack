package goblas

import (
	"math/cmplx"
)

// Zhpr2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 alpha
//       INTEGER incx,incy,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 ap(*),x(*),y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zhpr2  performs the hermitian rank 2 operation
//
//    a := alpha*x*y**H + conjg( alpha )*y*x**H + a,
//
// where alpha is a scalar, x and y are n element vectors and a is an
// n by n hermitian matrix, supplied in packed form.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the upper or lower
//           triangular part of the matrix a is supplied in the packed
//           array ap as follows:
//
//              uplo = 'U' or 'u'   The upper triangular part of a is
//                                  supplied in ap.
//
//              uplo = 'L' or 'L'   The lower triangular part of a is
//                                  supplied in ap.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix a.
//           n must be at least zero.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is COMPLEX*16
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is COMPLEX*16 array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ).
//           Before entry, the incremented array x must contain the n
//           element vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           x. incx must not be zero.
// \endverbatim
//
// \param[in] y
// \verbatim
//          y is COMPLEX*16 array, dimension at least
//           ( 1 + ( n - 1 )*abs( incy ) ).
//           Before entry, the incremented array y must contain the n
//           element vector y.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//           On entry, incy specifies the increment for the elements of
//           y. incy must not be zero.
// \endverbatim
//
// \param[in,out] ap
// \verbatim
//          ap is COMPLEX*16 array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  uplo = 'U' or 'u', the array ap must
//           contain the upper triangular part of the hermitian matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on. On exit, the array
//           ap is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with uplo = 'L' or 'L', the array ap must
//           contain the lower triangular part of the hermitian matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on. On exit, the array
//           ap is overwritten by the lower triangular part of the
//           updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set, they are assumed to be zero, and on exit they
//           are set to zero.
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
// \ingroup complex16_blas_level2
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
func Zhpr2(major, uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, ap *[]complex128) {
	var temp1, temp2 complex128
	var i, info, ix, iy, j, jx, jy, k, kk, kx, ky int
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     Test the input parameters.
	//
	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 2
	} else if *n < 0 {
		info = 3
	} else if *incx == 0 {
		info = 6
	} else if *incy == 0 {
		info = 8
	}
	if info != 0 {
		Xerbla(func() *string { y := "Zhpr2"; return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *n == 0 || *alpha == 0.0 {
		return
	}
	//
	//     Set up the start points in x and y if the increments are not both
	//     unity.
	//
	if *incx != 1 || *incy != 1 {
		if *incx > 0 {
			kx = 1
		} else {
			kx = 1 - ((*n)-1)*(*incx)
		}
		if *incy > 0 {
			ky = 1
		} else {
			ky = 1 - ((*n)-1)*(*incy)
		}
		jx = kx
		jy = ky
	}
	//
	//     Start the operations. In this version the elements of the array ap
	//     are accessed sequentially with one pass through ap.
	//
	kk = 1
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  a  when upper triangle is stored in ap.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 || (*y)[j-1] != 0.0 {
					temp1 = (*alpha) * cmplx.Conj((*y)[j-1])
					temp2 = cmplx.Conj((*alpha) * (*x)[j-1])
					k = kk
					for i = 1; i <= j-1; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
					(*ap)[kk+j-2] = complex(real((*ap)[kk+j-2])+real((*x)[j-1]*temp1+(*y)[j-1]*temp2), 0.0)
				} else {
					(*ap)[kk+j-2] = complex(real((*ap)[kk+j-2]), 0.0)
				}
				kk += j
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * cmplx.Conj((*y)[jy-1])
					temp2 = cmplx.Conj((*alpha) * (*x)[jx-1])
					ix = kx
					iy = ky
					for k = kk; k <= kk+j-2; k++ {
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += *incx
						iy += *incy
					}
					(*ap)[kk+j-2] = complex(real((*ap)[kk+j-2])+real((*x)[jx-1]*temp1+(*y)[jy-1]*temp2), 0.0)
				} else {
					(*ap)[kk+j-2] = complex(real((*ap)[kk+j-2]), 0.0)
				}
				jx += *incx
				jy += *incy
				kk += j
			}
		}
	} else {
		//
		//        Form  a  when lower triangle is stored in ap.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 || (*y)[j-1] != 0.0 {
					temp1 = (*alpha) * cmplx.Conj((*y)[j-1])
					temp2 = cmplx.Conj((*alpha) * (*x)[j-1])
					(*ap)[kk-1] = complex(real((*ap)[kk-1])+real((*x)[j-1]*temp1+(*y)[j-1]*temp2), 0.0)
					k = kk + 1
					for i = j + 1; i <= *n; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
				} else {
					(*ap)[kk-1] = complex(real((*ap)[kk-1]), 0.0)
				}
				kk += (*n) - j + 1
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * cmplx.Conj((*y)[jy-1])
					temp2 = cmplx.Conj((*alpha) * (*x)[jx-1])
					(*ap)[kk-1] = complex(real((*ap)[kk-1])+real((*x)[jx-1]*temp1+(*y)[jy-1]*temp2), 0.0)
					ix = jx
					iy = jy
					for k = kk + 1; k <= kk+(*n)-j; k++ {
						ix += *incx
						iy += *incy
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
					}
				} else {
					(*ap)[kk-1] = complex(real((*ap)[kk-1]), 0.0)
				}
				jx += *incx
				jy += *incy
				kk += (*n) - j + 1
			}
		}
	}
}
