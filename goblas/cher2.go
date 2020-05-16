package goblas

import (
	"math/cmplx"
)

// Cher2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cher2(uplo,n,alpha,x,incx,y,incy,a,lda)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha
//       INTEGER incx,incy,lda,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX a(lda,*),x(*),y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Cher2  performs the hermitian rank 2 operation
//
//    a := alpha*x*y**H + conjg( alpha )*y*x**H + a,
//
// where alpha is a scalar, x and y are n element vectors and a is an n
// by n hermitian matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the upper or lower
//           triangular part of the array a is to be referenced as
//           follows:
//
//              uplo = 'U' or 'u'   Only the upper triangular part of a
//                                  is to be referenced.
//
//              uplo = 'L' or 'l'   Only the lower triangular part of a
//                                  is to be referenced.
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
//          alpha is COMPLEX
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is COMPLEX array, dimension at least
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
//          y is COMPLEX array, dimension at least
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
// \param[in,out] a
// \verbatim
//          a is COMPLEX array, dimension ( lda, n )
//           Before entry with  uplo = 'U' or 'u', the leading n by n
//           upper triangular part of the array a must contain the upper
//           triangular part of the hermitian matrix and the strictly
//           lower triangular part of a is not referenced. On exit, the
//           upper triangular part of the array a is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry with uplo = 'L' or 'l', the leading n by n
//           lower triangular part of the array a must contain the lower
//           triangular part of the hermitian matrix and the strictly
//           upper triangular part of a is not referenced. On exit, the
//           lower triangular part of the array a is overwritten by the
//           lower triangular part of the updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set, they are assumed to be zero, and on exit they
//           are set to zero.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. lda must be at least
//           max( 1, n ).
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG ld.
//
// \date December 2016
//
// \ingroup complex_blas_level2
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  lvel 2 Blas routine.
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National lb.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National lbs.
// \endverbatim
//
//  =====================================================================
func Cher2(major, uplo *byte, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[][]complex64, lda *int) {
	var temp1, temp2 complex64
	var i, info, ix, iy, j, jx, jy, kx, ky int
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
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
	} else if *lda < max(1, *n) {
		info = 10
	}
	if info != 0 {
		name := "Cher2"
		if common.infoc.test {
			xerblaTest(&name, &info)
			return
		}
		Xerbla(&name, &info)
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
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with one pass through the triangular part
	//     of a.
	//
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  a  when a is stored in the upper triangle.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 || (*y)[j-1] != 0.0 {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*y)[j-1])))
					temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*x)[j-1])))
					for i = 1; i <= j-1; i++ {
						(*a)[i-1][j-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
					}
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0) + complex(real((*x)[j-1]*temp1+(*y)[j-1]*temp2), 0.0)
				} else {
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*y)[jy-1])))
					temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*x)[jx-1])))
					ix = kx
					iy = ky
					for i = 1; i <= j-1; i++ {
						(*a)[i-1][j-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += *incx
						iy += *incy
					}
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0) + complex(real((*x)[jx-1]*temp1+(*y)[jy-1]*temp2), 0.0)
				} else {
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
				}
				jx += *incx
				jy += *incy
			}
		}
	} else {
		//
		//        Form  a  when a is stored in the lower triangle.
		//
		if (*incx == 1) && (*incy == 1) {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 || (*y)[j-1] != 0.0 {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*y)[j-1])))
					temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*x)[j-1])))
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0) + complex(real((*x)[j-1]*temp1+(*y)[j-1]*temp2), 0.0)
					for i = j + 1; i <= *n; i++ {
						(*a)[i-1][j-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
					}
				} else {
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*y)[jy-1])))
					temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*x)[jx-1])))
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0) + complex(real((*x)[jx-1]*temp1+(*y)[jy-1]*temp2), 0.0)
					ix = jx
					iy = jy
					for i = j + 1; i <= *n; i++ {
						ix += *incx
						iy += *incy
						(*a)[i-1][j-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
					}
				} else {
					(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
				}
				jx += *incx
				jy += *incy
			}
		}
	}
}
