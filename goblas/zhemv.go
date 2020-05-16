package goblas

import (
	"math/cmplx"
)

// Zhemv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 alpha,beta
//       INTEGER incx,incy,lda,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 a(lda,*),x(*),y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zhemv  performs the matrix-vector  operation
//
//    y := alpha*a*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// a is an n by n hermitian matrix.
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
//              uplo = 'L' or 'L'   Only the lower triangular part of a
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
//          alpha is COMPLEX*16
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX*16 array, dimension ( lda, n )
//           Before entry with  uplo = 'U' or 'u', the leading n by n
//           upper triangular part of the array a must contain the upper
//           triangular part of the hermitian matrix and the strictly
//           lower triangular part of a is not referenced.
//           Before entry with uplo = 'L' or 'L', the leading n by n
//           lower triangular part of the array a must contain the lower
//           triangular part of the hermitian matrix and the strictly
//           upper triangular part of a is not referenced.
//           Note that the imaginary parts of the diagonal elements need
//           not be set and are assumed to be zero.
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
// \param[in] beta
// \verbatim
//          beta is COMPLEX*16
//           On entry, beta specifies the scalar beta. When beta is
//           supplied as zero then y need not be set on input.
// \endverbatim
//
// \param[in,out] y
// \verbatim
//          y is COMPLEX*16 array, dimension at least
//           ( 1 + ( n - 1 )*abs( incy ) ).
//           Before entry, the incremented array y must contain the n
//           element vector y. On exit, y is overwritten by the updated
//           vector y.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//           On entry, incy specifies the increment for the elements of
//           y. incy must not be zero.
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
//  The vector and matrix arguments are not referenced when n = 0, or m = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Zhemv(major, uplo *byte, n *int, alpha *complex128, a *[][]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	var temp1, temp2 complex128
	var i, info, ix, iy, j, jx, jy, kx, ky int
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
	} else if *lda < max(1, *n) {
		info = 6
	} else if *incx == 0 {
		info = 8
	} else if *incy == 0 {
		info = 11
	}
	if info != 0 {
		name := "Zhemv"
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
	if *n == 0 || (*alpha == 0.0 && *beta == 1.0) {
		return
	}
	//
	//     Set up the start points in  x  and  y.
	//
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
	//
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with one pass through the triangular part
	//     of a.
	//
	//     First form  y := beta*y.
	//
	if *beta != 1.0 {
		if *incy == 1 {
			if *beta == 0.0 {
				for i = 1; i <= *n; i++ {
					(*y)[i-1] = 0.0
				}
			} else {
				for i = 1; i <= *n; i++ {
					(*y)[i-1] = (*beta) * (*y)[i-1]
				}
			}
		} else {
			iy = ky
			if *beta == 0.0 {
				for i = 1; i <= *n; i++ {
					(*y)[iy-1] = 0.0
					iy += *incy
				}
			} else {
				for i = 1; i <= *n; i++ {
					(*y)[iy-1] = (*beta) * (*y)[iy-1]
					iy += *incy
				}
			}
		}
	}
	if *alpha == 0.0 {
		return
	}
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  y  when a is stored in upper triangle.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				for i = 1; i <= j-1; i++ {
					(*y)[i-1] += temp1 * (*a)[i-1][j-1]
					temp2 += cmplx.Conj((*a)[i-1][j-1]) * (*x)[i-1]
				}
				(*y)[j-1] += temp1*complex(real((*a)[j-1][j-1]), 0.0) + (*alpha)*temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				ix = kx
				iy = ky
				for i = 1; i <= j-1; i++ {
					(*y)[iy-1] += temp1 * (*a)[i-1][j-1]
					temp2 += cmplx.Conj((*a)[i-1][j-1]) * (*x)[ix-1]
					ix += *incx
					iy += *incy
				}
				(*y)[jy-1] += temp1*complex(real((*a)[j-1][j-1]), 0.0) + (*alpha)*temp2
				jx += *incx
				jy += *incy
			}
		}
	} else {
		//
		//        Form  y  when a is stored in lower triangle.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				(*y)[j-1] += temp1 * complex(real((*a)[j-1][j-1]), 0.0)
				for i = j + 1; i <= *n; i++ {
					(*y)[i-1] += temp1 * (*a)[i-1][j-1]
					temp2 += cmplx.Conj((*a)[i-1][j-1]) * (*x)[i-1]
				}
				(*y)[j-1] += (*alpha) * temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				(*y)[jy-1] += temp1 * complex(real((*a)[j-1][j-1]), 0.0)
				ix = jx
				iy = jy
				for i = j + 1; i <= *n; i++ {
					ix += *incx
					iy += *incy
					(*y)[iy-1] += temp1 * (*a)[i-1][j-1]
					temp2 += cmplx.Conj((*a)[i-1][j-1]) * (*x)[ix-1]
				}
				(*y)[jy-1] += (*alpha) * temp2
				jx += *incx
				jy += *incy
			}
		}
	}
}
