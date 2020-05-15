package goblas

import (
	"math/cmplx"
)

// Zgerc ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zgerc(m,n,alpha,x,incx,y,incy,a,lda)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 alpha
//       INTEGER incx,incy,lda,m,n
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
// Zgerc  performs the rank 1 operation
//
//    a := alpha*x*y**H + a,
//
// where alpha is a scalar, x is an m element vector, y is an n element
// vector and a is an m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] m
// \verbatim
//          m is INTEGER
//           On entry, m specifies the number of rows of the matrix a.
//           m must be at least zero.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the number of columns of the matrix a.
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
//           ( 1 + ( m - 1 )*abs( incx ) ).
//           Before entry, the incremented array x must contain the m
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
// \param[in,out] a
// \verbatim
//          a is COMPLEX*16 array, dimension ( lda, n )
//           Before entry, the leading m by n part of the array a must
//           contain the matrix of coefficients. On exit, a is
//           overwritten by the updated matrix.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. lda must be at least
//           max( 1, m ).
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
func Zgerc(major *byte, m, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[][]complex128, lda *int) {
	var temp complex128
	var i, info, ix, j, jy, kx int
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
	} else if *m < 0 {
		info = 2
	} else if *n < 0 {
		info = 3
	} else if *incx == 0 {
		info = 6
	} else if *incy == 0 {
		info = 8
	} else if *lda < max(1, *m) {
		info = 10
	}
	if info != 0 {
		Xerbla(func() *string { y := "Zgerc"; return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *m == 0 || *n == 0 || *alpha == 0.0 {
		return
	}
	//
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with one pass through a.
	//
	if *incy > 0 {
		jy = 1
	} else {
		jy = 1 - ((*n)-1)*(*incy)
	}
	if *incx == 1 {
		for j = 1; j <= *n; j++ {
			if (*y)[jy-1] != 0.0 {
				temp = (*alpha) * cmplx.Conj((*y)[jy-1])
				for i = 1; i <= *m; i++ {
					(*a)[i-1][j-1] += (*x)[i-1] * temp
				}
			}
			jy += *incy
		}
	} else {
		if *incx > 0 {
			kx = 1
		} else {
			kx = 1 - ((*m)-1)*(*incx)
		}
		for j = 1; j <= *n; j++ {
			if (*y)[jy-1] != 0.0 {
				temp = (*alpha) * cmplx.Conj((*y)[jy-1])
				ix = kx
				for i = 1; i <= *m; i++ {
					(*a)[i-1][j-1] += (*x)[ix-1] * temp
					ix += *incx
				}
			}
			jy += *incy
		}
	}
}
