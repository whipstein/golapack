package goblas

// Ssymv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYMV(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       REAL alpha,beta
//       INTEGER incx,incy,lda,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       REAL a(lda,*),x(*),y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSYMV  performs the matrix-vector  operation
//
//    y := alpha*a*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// a is an n by n symmetric matrix.
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
//           n must be at least 0.0.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is REAL
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is REAL array, dimension ( lda, n )
//           Before entry with  uplo = 'U' or 'u', the leading n by n
//           upper triangular part of the array a must contain the upper
//           triangular part of the symmetric matrix and the strictly
//           lower triangular part of a is not referenced.
//           Before entry with uplo = 'L' or 'l', the leading n by n
//           lower triangular part of the array a must contain the lower
//           triangular part of the symmetric matrix and the strictly
//           upper triangular part of a is not referenced.
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
//          x is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ).
//           Before entry, the incremented array x must contain the n
//           element vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           x. incx must not be 0.0.
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is REAL
//           On entry, beta specifies the scalar beta. When beta is
//           supplied as 0.0 then y need not be set on input.
// \endverbatim
//
// \param[in,out] y
// \verbatim
//          y is REAL array, dimension at least
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
//           y. incy must not be 0.0.
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
//  The vector and matrix arguments are not referenced when n = 0, or M = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Ssymv(major, uplo *byte, n *int, alpha *float32, a *[][]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	var temp1, temp2 float32
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
		name := "Ssymv"
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
	if (*n == 0) || ((*alpha == 0.0) && (*beta == 1.0)) {
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
	//     accessed sequentially with 1.0 pass through the triangular part
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
					(*y)[i-1] *= *beta
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
					(*y)[iy-1] *= *beta
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
		if (*incx == 1) && (*incy == 1) {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				for i = 1; i <= j-1; i++ {
					(*y)[i-1] += temp1 * (*a)[i-1][j-1]
					temp2 += (*a)[i-1][j-1] * (*x)[i-1]
				}
				(*y)[j-1] += temp1*(*a)[j-1][j-1] + (*alpha)*temp2
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
					temp2 += (*a)[i-1][j-1] * (*x)[ix-1]
					ix += *incx
					iy += *incy
				}
				(*y)[jy-1] += temp1*(*a)[j-1][j-1] + (*alpha)*temp2
				jx += *incx
				jy += *incy
			}
		}
	} else {
		//
		//        Form  y  when a is stored in lower triangle.
		//
		if (*incx == 1) && (*incy == 1) {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				(*y)[j-1] += temp1 * (*a)[j-1][j-1]
				for i = j + 1; i <= *n; i++ {
					(*y)[i-1] += temp1 * (*a)[i-1][j-1]
					temp2 += (*a)[i-1][j-1] * (*x)[i-1]
				}
				(*y)[j-1] += (*alpha) * temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				(*y)[jy-1] += temp1 * (*a)[j-1][j-1]
				ix = jx
				iy = jy
				for i = j + 1; i <= *n; i++ {
					ix += *incx
					iy += *incy
					(*y)[iy-1] += temp1 * (*a)[i-1][j-1]
					temp2 += (*a)[i-1][j-1] * (*x)[ix-1]
				}
				(*y)[jy-1] += (*alpha) * temp2
				jx += *incx
				jy += *incy
			}
		}
	}

	return
}
