package goblas

// Sgemv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SGEMV(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       REAL alpha,beta
//       INTEGER incx,incy,lda,m,n
//       CHARACTER trans
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
// SGEMV  performs 1.0 of the matrix-vector operations
//
//    y := alpha*a*x + beta*y,   or   y := alpha*a**T*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and a is an
// m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry, trans specifies the operation to be performed as
//           follows:
//
//              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
//
//              trans = 'T' or 't'   y := alpha*a**T*x + beta*y.
//
//              trans = 'C' or 'c'   y := alpha*a**T*x + beta*y.
// \endverbatim
//
// \param[in] m
// \verbatim
//          m is INTEGER
//           On entry, m specifies the number of rows of the matrix a.
//           m must be at least 0.0.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the number of columns of the matrix a.
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
//           Before entry, the leading m by n part of the array a must
//           contain the matrix of coefficients.
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
// \param[in] x
// \verbatim
//          x is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
//           and at least
//           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
//           Before entry, the incremented array x must contain the
//           vector x.
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
//           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
//           and at least
//           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
//           Before entry with beta non-0.0, the incremented array y
//           must contain the vector y. On exit, y is overwritten by the
//           updated vector y.
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
func Sgemv(major, trans *byte, m, n *int, alpha *float32, a *[][]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	var temp float32
	var i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny int
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
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 2
	} else if *m < 0 {
		info = 3
	} else if *n < 0 {
		info = 4
	} else if *lda < max(1, *m) {
		info = 7
	} else if *incx == 0 {
		info = 9
	} else if *incy == 0 {
		info = 12
	}
	if info != 0 {
		name := "Sgemv"
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
	if (*m == 0) || (*n == 0) || ((*alpha == 0.0) && (*beta == 1.0)) {
		return
	}
	//
	//     Set  lenx  and  leny, the lengths of the vectors x and y, and set
	//     up the start points in  x  and  y.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		lenx = *n
		leny = *m
	} else {
		lenx = *m
		leny = *n
	}
	if *incx > 0 {
		kx = 1
	} else {
		kx = 1 - (lenx-1)*(*incx)
	}
	if *incy > 0 {
		ky = 1
	} else {
		ky = 1 - (leny-1)*(*incy)
	}
	//
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with 1.0 pass through a.
	//
	//     First form  y := beta*y.
	//
	if *beta != 1.0 {
		if *incy == 1 {
			if *beta == 0.0 {
				for i = 1; i <= leny; i++ {
					(*y)[i-1] = 0.0
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[i-1] *= *beta
				}
			}
		} else {
			iy = ky
			if *beta == 0.0 {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1] = 0.0
					iy += *incy
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1] *= *beta
					iy += *incy
				}
			}
		}
	}
	if *alpha == 0.0 {
		return
	}
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  y := alpha*a*x + y.
		//
		jx = kx
		if *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp = (*alpha) * (*x)[jx-1]
				for i = 1; i <= *m; i++ {
					(*y)[i-1] += temp * (*a)[i-1][j-1]
				}
				jx += *incx
			}
		} else {
			for j = 1; j <= *n; j++ {
				temp = (*alpha) * (*x)[jx-1]
				iy = ky
				for i = 1; i <= *m; i++ {
					(*y)[iy-1] += temp * (*a)[i-1][j-1]
					iy += *incy
				}
				jx += *incx
			}
		}
	} else {
		//
		//        Form  y := alpha*a**T*x + y.
		//
		jy = ky
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				temp = 0.0
				for i = 1; i <= *m; i++ {
					temp += (*a)[i-1][j-1] * (*x)[i-1]
				}
				(*y)[jy-1] += (*alpha) * temp
				jy += *incy
			}
		} else {
			for j = 1; j <= *n; j++ {
				temp = 0.0
				ix = kx
				for i = 1; i <= *m; i++ {
					temp += (*a)[i-1][j-1] * (*x)[ix-1]
					ix += *incx
				}
				(*y)[jy-1] += (*alpha) * temp
				jy += *incy
			}
		}
	}

	return
}
