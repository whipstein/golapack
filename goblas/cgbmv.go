package goblas

import (
	"math/cmplx"
)

// Cgbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha,beta
//       INTEGER incx,incy,kl,ku,lda,m,n
//       CHARACTER trans
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
// Cgbmv  performs one of the matrix-vector operations
//
//    y := alpha*a*x + beta*y,   or   y := alpha*a**T*x + beta*y,   or
//
//    y := alpha*a**H*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and a is an
// m by n band matrix, with kl sub-diagonals and ku super-diagonals.
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
//              trans = 'N' or 'n'   y := alpha*a*x + beta*y.
//
//              trans = 'T' or 't'   y := alpha*a**T*x + beta*y.
//
//              trans = 'C' or 'c'   y := alpha*a**H*x + beta*y.
// \endverbatim
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
// \param[in] kl
// \verbatim
//          kl is INTEGER
//           On entry, kl specifies the number of sub-diagonals of the
//           matrix a. kl must satisfy  0 .le. kl.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is INTEGER
//           On entry, ku specifies the number of super-diagonals of the
//           matrix a. ku must satisfy  0 .le. ku.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is COMPLEX
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX array, dimension ( lda, n )
//           Before entry, the leading ( kl + ku + 1 ) by n part of the
//           array a must contain the matrix of coefficients, supplied
//           column by column, with the leading diagonal of the matrix in
//           row ( ku + 1 ) of the array, the first super-diagonal
//           starting at position 2 in row ku, the first sub-diagonal
//           starting at position 1 in row ( ku + 2 ), and so on.
//           Elements in the array a that do not correspond to elements
//           in the band matrix (such as the top left ku by ku triangle)
//           are not referenced.
//           The following program segment will transfer a band matrix
//           from conventional full matrix storage to band storage:
//
//                 DO 20, j = 1, n
//                    k = ku + 1 - j
//                    DO 10, I = MAX( 1, j - ku ), MIN( m, j + kl )
//                       a( k + I, j ) = matrix( I, j )
//              10    CONTINUE
//              20 CONTINUE
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. lda must be at least
//           ( kl + ku + 1 ).
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'N' or 'n'
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
//           x. incx must not be zero.
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is COMPLEX
//           On entry, beta specifies the scalar beta. When beta is
//           supplied as zero then y need not be set on input.
// \endverbatim
//
// \param[in,out] y
// \verbatim
//          y is COMPLEX array, dimension at least
//           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'N' or 'n'
//           and at least
//           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
//           Before entry, the incremented array y must contain the
//           vector y. On exit, y is overwritten by the updated vector y.
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
// \ingroup complex_blas_level2
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
func Cgbmv(major, trans *byte, m, n, kl, ku *int, alpha *complex64, a *[][]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	var temp complex64
	var i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx, leny int
	var noconj bool
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
	} else if *kl < 0 {
		info = 5
	} else if *ku < 0 {
		info = 6
	} else if *lda < (*kl)+(*ku)+1 {
		info = 9
	} else if *incx == 0 {
		info = 11
	} else if *incy == 0 {
		info = 14
	}
	if info != 0 {
		Xerbla(func() *string { y := "Cgbmv"; return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *m == 0 || *n == 0 || (*alpha == 0.0 && *beta == 1.0) {
		return
	}
	//
	noconj = Lsame(trans, func() *byte { y := byte('T'); return &y }())
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
	//     accessed sequentially with one pass through the band part of a.
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
					(*y)[i-1] = (*beta) * (*y)[i-1]
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
					(*y)[iy-1] = (*beta) * (*y)[iy-1]
					iy += *incy
				}
			}
		}
	}
	if *alpha == 0.0 {
		return
	}
	kup1 = (*ku) + 1
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  y := alpha*a*x + y.
		//
		jx = kx
		if *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp = (*alpha) * (*x)[jx-1]
				k = kup1 - j
				for i = max(1, j-(*ku)); i <= min(*m, j+(*kl)); i++ {
					(*y)[i-1] += temp * (*a)[k+i-1][j-1]
				}
				jx += *incx
			}
		} else {
			for j = 1; j <= *n; j++ {
				temp = (*alpha) * (*x)[jx-1]
				iy = ky
				k = kup1 - j
				for i = max(1, j-(*ku)); i <= min(*m, j+(*kl)); i++ {
					(*y)[iy-1] += temp * (*a)[k+i-1][j-1]
					iy += *incy
				}
				jx += *incx
				if j > *ku {
					ky += *incy
				}
			}
		}
	} else {
		//
		//        Form  y := alpha*a**T*x + y  or  y := alpha*a**H*x + y.
		//
		jy = ky
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				temp = 0.0
				k = kup1 - j
				if noconj {
					for i = max(1, j-(*ku)); i <= min(*m, j+(*kl)); i++ {
						temp += (*a)[k+i-1][j-1] * (*x)[i-1]
					}
				} else {
					for i = max(1, j-(*ku)); i <= min(*m, j+(*kl)); i++ {
						temp += complex64(cmplx.Conj(complex128((*a)[k+i-1][j-1]))) * (*x)[i-1]
					}
				}
				(*y)[jy-1] += (*alpha) * temp
				jy += *incy
			}
		} else {
			for j = 1; j <= *n; j++ {
				temp = 0.0
				ix = kx
				k = kup1 - j
				if noconj {
					for i = max(1, j-(*ku)); i <= min(*m, j+(*kl)); i++ {
						temp += (*a)[k+i-1][j-1] * (*x)[ix-1]
						ix += *incx
					}
				} else {
					for i = max(1, j-(*ku)); i <= (min(*m, j+(*kl))); i++ {
						temp += complex64(cmplx.Conj(complex128((*a)[k+i-1][j-1]))) * (*x)[ix-1]
						ix += *incx
					}
				}
				(*y)[jy-1] += (*alpha) * temp
				jy += *incy
				if j > *ku {
					kx += *incx
				}
			}
		}
	}
}
