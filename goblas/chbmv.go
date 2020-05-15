package goblas

import (
	"math/cmplx"
)

// Chbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Chbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha,beta
//       INTEGER incx,incy,k,lda,n
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
// Chbmv  performs the matrix-vector  operation
//
//    y := alpha*a*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// a is an n by n hermitian band matrix, with k super-diagonals.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the upper or lower
//           triangular part of the band matrix a is being supplied as
//           follows:
//
//              uplo = 'U' or 'u'   The upper triangular part of a is
//                                  being supplied.
//
//              uplo = 'L' or 'l'   The lower triangular part of a is
//                                  being supplied.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix a.
//           n must be at least zero.
// \endverbatim
//
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry, k specifies the number of super-diagonals of the
//           matrix a. k must satisfy  0 .le. k.
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
//           Before entry with uplo = 'U' or 'u', the leading ( k + 1 )
//           by n part of the array a must contain the upper triangular
//           band part of the hermitian matrix, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1 ) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array a is not referenced.
//           The following program segment will transfer the upper
//           triangular part of a hermitian band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, j = 1, n
//                    m = k + 1 - j
//                    DO 10, I = MAX( 1, j - k ), j
//                       a( m + I, j ) = matrix( I, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with uplo = 'L' or 'l', the leading ( k + 1 )
//           by n part of the array a must contain the lower triangular
//           band part of the hermitian matrix, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array a is not referenced.
//           The following program segment will transfer the lower
//           triangular part of a hermitian band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, j = 1, n
//                    m = 1 - j
//                    DO 10, I = j, MIN( n, j + k )
//                       a( m + I, j ) = matrix( I, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Note that the imaginary parts of the diagonal elements need
//           not be set and are assumed to be zero.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. lda must be at least
//           ( k + 1 ).
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ).
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
//           On entry, beta specifies the scalar beta.
// \endverbatim
//
// \param[in,out] y
// \verbatim
//          y is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( incy ) ).
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
//  The vector and matrix arguments are not referenced when n = 0, or m = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National lb.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National lbs.
// \endverbatim
//
//  =====================================================================
func Chbmv(major, uplo *byte, n, k *int, alpha *complex64, a *[][]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	var temp1, temp2 complex64
	var i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l int
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
	} else if *k < 0 {
		info = 4
	} else if *lda < (*k)+1 {
		info = 7
	} else if *incx == 0 {
		info = 9
	} else if *incy == 0 {
		info = 12
	}
	if info != 0 {
		Xerbla(func() *string { y := "Chbmv"; return &y }(), &info)
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
	//     Start the operations. In this version the elements of the array a
	//     are accessed sequentially with one pass through a.
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
		//        Form  y  when upper triangle of a is stored.
		//
		kplus1 = (*k) + 1
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				l = kplus1 - j
				for i = max(1, j-(*k)); i <= j-1; i++ {
					(*y)[i-1] = (*y)[i-1] + temp1*(*a)[l+i-1][j-1]
					temp2 += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
				}
				(*y)[j-1] += temp1*complex(real((*a)[kplus1-1][j-1]), 0.0) + (*alpha)*temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				ix = kx
				iy = ky
				l = kplus1 - j
				for i = max(1, j-(*k)); i <= j-1; i++ {
					(*y)[iy-1] = (*y)[iy-1] + temp1*(*a)[l+i-1][j-1]
					temp2 += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
					ix += *incx
					iy += *incy
				}
				(*y)[jy-1] += temp1*complex(real((*a)[kplus1-1][j-1]), 0.0) + (*alpha)*temp2
				jx += *incx
				jy += *incy
				if j > *k {
					kx += *incx
					ky += *incy
				}
			}
		}
	} else {
		//
		//        Form  y  when lower triangle of a is stored.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				(*y)[j-1] = (*y)[j-1] + temp1*complex(real((*a)[0][j-1]), 0.0)
				l = 1 - j
				for i = j + 1; i <= (min(*n, j+(*k))); i++ {
					(*y)[i-1] = (*y)[i-1] + temp1*(*a)[l+i-1][j-1]
					temp2 += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
				}
				(*y)[j-1] += (*alpha) * temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				(*y)[jy-1] += temp1 * complex(real((*a)[0][j-1]), 0.0)
				l = 1 - j
				ix = jx
				iy = jy
				for i = j + 1; i <= min(*n, j+(*k)); i++ {
					ix += *incx
					iy += *incy
					(*y)[iy-1] += temp1 * (*a)[l+i-1][j-1]
					temp2 += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
				}
				(*y)[jy-1] += (*alpha) * temp2
				jx += *incx
				jy += *incy
			}
		}
	}
}
