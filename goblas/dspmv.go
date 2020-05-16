package goblas

// Dspmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION alpha,beta
//       INTEGER incx,incy,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION ap(*),x(*),y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dspmv  performs the matrix-vector operation
//
//    y := alpha*a*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// a is an n by n symmetric matrix, supplied in packed form.
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
//              uplo = 'U' or 'U'   The upper triangular part of a is
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
//          alpha is DOUBLE PRECISION.
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] ap
// \verbatim
//          ap is DOUBLE PRECISION array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with uplo = 'U' or 'U', the array ap must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on.
//           Before entry with uplo = 'L' or 'L', the array ap must
//           contain the lower triangular part of the symmetric matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on.
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is DOUBLE PRECISION array, dimension at least
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
//          beta is DOUBLE PRECISION.
//           On entry, beta specifies the scalar beta. When beta is
//           supplied as zero then y need not be set on input.
// \endverbatim
//
// \param[in,out] y
// \verbatim
//          y is DOUBLE PRECISION array, dimension at least
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
// \ingroup double_blas_level2
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
func Dspmv(major, uplo *byte, n *int, alpha *float64, ap *[]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	var temp1, temp2 float64
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
		info = 7
	} else if *incy == 0 {
		info = 10
	}
	if info != 0 {
		name := "Dspmv"
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
	//     Start the operations. In this version the elements of the array ap
	//     are accessed sequentially with one pass through ap.
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
	kk = 1
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  y  when ap contains the upper triangle.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				k = kk
				for i = 1; i <= j-1; i++ {
					(*y)[i-1] += temp1 * (*ap)[k-1]
					temp2 += (*ap)[k-1] * (*x)[i-1]
					k++
				}
				(*y)[j-1] += temp1*(*ap)[kk+j-1-1] + (*alpha)*temp2
				kk += j
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				ix = kx
				iy = ky
				for k = kk; k <= kk+j-2; k++ {
					(*y)[iy-1] += temp1 * (*ap)[k-1]
					temp2 += (*ap)[k-1] * (*x)[ix-1]
					ix += *incx
					iy += *incy
				}
				(*y)[jy-1] += temp1*(*ap)[kk+j-1-1] + (*alpha)*temp2
				jx += *incx
				jy += *incy
				kk += j
			}
		}
	} else {
		//
		//        Form  y  when ap contains the lower triangle.
		//
		if *incx == 1 && *incy == 1 {
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[j-1]
				temp2 = 0.0
				(*y)[j-1] += temp1 * (*ap)[kk-1]
				k = kk + 1
				for i = j + 1; i <= *n; i++ {
					(*y)[i-1] += temp1 * (*ap)[k-1]
					temp2 += (*ap)[k-1] * (*x)[i-1]
					k++
				}
				(*y)[j-1] += (*alpha) * temp2
				kk += ((*n) - j + 1)
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= *n; j++ {
				temp1 = (*alpha) * (*x)[jx-1]
				temp2 = 0.0
				(*y)[jy-1] += temp1 * (*ap)[kk-1]
				ix = jx
				iy = jy
				for k = kk + 1; k <= kk+(*n)-j; k++ {
					ix += *incx
					iy += *incy
					(*y)[iy-1] += temp1 * (*ap)[k-1]
					temp2 += (*ap)[k-1] * (*x)[ix-1]
				}
				(*y)[jy-1] += (*alpha) * temp2
				jx += *incx
				jy += *incy
				kk += ((*n) - j + 1)
			}
		}
	}
}
