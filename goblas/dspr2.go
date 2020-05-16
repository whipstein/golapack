package goblas

// Dspr2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Dspr2(uplo,n,alpha,x,incx,y,incy,ap)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION alpha
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
// Dspr2  performs the symmetric rank 2 operation
//
//    a := alpha*x*y**T + alpha*y*x**T + a,
//
// where alpha is a scalar, x and y are n element vectors and a is an
// n by n symmetric matrix, supplied in packed form.
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
//              uplo = 'l' or 'l'   The lower triangular part of a is
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
// \param[in] y
// \verbatim
//          y is DOUBLE PRECISION array, dimension at least
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
//          ap is DOUBLE PRECISION array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  uplo = 'U' or 'u', the array ap must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on. On exit, the array
//           ap is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with uplo = 'l' or 'l', the array ap must
//           contain the lower triangular part of the symmetric matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on. On exit, the array
//           ap is overwritten by the lower triangular part of the
//           updated matrix.
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
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Dspr2(major, uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, ap *[]float64) {
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
	} else if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('l'); return &y }()) {
		info = 2
	} else if *n < 0 {
		info = 3
	} else if *incx == 0 {
		info = 6
	} else if *incy == 0 {
		info = 8
	}
	if info != 0 {
		name := "Dspr2"
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
					temp1 = (*alpha) * (*y)[j-1]
					temp2 = (*alpha) * (*x)[j-1]
					k = kk
					for i = 1; i <= j; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
				}
				kk += j
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * (*y)[jy-1]
					temp2 = (*alpha) * (*x)[jx-1]
					ix = kx
					iy = ky
					for k = kk; k <= kk+j-1; k++ {
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += *incx
						iy += *incy
					}
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
					temp1 = (*alpha) * (*y)[j-1]
					temp2 = (*alpha) * (*x)[j-1]
					k = kk
					for i = j; i <= *n; i++ {
						(*ap)[k-1] += (*x)[i-1]*temp1 + (*y)[i-1]*temp2
						k++
					}
				}
				kk += (*n) - j + 1
			}
		} else {
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 || (*y)[jy-1] != 0.0 {
					temp1 = (*alpha) * (*y)[jy-1]
					temp2 = (*alpha) * (*x)[jx-1]
					ix = jx
					iy = jy
					for k = kk; k <= kk+(*n)-j; k++ {
						(*ap)[k-1] += (*x)[ix-1]*temp1 + (*y)[iy-1]*temp2
						ix += *incx
						iy += *incy
					}
				}
				jx += *incx
				jy += *incy
				kk += (*n) - j + 1
			}
		}
	}
}
