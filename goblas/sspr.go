package goblas

// Sspr ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSPR(uplo,n,alpha,x,incx,ap)
//
//       .. Scalar Arguments ..
//       REAL alpha
//       INTEGER incx,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       REAL ap(*),x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSPR    performs the symmetric rank 1 operation
//
//    A := alpha*x*x**T + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
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
//           triangular part of the matrix A is supplied in the packed
//           array ap as follows:
//
//              uplo = 'U' or 'u'   The upper triangular part of A is
//                                  supplied in ap.
//
//              uplo = 'L' or 'l'   The lower triangular part of A is
//                                  supplied in ap.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix A.
//           n must be at least 0.0.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is REAL
//           On entry, alpha specifies the scalar alpha.
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
// \param[in,out] ap
// \verbatim
//          ap is REAL array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  uplo = 'U' or 'u', the array ap must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that ap( 1 )
//           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on. On exit, the array
//           ap is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with uplo = 'L' or 'l', the array ap must
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
// \ingroup single_blas_level2
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
func Sspr(major, uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, ap *[]float32) {
	var temp float32
	var i, info, ix, j, jx, k, kk, kx int
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
	}
	if info != 0 {
		name := "Sspr"
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
	if (*n == 0) || (*alpha == 0.0) {
		return
	}
	//
	//     Set the start point in x if the increment is not unity.
	//
	if *incx <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if *incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of the array ap
	//     are accessed sequentially with one pass through ap.
	//
	kk = 1
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  A  when upper triangle is stored in ap.
		//
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 {
					temp = (*alpha) * (*x)[j-1]
					k = kk
					for i = 1; i <= j; i++ {
						(*ap)[k-1] += (*x)[i-1] * temp
						k++
					}
				}
				kk += j
			}
		} else {
			jx = kx
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 {
					temp = (*alpha) * (*x)[jx-1]
					ix = kx
					for k = kk; k <= kk+j-1; k++ {
						(*ap)[k-1] += (*x)[ix-1] * temp
						ix += *incx
					}
				}
				jx += *incx
				kk += j
			}
		}
	} else {
		//
		//        Form  A  when lower triangle is stored in ap.
		//
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 {
					temp = (*alpha) * (*x)[j-1]
					k = kk
					for i = j; i <= *n; i++ {
						(*ap)[k-1] += (*x)[i-1] * temp
						k++
					}
				}
				kk += (*n) - j + 1
			}
		} else {
			jx = kx
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 {
					temp = (*alpha) * (*x)[jx-1]
					ix = jx
					for k = kk; k <= kk+(*n)-j; k++ {
						(*ap)[k-1] += (*x)[ix-1] * temp
						ix += *incx
					}
				}
				jx += *incx
				kk += (*n) - j + 1
			}
		}
	}

	return
}
