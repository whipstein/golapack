package goblas

// Ssyr ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYR(uplo,n,alpha,x,incx,a,lda)
//
//       .. Scalar Arguments ..
//       REAL alpha
//       INTEGER incx,lda,n
//       CHARACTER uplo
//       ..
//       .. Array Arguments ..
//       REAL a(lda,*),x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSYR   performs the symmetric rank 1 operation
//
//    a := alpha*x*x**T + a,
//
// where alpha is a real scalar, x is an n element vector and a is an
// n by n symmetric matrix.
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
// \param[in,out] a
// \verbatim
//          a is REAL array, dimension ( lda, n )
//           Before entry with  uplo = 'U' or 'u', the leading n by n
//           upper triangular part of the array a must contain the upper
//           triangular part of the symmetric matrix and the strictly
//           lower triangular part of a is not referenced. On exit, the
//           upper triangular part of the array a is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry with uplo = 'L' or 'l', the leading n by n
//           lower triangular part of the array a must contain the lower
//           triangular part of the symmetric matrix and the strictly
//           upper triangular part of a is not referenced. On exit, the
//           lower triangular part of the array a is overwritten by the
//           lower triangular part of the updated matrix.
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
func Ssyr(major, uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, a *[][]float32, lda *int) {
	var temp float32
	var i, info, ix, j, jx, kx int
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
	} else if *lda < max(1, *n) {
		info = 8
	}
	if info != 0 {
		name := "Ssyr"
		if common.infoc.test {
			xerblaTest(name, info)
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
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with one pass through the triangular part
	//     of a.
	//
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  a  when a is stored in upper triangle.
		//
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 {
					temp = (*alpha) * (*x)[j-1]
					for i = 1; i <= j; i++ {
						(*a)[i-1][j-1] += (*x)[i-1] * temp
					}
				}
			}
		} else {
			jx = kx
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 {
					temp = (*alpha) * (*x)[jx-1]
					ix = kx
					for i = 1; i <= j; i++ {
						(*a)[i-1][j-1] += (*x)[ix-1] * temp
						ix += *incx
					}
				}
				jx += *incx
			}
		}
	} else {
		//
		//        Form  a  when a is stored in lower triangle.
		//
		if *incx == 1 {
			for j = 1; j <= *n; j++ {
				if (*x)[j-1] != 0.0 {
					temp = (*alpha) * (*x)[j-1]
					for i = j; i <= *n; i++ {
						(*a)[i-1][j-1] += (*x)[i-1] * temp
					}
				}
			}
		} else {
			jx = kx
			for j = 1; j <= *n; j++ {
				if (*x)[jx-1] != 0.0 {
					temp = (*alpha) * (*x)[jx-1]
					ix = jx
					for i = j; i <= *n; i++ {
						(*a)[i-1][j-1] += (*x)[ix-1] * temp
						ix += *incx
					}
				}
				jx += *incx
			}
		}
	}

	return
}
