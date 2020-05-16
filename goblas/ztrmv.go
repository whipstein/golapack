package goblas

import (
	"math/cmplx"
)

// Ztrmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ztrmv(uplo,trans,diag,n,a,lda,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,lda,n
//       CHARACTER diag,trans,uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 a(lda,*),x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ztrmv  performs one of the matrix-vector operations
//
//    x := a*x,   or   x := a**T*x,   or   x := a**H*x,
//
// where x is an n element vector and  a is an n by n unit, or non-unit,
// upper or lower triangular matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              uplo = 'U' or 'u'   a is an upper triangular matrix.
//
//              uplo = 'L' or 'L'   a is a lower triangular matrix.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry, trans specifies the operation to be performed as
//           follows:
//
//              trans = 'N' or 'N'   x := a*x.
//
//              trans = 'T' or 't'   x := a**T*x.
//
//              trans = 'C' or 'C'   x := a**H*x.
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//           On entry, diag specifies whether or not a is unit
//           triangular as follows:
//
//              diag = 'U' or 'u'   a is assumed to be unit triangular.
//
//              diag = 'N' or 'N'   a is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix a.
//           n must be at least zero.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX*16 array, dimension ( lda, n ).
//           Before entry with  uplo = 'U' or 'u', the leading n by n
//           upper triangular part of the array a must contain the upper
//           triangular matrix and the strictly lower triangular part of
//           a is not referenced.
//           Before entry with uplo = 'L' or 'L', the leading n by n
//           lower triangular part of the array a must contain the lower
//           triangular matrix and the strictly upper triangular part of
//           a is not referenced.
//           Note that when  diag = 'U' or 'u', the diagonal elements of
//           a are not referenced either, but are assumed to be unity.
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
// \param[in,out] x
// \verbatim
//          x is COMPLEX*16 array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ).
//           Before entry, the incremented array x must contain the n
//           element vector x. On exit, x is overwritten with the
//           transformed vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           x. incx must not be zero.
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
func Ztrmv(major, uplo, trans, diag *byte, n *int, a *[][]complex128, lda *int, x *[]complex128, incx *int) {
	var temp complex128
	var i, info, ix, j, jx, kx int
	var noconj, nounit bool
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
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 3
	} else if !Lsame(diag, func() *byte { y := byte('U'); return &y }()) && !Lsame(diag, func() *byte { y := byte('N'); return &y }()) {
		info = 4
	} else if *n < 0 {
		info = 5
	} else if *lda < max(1, *n) {
		info = 7
	} else if *incx == 0 {
		info = 9
	}
	if info != 0 {
		name := "Ztrmv"
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
	if *n == 0 {
		return
	}
	//
	noconj = Lsame(trans, func() *byte { y := byte('T'); return &y }())
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	//
	//     Set up the start point in x if the increment is not unity. This
	//     will be  ( n - 1 )*incx  too small for descending loops.
	//
	if *incx <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if *incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of a are
	//     accessed sequentially with one pass through a.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  x := a*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						for i = 1; i <= j-1; i++ {
							(*x)[i-1] += temp * (*a)[i-1][j-1]
						}
						if nounit {
							(*x)[j-1] = (*x)[j-1] * (*a)[j-1][j-1]
						}
					}
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					if (*x)[jx-1] != 0.0 {
						temp = (*x)[jx-1]
						ix = kx
						for i = 1; i <= j-1; i++ {
							(*x)[ix-1] += temp * (*a)[i-1][j-1]
							ix += *incx
						}
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] * (*a)[j-1][j-1]
						}
					}
					jx += *incx
				}
			}
		} else {
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						for i = *n; i >= j+1; i-- {
							(*x)[i-1] += temp * (*a)[i-1][j-1]
						}
						if nounit {
							(*x)[j-1] = (*x)[j-1] * (*a)[j-1][j-1]
						}
					}
				}
			} else {
				kx += ((*n) - 1) * (*incx)
				jx = kx
				for j = *n; j >= 1; j-- {
					if (*x)[jx-1] != 0.0 {
						temp = (*x)[jx-1]
						ix = kx
						for i = *n; i >= j+1; i-- {
							(*x)[ix-1] += temp * (*a)[i-1][j-1]
							ix -= *incx
						}
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] * (*a)[j-1][j-1]
						}
					}
					jx -= *incx
				}
			}
		}
	} else {
		//
		//        Form  x := a**T*x  or  x := a**H*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					temp = (*x)[j-1]
					if noconj {
						if nounit {
							temp = temp * (*a)[j-1][j-1]
						}
						for i = j - 1; i >= 1; i-- {
							temp += (*a)[i-1][j-1] * (*x)[i-1]
						}
					} else {
						if nounit {
							temp = temp * cmplx.Conj((*a)[j-1][j-1])
						}
						for i = j - 1; i >= 1; i-- {
							temp += cmplx.Conj((*a)[i-1][j-1]) * (*x)[i-1]
						}
					}
					(*x)[j-1] = temp
				}
			} else {
				jx = kx + ((*n)-1)*(*incx)
				for j = *n; j >= 1; j-- {
					temp = (*x)[jx-1]
					ix = jx
					if noconj {
						if nounit {
							temp = temp * (*a)[j-1][j-1]
						}
						for i = j - 1; i >= 1; i-- {
							ix -= *incx
							temp += (*a)[i-1][j-1] * (*x)[ix-1]
						}
					} else {
						if nounit {
							temp = temp * cmplx.Conj((*a)[j-1][j-1])
						}
						for i = j - 1; i >= 1; i-- {
							ix -= *incx
							temp += cmplx.Conj((*a)[i-1][j-1]) * (*x)[ix-1]
						}
					}
					(*x)[jx-1] = temp
					jx -= *incx
				}
			}
		} else {
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					temp = (*x)[j-1]
					if noconj {
						if nounit {
							temp = temp * (*a)[j-1][j-1]
						}
						for i = j + 1; i <= *n; i++ {
							temp += (*a)[i-1][j-1] * (*x)[i-1]
						}
					} else {
						if nounit {
							temp = temp * cmplx.Conj((*a)[j-1][j-1])
						}
						for i = j + 1; i <= *n; i++ {
							temp += cmplx.Conj((*a)[i-1][j-1]) * (*x)[i-1]
						}
					}
					(*x)[j-1] = temp
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					temp = (*x)[jx-1]
					ix = jx
					if noconj {
						if nounit {
							temp = temp * (*a)[j-1][j-1]
						}
						for i = j + 1; i <= *n; i++ {
							ix += *incx
							temp += (*a)[i-1][j-1] * (*x)[ix-1]
						}
					} else {
						if nounit {
							temp = temp * cmplx.Conj((*a)[j-1][j-1])
						}
						for i = j + 1; i <= *n; i++ {
							ix += *incx
							temp += cmplx.Conj((*a)[i-1][j-1]) * (*x)[ix-1]
						}
					}
					(*x)[jx-1] = temp
					jx += *incx
				}
			}
		}
	}
}
