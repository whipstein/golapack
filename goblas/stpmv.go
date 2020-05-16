package goblas

// Stpmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE STPMV(uplo,trans,diag,n,ap,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       CHARACTER diag,trans,uplo
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
// STPMV  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix, supplied in packed form.
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
//              uplo = 'U' or 'u'   A is an upper triangular matrix.
//
//              uplo = 'L' or 'l'   A is a lower triangular matrix.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry, trans specifies the operation to be performed as
//           follows:
//
//              trans = 'N' or 'N'   x := A*x.
//
//              trans = 'T' or 't'   x := A**T*x.
//
//              trans = 'C' or 'c'   x := A**T*x.
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//           On entry, diag specifies whether or not A is unit
//           triangular as follows:
//
//              diag = 'U' or 'u'   A is assumed to be unit triangular.
//
//              diag = 'N' or 'N'   A is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix A.
//           n must be at least 0.0.
// \endverbatim
//
// \param[in] ap
// \verbatim
//          ap is REAL array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  uplo = 'U' or 'u', the array ap must
//           contain the upper triangular matrix packed sequentially,
//           column by column, so that ap( 1 ) contains a( 1, 1 ),
//           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
//           respectively, and so on.
//           Before entry with uplo = 'L' or 'l', the array ap must
//           contain the lower triangular matrix packed sequentially,
//           column by column, so that ap( 1 ) contains a( 1, 1 ),
//           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
//           respectively, and so on.
//           Note that when  diag = 'U' or 'u', the diagonal elements of
//           A are not referenced, but are assumed to be unity.
// \endverbatim
//
// \param[in,out] x
// \verbatim
//          x is REAL array, dimension at least
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
//           x. incx must not be 0.0.
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
func Stpmv(major, uplo, trans, diag *byte, n *int, ap *[]float32, x *[]float32, incx *int) {
	var temp float32
	var i, info, ix, j, jx, k, kk, kx int
	var nounit bool
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
	} else if *incx == 0 {
		info = 8
	}
	if info != 0 {
		name := "Stpmv"
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
	//     Start the operations. In this version the elements of ap are
	//     accessed sequentially with one pass through ap.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  x:= A*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kk = 1
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						k = kk
						for i = 1; i <= j-1; i++ {
							(*x)[i-1] += temp * (*ap)[k-1]
							k++
						}
						if nounit {
							(*x)[j-1] *= (*ap)[kk+j-1-1]
						}
					}
					kk += j
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					if (*x)[jx-1] != 0.0 {
						temp = (*x)[jx-1]
						ix = kx
						for k = kk; k <= kk+j-2; k++ {
							(*x)[ix-1] += temp * (*ap)[k-1]
							ix += *incx
						}
						if nounit {
							(*x)[jx-1] *= (*ap)[kk+j-1-1]
						}
					}
					jx += *incx
					kk += j
				}
			}
		} else {
			kk = ((*n) * ((*n) + 1)) / 2
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						k = kk
						for i = *n; i >= j+1; i-- {
							(*x)[i-1] += temp * (*ap)[k-1]
							k--
						}
						if nounit {
							(*x)[j-1] *= (*ap)[kk-(*n)+j-1]
						}
					}
					kk -= ((*n) - j + 1)
				}
			} else {
				kx += ((*n) - 1) * (*incx)
				jx = kx
				for j = *n; j >= 1; j-- {
					if (*x)[jx-1] != 0.0 {
						temp = (*x)[jx-1]
						ix = kx
						for k = kk; k >= kk-((*n)-(j+1)); k-- {
							(*x)[ix-1] += temp * (*ap)[k-1]
							ix -= *incx
						}
						if nounit {
							(*x)[jx-1] *= (*ap)[kk-(*n)+j-1]
						}
					}
					jx -= *incx
					kk -= ((*n) - j + 1)
				}
			}
		}
	} else {
		//
		//        Form  x := A**T*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kk = ((*n) * ((*n) + 1)) / 2
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					temp = (*x)[j-1]
					if nounit {
						temp *= (*ap)[kk-1]
					}
					k = kk - 1
					for i = j - 1; i >= 1; i-- {
						temp += (*ap)[k-1] * (*x)[i-1]
						k--
					}
					(*x)[j-1] = temp
					kk -= j
				}
			} else {
				jx = kx + ((*n)-1)*(*incx)
				for j = *n; j >= 1; j-- {
					temp = (*x)[jx-1]
					ix = jx
					if nounit {
						temp *= (*ap)[kk-1]
					}
					for k = kk - 1; k >= kk-j+1; k-- {
						ix -= *incx
						temp += (*ap)[k-1] * (*x)[ix-1]
					}
					(*x)[jx-1] = temp
					jx -= *incx
					kk -= j
				}
			}
		} else {
			kk = 1
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					temp = (*x)[j-1]
					if nounit {
						temp *= (*ap)[kk-1]
					}
					k = kk + 1
					for i = j + 1; i <= *n; i++ {
						temp += (*ap)[k-1] * (*x)[i-1]
						k++
					}
					(*x)[j-1] = temp
					kk += ((*n) - j + 1)
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					temp = (*x)[jx-1]
					ix = jx
					if nounit {
						temp *= (*ap)[kk-1]
					}
					for k = kk + 1; k <= kk+(*n)-j; k++ {
						ix += *incx
						temp += (*ap)[k-1] * (*x)[ix-1]
					}
					(*x)[jx-1] = temp
					jx += *incx
					kk += ((*n) - j + 1)
				}
			}
		}
	}

	return
}
