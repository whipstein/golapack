package goblas

// Dtbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Dtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,k,lda,n
//       CHARACTER diag,trans,uplo
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION a(lda,*),x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtbmv  performs one of the matrix-vector operations
//
//    x := a*x,   or   x := a**T*x,
//
// where x is an n element vector and  a is an n by n unit, or non-unit,
// upper or lower triangular band matrix, with ( k + 1 ) diagonals.
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
//              uplo = 'U' or 'U'   a is an upper triangular matrix.
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
//              trans = 'T' or 'T'   x := a**T*x.
//
//              trans = 'C' or 'C'   x := a**T*x.
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//           On entry, diag specifies whether or not a is unit
//           triangular as follows:
//
//              diag = 'U' or 'U'   a is assumed to be unit triangular.
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
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry with uplo = 'U' or 'U', k specifies the number of
//           super-diagonals of the matrix a.
//           On entry with uplo = 'L' or 'L', k specifies the number of
//           sub-diagonals of the matrix a.
//           k must satisfy  0 .le. k.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is DOUBLE PRECISION array, dimension ( lda, n )
//           Before entry with uplo = 'U' or 'U', the leading ( k + 1 )
//           by n part of the array a must contain the upper triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1 ) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array a is not referenced.
//           The following program segment will transfer an upper
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, n
//                    m = k + 1 - j
//                    DO 10, i = MAX( 1, j - k ), j
//                       a( m + i, j ) = matrix( i, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with uplo = 'L' or 'L', the leading ( k + 1 )
//           by n part of the array a must contain the lower triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array a is not referenced.
//           The following program segment will transfer a lower
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, n
//                    m = 1 - j
//                    DO 10, i = j, MIN( n, j + k )
//                       a( m + i, j ) = matrix( i, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Note that when diag = 'U' or 'U' the elements of the array a
//           corresponding to the diagonal elements of the matrix are not
//           referenced, but are assumed to be unity.
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
// \param[in,out] x
// \verbatim
//          x is DOUBLE PRECISION array, dimension at least
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
func Dtbmv(major, uplo, trans, diag *byte, n, k *int, a *[][]float64, lda *int, x *[]float64, incx *int) {
	var temp float64
	var i, info, ix, j, jx, kplus1, kx, l int
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
	} else if *k < 0 {
		info = 6
	} else if *lda < (*k)+1 {
		info = 8
	} else if *incx == 0 {
		info = 10
	}
	if info != 0 {
		name := "Dtbmv"
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
	//     will be  ( n - 1 )*incx   too small for descending loops.
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
		//         Form  x := a*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						l = kplus1 - j
						for i = max(1, j-(*k)); i <= j-1; i++ {
							(*x)[i-1] += temp * (*a)[l+i-1][j-1]
						}
						if nounit {
							(*x)[j-1] = (*x)[j-1] * (*a)[kplus1-1][j-1]
						}
					}
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					if (*x)[jx-1] != 0.0 {
						temp = (*x)[jx-1]
						ix = kx
						l = kplus1 - j
						for i = max(1, j-(*k)); i <= j-1; i++ {
							(*x)[ix-1] += temp * (*a)[l+i-1][j-1]
							ix += *incx
						}
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] * (*a)[kplus1-1][j-1]
						}
					}
					jx += *incx
					if j > *k {
						kx += *incx
					}
				}
			}
		} else {
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					if (*x)[j-1] != 0.0 {
						temp = (*x)[j-1]
						l = 1 - j
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							(*x)[i-1] += temp * (*a)[l+i-1][j-1]
						}
						if nounit {
							(*x)[j-1] = (*x)[j-1] * (*a)[1-1][j-1]
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
						l = 1 - j
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							(*x)[ix-1] += temp * (*a)[l+i-1][j-1]
							ix -= *incx
						}
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] * (*a)[1-1][j-1]
						}
					}
					jx -= *incx
					if ((*n) - j) >= *k {
						kx -= *incx
					}
				}
			}
		}
	} else {
		//
		//        Form  x := a**T*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					temp = (*x)[j-1]
					l = kplus1 - j
					if nounit {
						temp = temp * (*a)[kplus1-1][j-1]
					}
					for i = j - 1; i >= max(1, j-(*k)); i-- {
						temp += (*a)[l+i-1][j-1] * (*x)[i-1]
					}
					(*x)[j-1] = temp
				}
			} else {
				kx += ((*n) - 1) * (*incx)
				jx = kx
				for j = *n; j >= 1; j-- {
					temp = (*x)[jx-1]
					kx -= *incx
					ix = kx
					l = kplus1 - j
					if nounit {
						temp = temp * (*a)[kplus1-1][j-1]
					}
					for i = j - 1; i >= max(1, j-(*k)); i-- {
						temp += (*a)[l+i-1][j-1] * (*x)[ix-1]
						ix -= *incx
					}
					(*x)[jx-1] = temp
					jx -= *incx
				}
			}
		} else {
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					temp = (*x)[j-1]
					l = 1 - j
					if nounit {
						temp = temp * (*a)[1-1][j-1]
					}
					for i = j + 1; i <= min(*n, j+(*k)); i++ {
						temp += (*a)[l+i-1][j-1] * (*x)[i-1]
					}
					(*x)[j-1] = temp
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					temp = (*x)[jx-1]
					kx += *incx
					ix = kx
					l = 1 - j
					if nounit {
						temp = temp * (*a)[1-1][j-1]
					}
					for i = j + 1; i <= min(*n, j+(*k)); i++ {
						temp += (*a)[l+i-1][j-1] * (*x)[ix-1]
						ix += *incx
					}
					(*x)[jx-1] = temp
					jx += *incx
				}
			}
		}
	}
}
