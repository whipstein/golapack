package goblas

import "math/cmplx"

// Ctbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ctbmv(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,K,LDA,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX A(LDA,*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ctbmv  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular band matrix, with ( k + 1 ) diagonals.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//              UPLO = 'L' or 'l'   A is a lower triangular matrix.
// \endverbatim
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER*1
//           On entry, TRANS specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   x := A*x.
//
//              TRANS = 'T' or 't'   x := A**T*x.
//
//              TRANS = 'C' or 'c'   x := A**H*x.
// \endverbatim
//
// \param[in] DIAG
// \verbatim
//          DIAG is CHARACTER*1
//           On entry, DIAG specifies whether or not A is unit
//           triangular as follows:
//
//              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
//
//              DIAG = 'N' or 'n'   A is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry with UPLO = 'U' or 'u', K specifies the number of
//           super-diagonals of the matrix A.
//           On entry with UPLO = 'L' or 'l', K specifies the number of
//           sub-diagonals of the matrix A.
//           K must satisfy  0 .le. K.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX array, dimension ( LDA, N ).
//           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
//           by n part of the array A must contain the upper triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1 ) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array A is not referenced.
//           The following program segment will transfer an upper
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, N
//                    M = K + 1 - j
//                    DO 10, i = MAX( 1, j - K ), j
//                       A( M + i, j ) = matrix( i, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
//           by n part of the array A must contain the lower triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array A is not referenced.
//           The following program segment will transfer a lower
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, N
//                    M = 1 - j
//                    DO 10, i = j, MIN( N, j + K )
//                       A( M + i, j ) = matrix( i, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Note that when DIAG = 'U' or 'u' the elements of the array A
//           corresponding to the diagonal elements of the matrix are not
//           referenced, but are assumed to be unity.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           ( k + 1 ).
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the n
//           element vector x. On exit, X is overwritten with the
//           transformed vector x.
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//           On entry, INCX specifies the increment for the elements of
//           X. INCX must not be zero.
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
//  The vector and matrix arguments are not referenced when N = 0, or M = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Ctbmv(major, uplo, trans, diag *byte, n, k *int, a *[][]complex64, lda *int, x *[]complex64, incx *int) {
	var temp complex64
	var i, info, ix, j, jx, kplus1, kx, l int
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
	} else if *k < 0 {
		info = 6
	} else if *lda < ((*k) + 1) {
		info = 8
	} else if *incx == 0 {
		info = 10
	}
	if info != 0 {
		Xerbla(func() *string { y := "Ctbmv"; return &y }(), &info)
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
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX   too small for descending loops.
	//
	if *incx <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if *incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through A.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//         Form  x := A*x.
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
		//        Form  x := A**T*x  or  x := A**H*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					temp = (*x)[j-1]
					l = kplus1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[kplus1-1][j-1]
						}
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							temp += (*a)[l+i-1][j-1] * (*x)[i-1]
						}
					} else {
						if nounit {
							temp = temp * complex64(cmplx.Conj(complex128((*a)[kplus1-1][j-1])))
						}
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							temp += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
						}
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
					if noconj {
						if nounit {
							temp = temp * (*a)[kplus1-1][j-1]
						}
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							temp += (*a)[l+i-1][j-1] * (*x)[ix-1]
							ix -= *incx
						}
					} else {
						if nounit {
							temp = temp * complex64(cmplx.Conj(complex128((*a)[kplus1-1][j-1])))
						}
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							temp += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
							ix -= *incx
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
					l = 1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[1-1][j-1]
						}
						for i = j + 1; i <= min(*n, j+(*k)); i++ {
							temp += (*a)[l+i-1][j-1] * (*x)[i-1]
						}
					} else {
						if nounit {
							temp = temp * complex64(cmplx.Conj(complex128((*a)[1-1][j-1])))
						}
						for i = j + 1; i <= (min(*n, j+(*k))); i++ {
							temp += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
						}
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
					if noconj {
						if nounit {
							temp = temp * (*a)[1-1][j-1]
						}
						for i = j + 1; i <= min(*n, j+(*k)); i++ {
							temp += (*a)[l+i-1][j-1] * (*x)[ix-1]
							ix += *incx
						}
					} else {
						if nounit {
							temp = temp * complex64(cmplx.Conj(complex128((*a)[1-1][j-1])))
						}
						for i = j + 1; i <= (min(*n, j+(*k))); i++ {
							temp += complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
							ix += *incx
						}
					}
					(*x)[jx-1] = temp
					jx += *incx
				}
			}
		}
	}
}
