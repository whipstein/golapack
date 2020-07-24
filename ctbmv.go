package golapack

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
//       SUBROUTINE CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
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
// CTBMV  performs one of the matrix-vector operations
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
//                 DO 20, J = 1, N
//                    M = K + 1 - J
//                    DO 10, I = MAX( 1, J - K ), J
//                       A( M + I, J ) = matrix( I, J )
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
//                 DO 20, J = 1, N
//                    M = 1 - J
//                    DO 10, I = J, MIN( N, J + K )
//                       A( M + I, J ) = matrix( I, J )
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
func Ctbmv(uplo *byte, trans *byte, diag *byte, n *int, k *int, a *[]complex64, aoff, lda *int, x *[]complex64, xoff, incx *int) {
	var noconj, nounit bool
	var temp complex64
	var i, info, ix, j, jx, kplus1, kx, l int
	var zero complex64 = (0.0 + 0.0*1i)

	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 1
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 2
	} else if !Lsame(diag, func() *byte { y := byte('U'); return &y }()) && !Lsame(diag, func() *byte { y := byte('N'); return &y }()) {
		info = 3
	} else if (*n) < 0 {
		info = 4
	} else if (*k) < 0 {
		info = 5
	} else if (*lda) < ((*k) + 1) {
		info = 7
	} else if (*incx) == 0 {
		info = 9
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("CTBMV "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (*n) == 0 {
		return
	}

	noconj = Lsame(trans, func() *byte { y := byte('T'); return &y }())
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	//
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX   too small for descending loops.
	//
	if (*incx) <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if (*incx) != 1 {
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
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					if (*x)[j-1+(*xoff)] != zero {
						temp = (*x)[j-1+(*xoff)]
						l = kplus1 - j
						for i = maxint(1, j-(*k)); i <= j-1; i++ {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] + temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
						}
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] * (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
					}
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					if (*x)[jx-1+(*xoff)] != zero {
						temp = (*x)[jx-1+(*xoff)]
						ix = kx
						l = kplus1 - j
						for i = maxint(1, j-(*k)); i <= j-1; i++ {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] + temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix + (*incx)
						}
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] * (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
					}
					jx = jx + (*incx)
					if j > (*k) {
						kx = kx + (*incx)
					}
				}
			}
		} else {
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					if (*x)[j-1+(*xoff)] != zero {
						temp = (*x)[j-1+(*xoff)]
						l = 1 - j
						for i = minint(*n, j+(*k)); i >= j+1; i-- {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] + temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
						}
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] * (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
					}
				}
			} else {
				kx = kx + ((*n)-1)*(*incx)
				jx = kx
				for j = (*n); j >= 1; j-- {
					if (*x)[jx-1+(*xoff)] != zero {
						temp = (*x)[jx-1+(*xoff)]
						ix = kx
						l = 1 - j
						for i = minint(*n, j+(*k)); i >= j+1; i-- {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] + temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix - (*incx)
						}
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] * (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
					}
					jx = jx - (*incx)
					if ((*n) - j) >= (*k) {
						kx = kx - (*incx)
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
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					temp = (*x)[j-1+(*xoff)]
					l = kplus1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							temp = temp + (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
						}
					} else {
						if nounit {
							temp = temp * conjc64((*a)[kplus1-1+(j-1)*(*lda)+(*aoff)])
						}
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							temp = temp + conjc64((*a)[l+i-1+(j-1)*(*lda)+(*aoff)])*(*x)[i-1+(*xoff)]
						}
					}
					(*x)[j-1+(*xoff)] = temp
				}
			} else {
				kx = kx + ((*n)-1)*(*incx)
				jx = kx
				for j = (*n); j >= 1; j-- {
					temp = (*x)[jx-1+(*xoff)]
					kx = kx - (*incx)
					ix = kx
					l = kplus1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							temp = temp + (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
							ix = ix - (*incx)
						}
					} else {
						if nounit {
							temp = temp * conjc64((*a)[kplus1-1+(j-1)*(*lda)+(*aoff)])
						}
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							temp = temp + conjc64((*a)[l+i-1+(j-1)*(*lda)+(*aoff)])*(*x)[ix-1+(*xoff)]
							ix = ix - (*incx)
						}
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx - (*incx)
				}
			}
		} else {
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					temp = (*x)[j-1+(*xoff)]
					l = 1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							temp = temp + (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
						}
					} else {
						if nounit {
							temp = temp * conjc64((*a)[0+(j-1)*(*lda)+(*aoff)])
						}
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							temp = temp + conjc64((*a)[l+i-1+(j-1)*(*lda)+(*aoff)])*(*x)[i-1+(*xoff)]
						}
					}
					(*x)[j-1+(*xoff)] = temp
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					temp = (*x)[jx-1+(*xoff)]
					kx = kx + (*incx)
					ix = kx
					l = 1 - j
					if noconj {
						if nounit {
							temp = temp * (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							temp = temp + (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
							ix = ix + (*incx)
						}
					} else {
						if nounit {
							temp = temp * conjc64((*a)[0+(j-1)*(*lda)+(*aoff)])
						}
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							temp = temp + conjc64((*a)[l+i-1+(j-1)*(*lda)+(*aoff)])*(*x)[ix-1+(*xoff)]
							ix = ix + (*incx)
						}
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx + (*incx)
				}
			}
		}
	}
}
