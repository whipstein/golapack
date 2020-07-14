package golapack

// Strmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,LDA,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// STRMV  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix.
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
//              TRANS = 'C' or 'c'   x := A**T*x.
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
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, N )
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular matrix and the strictly lower triangular part of
//           A is not referenced.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular matrix and the strictly upper triangular part of
//           A is not referenced.
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced either, but are assumed to be unity.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           max( 1, n ).
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is REAL array, dimension at least
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
// \ingroup single_blas_level2
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
func Strmv(uplo *byte, trans *byte, diag *byte, n *int, a *[]float32, aoff, lda *int, x *[]float32, xoff, incx *int) {
	var nounit bool
	var temp float32
	var i, info, ix, j, jx, kx int
	var zero float32 = 0.0

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
	} else if (*lda) < maxint(1, *n) {
		info = 6
	} else if (*incx) == 0 {
		info = 8
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("STRMV "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (*n) == 0 {
		return
	}
	//
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	//
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX  too small for descending loops.
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
		//        Form  x := A*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					if (*x)[j-1+(*xoff)] != zero {
						temp = (*x)[j-1+(*xoff)]
						for i = 1; i <= j-1; i++ {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
						}
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
						}
					}
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					if (*x)[jx-1+(*xoff)] != zero {
						temp = (*x)[jx-1+(*xoff)]
						ix = kx
						for i = 1; i <= j-1; i++ {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix + (*incx)
						}
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
						}
					}
					jx = jx + (*incx)
				}
			}
		} else {
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					if (*x)[j-1+(*xoff)] != zero {
						temp = (*x)[j-1+(*xoff)]
						for i = (*n); i >= j+1; i-- {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
						}
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
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
						for i = (*n); i >= j+1; i-- {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix - (*incx)
						}
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
						}
					}
					jx = jx - (*incx)
				}
			}
		}
	} else {
		//
		//        Form  x := A**T*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					temp = (*x)[j-1+(*xoff)]
					if nounit {
						temp = temp * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
					}
					for i = j - 1; i >= 1; i-- {
						temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
					}
					(*x)[j-1+(*xoff)] = temp
				}
			} else {
				jx = kx + ((*n)-1)*(*incx)
				for j = (*n); j >= 1; j-- {
					temp = (*x)[jx-1+(*xoff)]
					ix = jx
					if nounit {
						temp = temp * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
					}
					for i = j - 1; i >= 1; i-- {
						ix = ix - (*incx)
						temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx - (*incx)
				}
			}
		} else {
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					temp = (*x)[j-1+(*xoff)]
					if nounit {
						temp = temp * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
					}
					for i = j + 1; i <= (*n); i++ {
						temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
					}
					(*x)[j-1+(*xoff)] = temp
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					temp = (*x)[jx-1+(*xoff)]
					ix = jx
					if nounit {
						temp = temp * (*a)[j-1+(j-1)*(*lda)+(*aoff)]
					}
					for i = j + 1; i <= (*n); i++ {
						ix = ix + (*incx)
						temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx + (*incx)
				}
			}
		}
	}
}
