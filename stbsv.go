package golapack

// Stbsv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,K,LDA,N
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
// STBSV  solves one of the systems of equations
//
//    A*x = b,   or   A**T*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular band matrix, with ( k + 1 )
// diagonals.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.
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
//           On entry, TRANS specifies the equations to be solved as
//           follows:
//
//              TRANS = 'N' or 'n'   A*x = b.
//
//              TRANS = 'T' or 't'   A**T*x = b.
//
//              TRANS = 'C' or 'c'   A**T*x = b.
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
//          A is REAL array, dimension ( LDA, N )
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
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the n
//           element right-hand side vector b. On exit, X is overwritten
//           with the solution vector x.
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
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Stbsv(uplo *byte, trans *byte, diag *byte, n *int, k *int, a *[]float32, aoff, lda *int, x *[]float32, xoff, incx *int) {
	var nounit bool
	var temp float32
	var i, info, ix, j, jx, kplus1, kx, l int
	var zero float32 = 0.0

	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
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
		Xerbla(func() *[]byte { y := []byte("STBSV "); return &y }(), &info)
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
	//     accessed by sequentially with one pass through A.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  x := inv( A )*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					if (*x)[j-1+(*xoff)] != zero {
						l = kplus1 - j
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] / (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
						temp = (*x)[j-1+(*xoff)]
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] - temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							//Label10:
						}
					}
					//Label20:
				}
			} else {
				kx = kx + ((*n)-1)*(*incx)
				jx = kx
				for j = (*n); j >= 1; j-- {
					kx = kx - (*incx)
					if (*x)[jx-1+(*xoff)] != zero {
						ix = kx
						l = kplus1 - j
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] / (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
						}
						temp = (*x)[jx-1+(*xoff)]
						for i = j - 1; i >= maxint(1, j-(*k)); i-- {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] - temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix - (*incx)
							//Label30:
						}
					}
					jx = jx - (*incx)
					//Label40:
				}
			}
		} else {
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					if (*x)[j-1+(*xoff)] != zero {
						l = 1 - j
						if nounit {
							(*x)[j-1+(*xoff)] = (*x)[j-1+(*xoff)] / (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
						temp = (*x)[j-1+(*xoff)]
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							(*x)[i-1+(*xoff)] = (*x)[i-1+(*xoff)] - temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							//Label50:
						}
					}
					//Label60:
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					kx = kx + (*incx)
					if (*x)[jx-1+(*xoff)] != zero {
						ix = kx
						l = 1 - j
						if nounit {
							(*x)[jx-1+(*xoff)] = (*x)[jx-1+(*xoff)] / (*a)[0+(j-1)*(*lda)+(*aoff)]
						}
						temp = (*x)[jx-1+(*xoff)]
						for i = j + 1; i <= minint(*n, j+(*k)); i++ {
							(*x)[ix-1+(*xoff)] = (*x)[ix-1+(*xoff)] - temp*(*a)[l+i-1+(j-1)*(*lda)+(*aoff)]
							ix = ix + (*incx)
							//Label70:
						}
					}
					jx = jx + (*incx)
					//Label80:
				}
			}
		}
	} else {
		//
		//        Form  x := inv( A**T)*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if (*incx) == 1 {
				for j = 1; j <= (*n); j++ {
					temp = (*x)[j-1+(*xoff)]
					l = kplus1 - j
					for i = maxint(1, j-(*k)); i <= j-1; i++ {
						temp = temp - (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
						//Label90:
					}
					if nounit {
						temp = temp / (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
					}
					(*x)[j-1+(*xoff)] = temp
					//Label100:
				}
			} else {
				jx = kx
				for j = 1; j <= (*n); j++ {
					temp = (*x)[jx-1+(*xoff)]
					ix = kx
					l = kplus1 - j
					for i = maxint(1, j-(*k)); i <= j-1; i++ {
						temp = temp - (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
						ix = ix + (*incx)
						//Label110:
					}
					if nounit {
						temp = temp / (*a)[kplus1-1+(j-1)*(*lda)+(*aoff)]
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx + (*incx)
					if j > (*k) {
						kx = kx + (*incx)
					}
					//Label120:
				}
			}
		} else {
			if (*incx) == 1 {
				for j = (*n); j >= 1; j-- {
					temp = (*x)[j-1+(*xoff)]
					l = 1 - j
					for i = minint(*n, j+(*k)); i >= j+1; i-- {
						temp = temp - (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
						//Label130:
					}
					if nounit {
						temp = temp / (*a)[0+(j-1)*(*lda)+(*aoff)]
					}
					(*x)[j-1+(*xoff)] = temp
					//Label140:
				}
			} else {
				kx = kx + ((*n)-1)*(*incx)
				jx = kx
				for j = (*n); j >= 1; j-- {
					temp = (*x)[jx-1+(*xoff)]
					ix = kx
					l = 1 - j
					for i = minint(*n, j+(*k)); i >= j+1; i-- {
						temp = temp - (*a)[l+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
						ix = ix - (*incx)
						//Label150:
					}
					if nounit {
						temp = temp / (*a)[0+(j-1)*(*lda)+(*aoff)]
					}
					(*x)[jx-1+(*xoff)] = temp
					jx = jx - (*incx)
					if ((*n) - j) >= (*k) {
						kx = kx - (*incx)
					}
					//Label160:
				}
			}
		}
	}
	//
	//
	//     End of STBSV .
	//
}
