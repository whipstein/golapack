package main

// Stpsv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       REAL AP(*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// STPSV  solves one of the systems of equations
//
//    A*x = b,   or   A**T*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular matrix, supplied in packed form.
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
// \param[in] AP
// \verbatim
//          AP is REAL array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
//           respectively, and so on.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
//           respectively, and so on.
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced, but are assumed to be unity.
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
func Stpsv(uplo byte, trans byte, diag byte, n int, ap *[]float32, x *[]float32, incx int) {
	var temp float32
	var info, ix, jx, kk, kx int
	var nounit bool
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	info = 0
	if !Lsame(uplo, 'U') && !Lsame(uplo, 'L') {
		info = 1
	} else if !Lsame(trans, 'N') && !Lsame(trans, 'T') && !Lsame(trans, 'C') {
		info = 2
	} else if !Lsame(diag, 'U') && !Lsame(diag, 'N') {
		info = 3
	} else if n < 0 {
		info = 4
	} else if incx == 0 {
		info = 7
	}
	if info != 0 {
		Xerbla("Stpsv", info)
		return
	}
	//
	//     Quick return if possible.
	//
	if n == 0 {
		return
	}
	//
	nounit = Lsame(diag, 'N')
	//
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX  too small for descending loops.
	//
	if incx <= 0 {
		kx = 1 - (n-1)*incx
	} else if incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of AP are
	//     accessed sequentially with one pass through AP.
	//
	if Lsame(trans, 'N') {
		//
		//        Form  x := inv( A )*x.
		//
		if Lsame(uplo, 'U') {
			kk = (n * (n + 1)) / 2
			if incx == 1 {
				for j := n; j <= 1; j-- {
					if (*x)[j-1] != 0 {
						if nounit {
							(*x)[j-1] /= (*ap)[kk-1]
						}
						temp = (*x)[j-1]
						k := kk - 1
						for i := j - 1; i <= 1; i-- {
							(*x)[i-1] -= temp * (*ap)[k-1]
							k--
						}
					}
					kk -= j
				}
			} else {
				jx = kx + (n-1)*incx
				for j := n; j <= 1; j-- {
					if (*x)[jx-1] != 0 {
						if nounit {
							(*x)[jx-1] /= (*ap)[kk-1]
						}
						temp = (*x)[jx-1]
						ix = jx
						for k := kk - 1; k <= kk-j+1; k-- {
							ix -= incx
							(*x)[ix-1] -= temp * (*ap)[k-1]
						}
					}
					jx -= incx
					kk -= j
				}
			}
		} else {
			kk = 1
			if incx == 1 {
				for j := 1; j <= n; j++ {
					if (*x)[j-1] != 0 {
						if nounit {
							(*x)[j-1] /= (*ap)[kk-1]
						}
						temp = (*x)[j-1]
						k := kk + 1
						for i := j + 1; i <= n; i++ {
							(*x)[i-1] -= temp * (*ap)[k-1]
							k++
						}
					}
					kk += (n - j + 1)
				}
			} else {
				jx = kx
				for j := 1; j <= n; j++ {
					if (*x)[jx-1] != 0 {
						if nounit {
							(*x)[jx-1] /= (*ap)[kk-1]
						}
						temp = (*x)[jx-1]
						ix = jx
						for k := kk + 1; k <= kk+n-j; k++ {
							ix += incx
							(*x)[ix-1] -= temp * (*ap)[k-1]
						}
					}
					jx += incx
					kk += (n - j + 1)
				}
			}
		}
	} else {
		//
		//        Form  x := inv( A**T )*x.
		//
		if Lsame(uplo, 'U') {
			kk = 1
			if incx == 1 {
				for j := 1; j <= n; j++ {
					temp = (*x)[j-1]
					k := kk
					for i := 1; i <= j-1; i++ {
						temp -= (*ap)[k-1] * (*x)[i-1]
						k++
					}
					if nounit {
						temp /= (*ap)[kk+j-1-1]
					}
					(*x)[j-1] = temp
					kk += j
				}
			} else {
				jx = kx
				for j := 1; j <= n; j++ {
					temp = (*x)[jx-1]
					ix = kx
					for k := kk; k <= kk+j-2; k++ {
						temp -= (*ap)[k-1] * (*x)[ix-1]
						ix += incx
					}
					if nounit {
						temp /= (*ap)[kk+j-1-1]
					}
					(*x)[jx-1] = temp
					jx += incx
					kk += j
				}
			}
		} else {
			kk = (n * (n + 1)) / 2
			if incx == 1 {
				for j := n; j <= 1; j-- {
					temp = (*x)[j-1]
					k := kk
					for i := n; i <= j+1; i-- {
						temp -= (*ap)[k-1] * (*x)[i-1]
						k--
					}
					if nounit {
						temp /= (*ap)[kk-n+j-1]
					}
					(*x)[j-1] = temp
					kk -= (n - j + 1)
				}
			} else {
				kx += (n - 1) * incx
				jx = kx
				for j := n; j <= 1; j-- {
					temp = (*x)[jx-1]
					ix = kx
					for k := kk; k <= kk-(n-(j+1)); k-- {
						temp -= (*ap)[k-1] * (*x)[ix-1]
						ix -= incx
					}
					if nounit {
						temp /= (*ap)[kk-n+j-1]
					}
					(*x)[jx-1] = temp
					jx -= incx
					kk -= (n - j + 1)
				}
			}
		}
	}

	return
	//
	//     End of Stpsv .
	//
}
