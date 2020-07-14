package golapack

// Dspr ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION ALPHA
//       INTEGER INCX,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION AP(*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// DSPR    performs the symmetric rank 1 operation
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
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the matrix A is supplied in the packed
//           array AP as follows:
//
//              UPLO = 'U' or 'u'   The upper triangular part of A is
//                                  supplied in AP.
//
//              UPLO = 'L' or 'l'   The lower triangular part of A is
//                                  supplied in AP.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is DOUBLE PRECISION.
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the n
//           element vector x.
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//           On entry, INCX specifies the increment for the elements of
//           X. INCX must not be zero.
// \endverbatim
//
// \param[in,out] AP
// \verbatim
//          AP is DOUBLE PRECISION array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
//           and a( 2, 2 ) respectively, and so on. On exit, the array
//           AP is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular part of the symmetric matrix
//           packed sequentially, column by column, so that AP( 1 )
//           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
//           and a( 3, 1 ) respectively, and so on. On exit, the array
//           AP is overwritten by the lower triangular part of the
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
// \ingroup double_blas_level2
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
func Dspr(uplo *byte, n *int, alpha *float64, x *[]float64, xoff, incx *int, ap *[]float64, apoff *int) {
	var temp float64
	var i, info, ix, j, jx, k, kk, kx int
	var zero float64 = 0.0

	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 1
	} else if (*n) < 0 {
		info = 2
	} else if (*incx) == 0 {
		info = 5
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("DSPR  "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || ((*alpha) == zero) {
		return
	}
	//
	//     Set the start point in X if the increment is not unity.
	//
	if (*incx) <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if (*incx) != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of the array AP
	//     are accessed sequentially with one pass through AP.
	//
	kk = 1
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  A  when upper triangle is stored in AP.
		//
		if (*incx) == 1 {
			for j = 1; j <= (*n); j++ {
				if (*x)[j-1+(*xoff)] != zero {
					temp = (*alpha) * (*x)[j-1+(*xoff)]
					k = kk
					for i = 1; i <= j; i++ {
						(*ap)[k-1+(*apoff)] = (*ap)[k-1+(*apoff)] + (*x)[i-1+(*xoff)]*temp
						k = k + 1
					}
				}
				kk = kk + j
			}
		} else {
			jx = kx
			for j = 1; j <= (*n); j++ {
				if (*x)[jx-1+(*xoff)] != zero {
					temp = (*alpha) * (*x)[jx-1+(*xoff)]
					ix = kx
					for k = kk; k <= kk+j-1; k++ {
						(*ap)[k-1+(*apoff)] = (*ap)[k-1+(*apoff)] + (*x)[ix-1+(*xoff)]*temp
						ix = ix + (*incx)
					}
				}
				jx = jx + (*incx)
				kk = kk + j
			}
		}
	} else {
		//
		//        Form  A  when lower triangle is stored in AP.
		//
		if (*incx) == 1 {
			for j = 1; j <= (*n); j++ {
				if (*x)[j-1+(*xoff)] != zero {
					temp = (*alpha) * (*x)[j-1+(*xoff)]
					k = kk
					for i = j; i <= (*n); i++ {
						(*ap)[k-1+(*apoff)] = (*ap)[k-1+(*apoff)] + (*x)[i-1+(*xoff)]*temp
						k = k + 1
					}
				}
				kk = kk + (*n) - j + 1
			}
		} else {
			jx = kx
			for j = 1; j <= (*n); j++ {
				if (*x)[jx-1+(*xoff)] != zero {
					temp = (*alpha) * (*x)[jx-1+(*xoff)]
					ix = jx
					for k = kk; k <= kk+(*n)-j; k++ {
						(*ap)[k-1+(*apoff)] = (*ap)[k-1+(*apoff)] + (*x)[ix-1+(*xoff)]*temp
						ix = ix + (*incx)
					}
				}
				jx = jx + (*incx)
				kk = kk + (*n) - j + 1
			}
		}
	}
}
