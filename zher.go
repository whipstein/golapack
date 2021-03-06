package golapack

// Zher ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION ALPHA
//       INTEGER INCX,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 A(LDA,*),X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// ZHER   performs the hermitian rank 1 operation
//
//    A := alpha*x*x**H + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
// n by n hermitian matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the array A is to be referenced as
//           follows:
//
//              UPLO = 'U' or 'u'   Only the upper triangular part of A
//                                  is to be referenced.
//
//              UPLO = 'L' or 'l'   Only the lower triangular part of A
//                                  is to be referenced.
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
//          X is COMPLEX*16 array, dimension at least
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
// \param[in,out] A
// \verbatim
//          A is COMPLEX*16 array, dimension ( LDA, N )
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular part of the hermitian matrix and the strictly
//           lower triangular part of A is not referenced. On exit, the
//           upper triangular part of the array A is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular part of the hermitian matrix and the strictly
//           upper triangular part of A is not referenced. On exit, the
//           lower triangular part of the array A is overwritten by the
//           lower triangular part of the updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set, they are assumed to be zero, and on exit they
//           are set to zero.
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
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Zher(uplo *byte, n *int, alpha *float64, x *[]complex128, xoff, incx *int, a *[]complex128, aoff, lda *int) {
	var temp complex128
	var i, info, ix, j, jx, kx int
	var zero complex128 = (0.0 + 0.0*1i)
	var calpha complex128 = complex(*alpha, 0.0)

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
	} else if (*lda) < maxint(1, *n) {
		info = 7
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("ZHER  "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || ((*alpha) == real(zero)) {
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
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through the triangular part
	//     of A.
	//
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  A  when A is stored in upper triangle.
		//
		if (*incx) == 1 {
			for j = 1; j <= (*n); j++ {
				if (*x)[j-1+(*xoff)] != zero {
					temp = calpha * conjc128((*x)[j-1+(*xoff)])
					for i = 1; i <= j-1; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp
					}
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[j-1+(*xoff)]*temp), 0.0)
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
			}
		} else {
			jx = kx
			for j = 1; j <= (*n); j++ {
				if (*x)[jx-1+(*xoff)] != zero {
					temp = calpha * conjc128((*x)[jx-1+(*xoff)])
					ix = kx
					for i = 1; i <= j-1; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp
						ix = ix + (*incx)
					}
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[jx-1+(*xoff)]*temp), 0.0)
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
				jx = jx + (*incx)
			}
		}
	} else {
		//
		//        Form  A  when A is stored in lower triangle.
		//
		if (*incx) == 1 {
			for j = 1; j <= (*n); j++ {
				if (*x)[j-1+(*xoff)] != zero {
					temp = calpha * conjc128((*x)[j-1+(*xoff)])
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real(temp*(*x)[j-1+(*xoff)]), 0.0)
					for i = j + 1; i <= (*n); i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp
					}
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
			}
		} else {
			jx = kx
			for j = 1; j <= (*n); j++ {
				if (*x)[jx-1+(*xoff)] != zero {
					temp = calpha * conjc128((*x)[jx-1+(*xoff)])
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real(temp*(*x)[jx-1+(*xoff)]), 0.0)
					ix = jx
					for i = j + 1; i <= (*n); i++ {
						ix = ix + (*incx)
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp
					}
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
				jx = jx + (*incx)
			}
		}
	}
}
