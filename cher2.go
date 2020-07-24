package golapack

// Cher2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
//
//       .. Scalar Arguments ..
//       COMPLEX ALPHA
//       INTEGER INCX,INCY,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX A(LDA,*),X(*),Y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// CHER2  performs the hermitian rank 2 operation
//
//    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an n
// by n hermitian matrix.
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
//          ALPHA is COMPLEX
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is COMPLEX array, dimension at least
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
// \param[in] Y
// \verbatim
//          Y is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCY ) ).
//           Before entry, the incremented array Y must contain the n
//           element vector y.
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//           On entry, INCY specifies the increment for the elements of
//           Y. INCY must not be zero.
// \endverbatim
//
// \param[in,out] A
// \verbatim
//          A is COMPLEX array, dimension ( LDA, N )
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
// \ingroup complex_blas_level2
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
func Cher2(uplo *byte, n *int, alpha *complex64, x *[]complex64, xoff, incx *int, y *[]complex64, yoff, incy *int, a *[]complex64, aoff, lda *int) {
	var temp1, temp2 complex64
	var i, info, ix, iy, j, jx, jy, kx, ky int
	var zero complex64 = (0.0 + 0.0*1i)

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
	} else if (*incy) == 0 {
		info = 7
	} else if (*lda) < maxint(1, *n) {
		info = 9
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("CHER2 "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || ((*alpha) == zero) {
		return
	}
	//
	//     Set up the start points in X and Y if the increments are not both
	//     unity.
	//
	if ((*incx) != 1) || ((*incy) != 1) {
		if (*incx) > 0 {
			kx = 1
		} else {
			kx = 1 - ((*n)-1)*(*incx)
		}
		if (*incy) > 0 {
			ky = 1
		} else {
			ky = 1 - ((*n)-1)*(*incy)
		}
		jx = kx
		jy = ky
	}
	//
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through the triangular part
	//     of A.
	//
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  A  when A is stored in the upper triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				if ((*x)[j-1+(*xoff)] != zero) || ((*y)[j-1+(*yoff)] != zero) {
					temp1 = (*alpha) * conjc64((*y)[j-1+(*yoff)])
					temp2 = conjc64((*alpha) * (*x)[j-1+(*xoff)])
					for i = 1; i <= j-1; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp1 + (*y)[i-1+(*yoff)]*temp2
					}
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[j-1+(*xoff)]*temp1+(*y)[j-1+(*yoff)]*temp2), 0.0)
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				if ((*x)[jx-1+(*xoff)] != zero) || ((*y)[jy-1+(*yoff)] != zero) {
					temp1 = (*alpha) * conjc64((*y)[jy-1+(*yoff)])
					temp2 = conjc64((*alpha) * (*x)[jx-1+(*xoff)])
					ix = kx
					iy = ky
					for i = 1; i <= j-1; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp1 + (*y)[iy-1+(*yoff)]*temp2
						ix = ix + (*incx)
						iy = iy + (*incy)
					}
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[jx-1+(*xoff)]*temp1+(*y)[jy-1+(*yoff)]*temp2), 0.0)
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
				jx = jx + (*incx)
				jy = jy + (*incy)
			}
		}
	} else {
		//
		//        Form  A  when A is stored in the lower triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				if ((*x)[j-1+(*xoff)] != zero) || ((*y)[j-1+(*yoff)] != zero) {
					temp1 = (*alpha) * conjc64((*y)[j-1+(*yoff)])
					temp2 = conjc64((*alpha) * (*x)[j-1+(*xoff)])
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[j-1+(*xoff)]*temp1+(*y)[j-1+(*yoff)]*temp2), 0.0)
					for i = j + 1; i <= (*n); i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp1 + (*y)[i-1+(*yoff)]*temp2
					}
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				if ((*x)[jx-1+(*xoff)] != zero) || ((*y)[jy-1+(*yoff)] != zero) {
					temp1 = (*alpha) * conjc64((*y)[jy-1+(*yoff)])
					temp2 = conjc64((*alpha) * (*x)[jx-1+(*xoff)])
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)])+real((*x)[jx-1+(*xoff)]*temp1+(*y)[jy-1+(*yoff)]*temp2), 0.0)
					ix = jx
					iy = jy
					for i = j + 1; i <= (*n); i++ {
						ix = ix + (*incx)
						iy = iy + (*incy)
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp1 + (*y)[iy-1+(*yoff)]*temp2
					}
				} else {
					(*a)[j-1+(j-1)*(*lda)+(*aoff)] = complex(real((*a)[j-1+(j-1)*(*lda)+(*aoff)]), 0.0)
				}
				jx = jx + (*incx)
				jy = jy + (*incy)
			}
		}
	}
}
