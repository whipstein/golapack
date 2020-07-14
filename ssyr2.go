package golapack

// Ssyr2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
//       INTEGER INCX,INCY,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,*),X(*),Y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSYR2  performs the symmetric rank 2 operation
//
//    A := alpha*x*y**T + alpha*y*x**T + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an n
// by n symmetric matrix.
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
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
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
//          Y is REAL array, dimension at least
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
//          A is REAL array, dimension ( LDA, N )
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular part of the symmetric matrix and the strictly
//           lower triangular part of A is not referenced. On exit, the
//           upper triangular part of the array A is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular part of the symmetric matrix and the strictly
//           upper triangular part of A is not referenced. On exit, the
//           lower triangular part of the array A is overwritten by the
//           lower triangular part of the updated matrix.
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
func Ssyr2(uplo *byte, n *int, alpha *float32, x *[]float32, xoff, incx *int, y *[]float32, yoff, incy *int, a *[]float32, aoff, lda *int) {
	var temp1, temp2 float32
	var i, info, ix, iy, j, jx, jy, kx, ky int
	var zero float32 = 0.0

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
		Xerbla(func() *[]byte { y := []byte("SSYR2 "); return &y }(), &info)
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
					temp1 = (*alpha) * (*y)[j-1+(*yoff)]
					temp2 = (*alpha) * (*x)[j-1+(*xoff)]
					for i = 1; i <= j; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp1 + (*y)[i-1+(*yoff)]*temp2
					}
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				if ((*x)[jx-1+(*xoff)] != zero) || ((*y)[jy-1+(*yoff)] != zero) {
					temp1 = (*alpha) * (*y)[jy-1+(*yoff)]
					temp2 = (*alpha) * (*x)[jx-1+(*xoff)]
					ix = kx
					iy = ky
					for i = 1; i <= j; i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp1 + (*y)[iy-1+(*yoff)]*temp2
						ix = ix + (*incx)
						iy = iy + (*incy)
					}
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
					temp1 = (*alpha) * (*y)[j-1+(*yoff)]
					temp2 = (*alpha) * (*x)[j-1+(*xoff)]
					for i = j; i <= (*n); i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp1 + (*y)[i-1+(*yoff)]*temp2
					}
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				if ((*x)[jx-1+(*xoff)] != zero) || ((*y)[jy-1+(*yoff)] != zero) {
					temp1 = (*alpha) * (*y)[jy-1+(*yoff)]
					temp2 = (*alpha) * (*x)[jx-1+(*xoff)]
					ix = jx
					iy = jy
					for i = j; i <= (*n); i++ {
						(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp1 + (*y)[iy-1+(*yoff)]*temp2
						ix = ix + (*incx)
						iy = iy + (*incy)
					}
				}
				jx = jx + (*incx)
				jy = jy + (*incy)
			}
		}
	}
}
