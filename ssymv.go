package golapack

// Ssymv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
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
// SSYMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n symmetric matrix.
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
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, N )
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular part of the symmetric matrix and the strictly
//           lower triangular part of A is not referenced.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular part of the symmetric matrix and the strictly
//           upper triangular part of A is not referenced.
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
// \param[in] BETA
// \verbatim
//          BETA is REAL
//           On entry, BETA specifies the scalar beta. When BETA is
//           supplied as zero then Y need not be set on input.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCY ) ).
//           Before entry, the incremented array Y must contain the n
//           element vector y. On exit, Y is overwritten by the updated
//           vector y.
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//           On entry, INCY specifies the increment for the elements of
//           Y. INCY must not be zero.
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
func Ssymv(uplo *byte, n *int, alpha *float32, a *[]float32, aoff, lda *int, x *[]float32, xoff, incx *int, beta *float32, y *[]float32, yoff, incy *int) {
	var temp1, temp2 float32
	var i, info, ix, iy, j, jx, jy, kx, ky int
	var one float32 = 1.0
	var zero float32 = 0.0

	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 1
	} else if (*n) < 0 {
		info = 2
	} else if (*lda) < maxint(1, *n) {
		info = 5
	} else if (*incx) == 0 {
		info = 7
	} else if (*incy) == 0 {
		info = 10
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("SSYMV "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || (((*alpha) == zero) && ((*beta) == one)) {
		return
	}
	//
	//     Set up the start points in  X  and  Y.
	//
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
	//
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through the triangular part
	//     of A.
	//
	//     First form  y := beta*y.
	//
	if (*beta) != one {
		if (*incy) == 1 {
			if (*beta) == zero {
				for i = 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = zero
				}
			} else {
				for i = 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = (*beta) * (*y)[i-1+(*yoff)]
				}
			}
		} else {
			iy = ky
			if (*beta) == zero {
				for i = 1; i <= (*n); i++ {
					(*y)[iy-1+(*yoff)] = zero
					iy = iy + (*incy)
				}
			} else {
				for i = 1; i <= (*n); i++ {
					(*y)[iy-1+(*yoff)] = (*beta) * (*y)[iy-1+(*yoff)]
					iy = iy + (*incy)
				}
			}
		}
	}
	if (*alpha) == zero {
		return
	}
	if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
		//
		//        Form  y  when A is stored in upper triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[j-1+(*xoff)]
				temp2 = zero
				for i = 1; i <= j-1; i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp1*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
					temp2 = temp2 + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
				}
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + temp1*(*a)[j-1+(j-1)*(*lda)+(*aoff)] + (*alpha)*temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[jx-1+(*xoff)]
				temp2 = zero
				ix = kx
				iy = ky
				for i = 1; i <= j-1; i++ {
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp1*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
					temp2 = temp2 + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
					ix = ix + (*incx)
					iy = iy + (*incy)
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + temp1*(*a)[j-1+(j-1)*(*lda)+(*aoff)] + (*alpha)*temp2
				jx = jx + (*incx)
				jy = jy + (*incy)
			}
		}
	} else {
		//
		//        Form  y  when A is stored in lower triangle.
		//
		if ((*incx) == 1) && ((*incy) == 1) {
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[j-1+(*xoff)]
				temp2 = zero
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + temp1*(*a)[j-1+(j-1)*(*lda)+(*aoff)]
				for i = j + 1; i <= (*n); i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp1*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
					temp2 = temp2 + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
				}
				(*y)[j-1+(*yoff)] = (*y)[j-1+(*yoff)] + (*alpha)*temp2
			}
		} else {
			jx = kx
			jy = ky
			for j = 1; j <= (*n); j++ {
				temp1 = (*alpha) * (*x)[jx-1+(*xoff)]
				temp2 = zero
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + temp1*(*a)[j-1+(j-1)*(*lda)+(*aoff)]
				ix = jx
				iy = jy
				for i = j + 1; i <= (*n); i++ {
					ix = ix + (*incx)
					iy = iy + (*incy)
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp1*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
					temp2 = temp2 + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp2
				jx = jx + (*incx)
				jy = jy + (*incy)
			}
		}
	}
}
