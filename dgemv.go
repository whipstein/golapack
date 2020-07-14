package golapack

// Dgemv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION ALPHA,BETA
//       INTEGER INCX,INCY,LDA,M,N
//       CHARACTER TRANS
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// DGEMV  performs one of the matrix-vector operations
//
//    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and A is an
// m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER*1
//           On entry, TRANS specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
//
//              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
//
//              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is INTEGER
//           On entry, M specifies the number of rows of the matrix A.
//           M must be at least zero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the number of columns of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is DOUBLE PRECISION.
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension ( LDA, N )
//           Before entry, the leading m by n part of the array A must
//           contain the matrix of coefficients.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           max( 1, m ).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension at least
//           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
//           and at least
//           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
//           Before entry, the incremented array X must contain the
//           vector x.
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
//          BETA is DOUBLE PRECISION.
//           On entry, BETA specifies the scalar beta. When BETA is
//           supplied as zero then Y need not be set on input.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is DOUBLE PRECISION array, dimension at least
//           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
//           and at least
//           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
//           Before entry with BETA non-zero, the incremented array Y
//           must contain the vector y. On exit, Y is overwritten by the
//           updated vector y.
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
// \ingroup double_blas_level2
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
func Dgemv(trans *byte, m *int, n *int, alpha *float64, a *[]float64, aoff, lda *int, x *[]float64, xoff, incx *int, beta *float64, y *[]float64, yoff, incy *int) {
	var temp float64
	var i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny int
	var one float64 = 1.0
	var zero float64 = 0.0

	//
	//     Test the input parameters.
	//
	info = 0
	if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 1
	} else if (*m) < 0 {
		info = 2
	} else if (*n) < 0 {
		info = 3
	} else if (*lda) < maxint(1, *m) {
		info = 6
	} else if (*incx) == 0 {
		info = 8
	} else if (*incy) == 0 {
		info = 11
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("DGEMV "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*m) == 0) || ((*n) == 0) || (((*alpha) == zero) && ((*beta) == one)) {
		return
	}
	//
	//     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
	//     up the start points in  X  and  Y.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		lenx = (*n)
		leny = (*m)
	} else {
		lenx = (*m)
		leny = (*n)
	}
	if (*incx) > 0 {
		kx = 1
	} else {
		kx = 1 - (lenx-1)*(*incx)
	}
	if (*incy) > 0 {
		ky = 1
	} else {
		ky = 1 - (leny-1)*(*incy)
	}
	//
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through A.
	//
	//     First form  y := beta*y.
	//
	if (*beta) != one {
		if (*incy) == 1 {
			if (*beta) == zero {
				for i = 1; i <= leny; i++ {
					(*y)[i-1+(*yoff)] = zero
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[i-1+(*yoff)] = (*beta) * (*y)[i-1+(*yoff)]
				}
			}
		} else {
			iy = ky
			if (*beta) == zero {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1+(*yoff)] = zero
					iy = iy + (*incy)
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1+(*yoff)] = (*beta) * (*y)[iy-1+(*yoff)]
					iy = iy + (*incy)
				}
			}
		}
	}
	if (*alpha) == zero {
		return
	}
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  y := alpha*A*x + y.
		//
		jx = kx
		if (*incy) == 1 {
			for j = 1; j <= (*n); j++ {
				temp = (*alpha) * (*x)[jx-1+(*xoff)]
				for i = 1; i <= (*m); i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
				}
				jx = jx + (*incx)
			}
		} else {
			for j = 1; j <= (*n); j++ {
				temp = (*alpha) * (*x)[jx-1+(*xoff)]
				iy = ky
				for i = 1; i <= (*m); i++ {
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp*(*a)[i-1+(j-1)*(*lda)+(*aoff)]
					iy = iy + (*incy)
				}
				jx = jx + (*incx)
			}
		}
	} else {
		//
		//        Form  y := alpha*A**T*x + y.
		//
		jy = ky
		if (*incx) == 1 {
			for j = 1; j <= (*n); j++ {
				temp = zero
				for i = 1; i <= (*m); i++ {
					temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp
				jy = jy + (*incy)
			}
		} else {
			for j = 1; j <= (*n); j++ {
				temp = zero
				ix = kx
				for i = 1; i <= (*m); i++ {
					temp = temp + (*a)[i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
					ix = ix + (*incx)
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp
				jy = jy + (*incy)
			}
		}
	}
}
