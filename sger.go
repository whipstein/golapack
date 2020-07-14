package golapack

// Sger ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
//       INTEGER INCX,INCY,LDA,M,N
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
// SGER   performs the rank 1 operation
//
//    A := alpha*x*y**T + A,
//
// where alpha is a scalar, x is an m element vector, y is an n element
// vector and A is an m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
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
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( m - 1 )*abs( INCX ) ).
//           Before entry, the incremented array X must contain the m
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
//           Before entry, the leading m by n part of the array A must
//           contain the matrix of coefficients. On exit, A is
//           overwritten by the updated matrix.
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
func Sger(m *int, n *int, alpha *float32, x *[]float32, xoff, incx *int, y *[]float32, yoff, incy *int, a *[]float32, aoff, lda *int) {
	var temp float32
	var i, info, ix, j, jy, kx int
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
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//
	//     Test the input parameters.
	//
	info = 0
	if (*m) < 0 {
		info = 1
	} else if (*n) < 0 {
		info = 2
	} else if (*incx) == 0 {
		info = 5
	} else if (*incy) == 0 {
		info = 7
	} else if (*lda) < maxint(1, *m) {
		info = 9
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("SGER  "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*m) == 0) || ((*n) == 0) || ((*alpha) == zero) {
		return
	}
	//
	//     Start the operations. In this version the elements of A are
	//     accessed sequentially with one pass through A.
	//
	if (*incy) > 0 {
		jy = 1
	} else {
		jy = 1 - ((*n)-1)*(*incy)
	}
	if (*incx) == 1 {
		for j = 1; j <= (*n); j++ {
			if (*y)[jy-1+(*yoff)] != zero {
				temp = (*alpha) * (*y)[jy-1+(*yoff)]
				for i = 1; i <= (*m); i++ {
					(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[i-1+(*xoff)]*temp
					//Label10:
				}
			}
			jy = jy + (*incy)
			//Label20:
		}
	} else {
		if (*incx) > 0 {
			kx = 1
		} else {
			kx = 1 - ((*m)-1)*(*incx)
		}
		for j = 1; j <= (*n); j++ {
			if (*y)[jy-1+(*yoff)] != zero {
				temp = (*alpha) * (*y)[jy-1+(*yoff)]
				ix = kx
				for i = 1; i <= (*m); i++ {
					(*a)[i-1+(j-1)*(*lda)+(*aoff)] = (*a)[i-1+(j-1)*(*lda)+(*aoff)] + (*x)[ix-1+(*xoff)]*temp
					ix = ix + (*incx)
					//Label30:
				}
			}
			jy = jy + (*incy)
			//Label40:
		}
	}
	//
	//
	//     End of SGER  .
	//
}
