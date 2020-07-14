package golapack

// Sgbmv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER INCX,INCY,KL,KU,LDA,M,N
//       CHARACTER TRANS
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
// SGBMV  performs one of the matrix-vector operations
//
//    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and A is an
// m by n band matrix, with kl sub-diagonals and ku super-diagonals.
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
// \param[in] KL
// \verbatim
//          KL is INTEGER
//           On entry, KL specifies the number of sub-diagonals of the
//           matrix A. KL must satisfy  0 .le. KL.
// \endverbatim
//
// \param[in] KU
// \verbatim
//          KU is INTEGER
//           On entry, KU specifies the number of super-diagonals of the
//           matrix A. KU must satisfy  0 .le. KU.
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
//           Before entry, the leading ( kl + ku + 1 ) by n part of the
//           array A must contain the matrix of coefficients, supplied
//           column by column, with the leading diagonal of the matrix in
//           row ( ku + 1 ) of the array, the first super-diagonal
//           starting at position 2 in row ku, the first sub-diagonal
//           starting at position 1 in row ( ku + 2 ), and so on.
//           Elements in the array A that do not correspond to elements
//           in the band matrix (such as the top left ku by ku triangle)
//           are not referenced.
//           The following program segment will transfer a band matrix
//           from conventional full matrix storage to band storage:
//
//                 DO 20, J = 1, N
//                    K = KU + 1 - J
//                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
//                       A( K + I, J ) = matrix( I, J )
//              10    CONTINUE
//              20 CONTINUE
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           ( kl + ku + 1 ).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
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
//          BETA is REAL
//           On entry, BETA specifies the scalar beta. When BETA is
//           supplied as zero then Y need not be set on input.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
//           and at least
//           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
//           Before entry, the incremented array Y must contain the
//           vector y. On exit, Y is overwritten by the updated vector y.
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
func Sgbmv(trans *byte, m *int, n *int, kl *int, ku *int, alpha *float32, a *[]float32, aoff, lda *int, x *[]float32, xoff, incx *int, beta *float32, y *[]float32, yoff, incy *int) {
	var temp float32
	var i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx, leny int
	var one float32 = 1.0
	var zero float32 = 0.0

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
	} else if (*kl) < 0 {
		info = 4
	} else if (*ku) < 0 {
		info = 5
	} else if (*lda) < ((*kl) + (*ku) + 1) {
		info = 8
	} else if (*incx) == 0 {
		info = 10
	} else if (*incy) == 0 {
		info = 13
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("SGBMV "); return &y }(), &info)
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
	//     accessed sequentially with one pass through the band part of A.
	//
	//     First form  y := beta*y.
	//
	if (*beta) != one {
		if (*incy) == 1 {
			if (*beta) == zero {
				for i = 1; i <= leny; i++ {
					(*y)[i-1+(*yoff)] = zero
					//Label10:
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[i-1+(*yoff)] = (*beta) * (*y)[i-1+(*yoff)]
					//Label20:
				}
			}
		} else {
			iy = ky
			if (*beta) == zero {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1+(*yoff)] = zero
					iy = iy + (*incy)
					//Label30:
				}
			} else {
				for i = 1; i <= leny; i++ {
					(*y)[iy-1+(*yoff)] = (*beta) * (*y)[iy-1+(*yoff)]
					iy = iy + (*incy)
					//Label40:
				}
			}
		}
	}
	if (*alpha) == zero {
		return
	}
	kup1 = (*ku) + 1
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  y := alpha*A*x + y.
		//
		jx = kx
		if (*incy) == 1 {
			for j = 1; j <= (*n); j++ {
				temp = (*alpha) * (*x)[jx-1+(*xoff)]
				k = kup1 - j
				for i = maxint(1, j-(*ku)); i <= minint(*m, j+(*kl)); i++ {
					(*y)[i-1+(*yoff)] = (*y)[i-1+(*yoff)] + temp*(*a)[k+i-1+(j-1)*(*lda)+(*aoff)]
					//Label50:
				}
				jx = jx + (*incx)
				//Label60:
			}
		} else {
			for j = 1; j <= (*n); j++ {
				temp = (*alpha) * (*x)[jx-1+(*xoff)]
				iy = ky
				k = kup1 - j
				for i = maxint(1, j-(*ku)); i <= minint(*m, j+(*kl)); i++ {
					(*y)[iy-1+(*yoff)] = (*y)[iy-1+(*yoff)] + temp*(*a)[k+i-1+(j-1)*(*lda)+(*aoff)]
					iy = iy + (*incy)
					//Label70:
				}
				jx = jx + (*incx)
				if j > (*ku) {
					ky = ky + (*incy)
				}
				//Label80:
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
				k = kup1 - j
				for i = maxint(1, j-(*ku)); i <= minint(*m, j+(*kl)); i++ {
					temp = temp + (*a)[k+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[i-1+(*xoff)]
					//Label90:
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp
				jy = jy + (*incy)
				//Label100:
			}
		} else {
			for j = 1; j <= (*n); j++ {
				temp = zero
				ix = kx
				k = kup1 - j
				for i = maxint(1, j-(*ku)); i <= minint(*m, j+(*kl)); i++ {
					temp = temp + (*a)[k+i-1+(j-1)*(*lda)+(*aoff)]*(*x)[ix-1+(*xoff)]
					ix = ix + (*incx)
					//Label110:
				}
				(*y)[jy-1+(*yoff)] = (*y)[jy-1+(*yoff)] + (*alpha)*temp
				jy = jy + (*incy)
				if j > (*ku) {
					kx = kx + (*incx)
				}
				//Label120:
			}
		}
	}
	//
	//
	//     End of SGBMV .
	//
}
