package goblas

// \brief \b Sgemv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Sgemv(TRANS,M,N,ALPHA,A,LDA,X,incx,BETA,Y,incy)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER incx,incy,LDA,M,N
//       CHARACTER TRANS
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,//),X(//),Y(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Sgemv  performs one of the matrix-vector operations
//
//    y := alpha//A//x + beta//y,   or   y := alpha//A////T//x + beta//y,
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
//          TRANS is CHARACTER//1
//           On entry, TRANS specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   y := alpha//A//x + beta//y.
//
//              TRANS = 'T' or 't'   y := alpha//A////T//x + beta//y.
//
//              TRANS = 'C' or 'c'   y := alpha//A////T//x + beta//y.
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
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, N)
//           Before entry, the leading m by n part of the array A must
//           contain the matrix of coefficients.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           max( 1, m).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)) when TRANS = 'N' or 'n'
//           and at least
//           ( 1 + ( m - 1)//abs( incx)) otherwise.
//           Before entry, the incremented array X must contain the
//           vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           X. incx must not be zero.
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
//           ( 1 + ( m - 1)//abs( incy)) when TRANS = 'N' or 'n'
//           and at least
//           ( 1 + ( n - 1)//abs( incy)) otherwise.
//           Before entry with BETA non-zero, the incremented array Y
//           must contain the vector y. On exit, Y is overwritten by the
//           updated vector y.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//           On entry, incy specifies the increment for the elements of
//           Y. incy must not be zero.
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
//  The vector and matrix arguments are not referenced when n = 0, or M = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Sgemv(trans *byte, m *int, n *int, alpha *float64, a *[][]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	one := new(float64)
	zero := new(float64)
	temp := new(float64)
	i := new(int)
	info := new(int)
	ix := new(int)
	iy := new(int)
	j := new(int)
	jx := new(int)
	jy := new(int)
	kx := new(int)
	ky := new(int)
	lenx := new(int)
	leny := new(int)
	//*
	//*  -- Reference BLAS level2 routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     December 2016
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. Parameters ..
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
	//*     ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. External Functions ..
	//*     ..
	//*     .. External Subroutines ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	//*
	//*     Test the input parameters.
	//*
	(*info) = 0
	if !Lsame((*trans), "n") && !Lsame((*trans), "t") && !Lsame((*trans), "c") {
		(*info) = 1
	} else if (*m) < 0 {
		(*info) = 2
	} else if (*n) < 0 {
		(*info) = 3
	} else if (*lda) < max(func() *int { y := 1; return &y }(), m) {
		(*info) = 6
	} else if (*incx) == 0 {
		(*info) = 8
	} else if (*incy) == 0 {
		(*info) = 11
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte { y := []byte("sgemv "); return &y }(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*m) == 0) || ((*n) == 0) || (((*alpha) == (*zero)) && ((*beta) == (*one))) {
		return
	}
	//*
	//*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
	//*     up the start points in  X  and  Y.
	//*
	if Lsame(trans, func() *byte { y := byte('n'); return &y }()) {
		(*lenx) = (*n)
		(*leny) = (*m)
	} else {
		(*lenx) = (*m)
		(*leny) = (*n)
	}
	if (*incx) > 0 {
		(*kx) = 1
	} else {
		(*kx) = 1 - ((*lenx)-1)*(*incx)
	}
	if (*incy) > 0 {
		(*ky) = 1
	} else {
		(*ky) = 1 - ((*leny)-1)*(*incy)
	}
	//*
	//*     Start the operations. In this version the elements of A are
	//*     accessed sequentially with one pass through A.
	//*
	//*     First form  y := beta*y.
	//*
	if (*beta) != (*one) {
		if (*incy) == 1 {
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*leny); (*i)++ {
					(*y)[(*i)-1] = (*zero)
					//Label10:
				}
			} else {
				for (*i) = 1; (*i) <= (*leny); (*i)++ {
					(*y)[(*i)-1] = (*beta) * (*y)[(*i)-1]
					//Label20:
				}
			}
		} else {
			(*iy) = (*ky)
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*leny); (*i)++ {
					(*y)[(*iy)-1] = (*zero)
					(*iy) = (*iy) + (*incy)
					//Label30:
				}
			} else {
				for (*i) = 1; (*i) <= (*leny); (*i)++ {
					(*y)[(*iy)-1] = (*beta) * (*y)[(*iy)-1]
					(*iy) = (*iy) + (*incy)
					//Label40:
				}
			}
		}
	}
	if (*alpha) == (*zero) {
		return
	}
	if Lsame(trans, func() *byte { y := byte('n'); return &y }()) {
		//*
		//*        Form  y := alpha*A*x + y.
		//*
		(*jx) = (*kx)
		if (*incy) == 1 {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp) = (*alpha) * (*x)[(*jx)-1]
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*y)[(*i)-1] = (*y)[(*i)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
					//Label50:
				}
				(*jx) = (*jx) + (*incx)
				//Label60:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp) = (*alpha) * (*x)[(*jx)-1]
				(*iy) = (*ky)
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*y)[(*iy)-1] = (*y)[(*iy)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
					(*iy) = (*iy) + (*incy)
					//Label70:
				}
				(*jx) = (*jx) + (*incx)
				//Label80:
			}
		}
	} else {
		//*
		//*        Form  y := alpha*A**T*x + y.
		//*
		(*jy) = (*ky)
		if (*incx) == 1 {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp) = (*zero)
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*i)-1]
					//Label90:
				}
				(*y)[(*jy)-1] = (*y)[(*jy)-1] + (*alpha)*(*temp)
				(*jy) = (*jy) + (*incy)
				//Label100:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp) = (*zero)
				(*ix) = (*kx)
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*ix)-1]
					(*ix) = (*ix) + (*incx)
					//Label110:
				}
				(*y)[(*jy)-1] = (*y)[(*jy)-1] + (*alpha)*(*temp)
				(*jy) = (*jy) + (*incy)
				//Label120:
			}
		}
	}
	//*
	return
	//*
	//*     End of Sgemv .
	//*
}
