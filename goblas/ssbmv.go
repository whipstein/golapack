package goblas

import 
// \brief \b Ssbmv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ssbmv(UPLO,N,K,ALPHA,A,LDA,X,incx,BETA,Y,incy)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER incx,incy,K,LDA,N
//       CHARACTER UPLO
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
// Ssbmv  performs the matrix-vector  operation
//
//    y := alpha//A//x + beta//y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n symmetric band matrix, with k super-diagonals.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the band matrix A is being supplied as
//           follows:
//
//              UPLO = 'U' or 'u'   The upper triangular part of A is
//                                  being supplied.
//
//              UPLO = 'L' or 'l'   The lower triangular part of A is
//                                  being supplied.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry, K specifies the number of super-diagonals of the
//           matrix A. K must satisfy  0 .le. K.
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
//           Before entry with UPLO = 'U' or 'u', the leading ( k + 1)
//           by n part of the array A must contain the upper triangular
//           band part of the symmetric matrix, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array A is not referenced.
//           The following program segment will transfer the upper
//           triangular part of a symmetric band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, J = 1, N
//                    M = K + 1 - J
//                    DO 10, I = MAX( 1, J - K), J
//                       A( M + I, J) = matrix( I, J)
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with UPLO = 'L' or 'l', the leading ( k + 1)
//           by n part of the array A must contain the lower triangular
//           band part of the symmetric matrix, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array A is not referenced.
//           The following program segment will transfer the lower
//           triangular part of a symmetric band matrix from conventional
//           full matrix storage to band storage:
//
//                 DO 20, J = 1, N
//                    M = 1 - J
//                    DO 10, I = J, MIN( N, J + K)
//                       A( M + I, J) = matrix( I, J)
//              10    CONTINUE
//              20 CONTINUE
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           ( k + 1).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)).
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
//           On entry, BETA specifies the scalar beta.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is REAL array, dimension at least
//           ( 1 + ( n - 1)//abs( incy)).
//           Before entry, the incremented array Y must contain the
//           vector y. On exit, Y is overwritten by the updated vector y.
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
func Ssbmv(uplo *byte, n *int, k *int, alpha *float64, a *[][]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	one := new(float64)
	zero := new(float64)
	temp1 := new(float64)
	temp2 := new(float64)
	i := new(int)
	info := new(int)
	ix := new(int)
	iy := new(int)
	j := new(int)
	jx := new(int)
	jy := new(int)
	kplus1 := new(int)
	kx := new(int)
	ky := new(int)
	l := new(int)
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
	if !Lsame((*uplo), "u") && . !Lsame((*uplo), "l") {
		(*info) = 1
	} else if (*n) < 0 {
		(*info) = 2
	} else if (*k) < 0 {
		(*info) = 3
	} else if (*lda) < ((*k) + 1) {
		(*info) = 6
	} else if (*incx) == 0 {
		(*info) = 8
	} else if (*incy) == 0 {
		(*info) = 11
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte{y := []byte("ssbmv "); return &y }(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*n) == 0) || ( ((*alpha) == (*zero)) && ((*beta) == (*one))) {
		return
	}
	//*
	//*     Set up the start points in  X  and  Y.
	//*
	if (*incx) > 0 {
		(*kx) = 1
	} else {
		(*kx) = 1 - ((*n)-1)*(*incx)
	}
	if (*incy) > 0 {
		(*ky) = 1
	} else {
		(*ky) = 1 - ((*n)-1)*(*incy)
	}
	//*
	//*     Start the operations. In this version the elements of the array A
	//*     are accessed sequentially with one pass through A.
	//*
	//*     First form  y := beta*y.
	//*
	if (*beta) != (*one) {
		if (*incy) == 1 {
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*i)-(1)] = (*zero)
				//Label10:
				}
			} else {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*i)-(1)] = (*beta) * (*y)[(*i)-(1)]
				//Label20:
				}
			}
		} else {
			(*iy) = (*ky)
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*iy)-(1)] = (*zero)
					(*iy) = (*iy) + (*incy)
				//Label30:
				}
			} else {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*iy)-(1)] = (*beta) * (*y)[(*iy)-(1)]
					(*iy) = (*iy) + (*incy)
				//Label40:
				}
			}
		}
	}
	if (*alpha) == (*zero) {
		return
	}
	if Lsame(uplo, func() *byte{y := byte('u'); return &y }()) {
		//*
		//*        Form  y  when upper triangle of A is stored.
		//*
		(*kplus1) = (*k) + 1
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*j)-(1)]
				(*temp2) = (*zero)
				(*l) = (*kplus1) - (*j)
				for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
					(*y)[(*i)-(1)] = (*y)[(*i)-(1)] + (*temp1)*(*a)[(*l)+(*i)-(1)][(*j)-(1)]
					(*temp2) = (*temp2) + (*a)[(*l)+(*i)-(1)][(*j)-(1)]*(*x)[(*i)-(1)]
				//Label50:
				}
				(*y)[(*j)-(1)] = (*y)[(*j)-(1)] + (*temp1)*(*a)[(*kplus1)-(1)][(*j)-(1)] + (*alpha)*(*temp2)
			//Label60:
			}
		} else {
			(*jx) = (*kx)
			(*jy) = (*ky)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*jx)-(1)]
				(*temp2) = (*zero)
				(*ix) = (*kx)
				(*iy) = (*ky)
				(*l) = (*kplus1) - (*j)
				for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
					(*y)[(*iy)-(1)] = (*y)[(*iy)-(1)] + (*temp1)*(*a)[(*l)+(*i)-(1)][(*j)-(1)]
					(*temp2) = (*temp2) + (*a)[(*l)+(*i)-(1)][(*j)-(1)]*(*x)[(*ix)-(1)]
					(*ix) = (*ix) + (*incx)
					(*iy) = (*iy) + (*incy)
				//Label70:
				}
				(*y)[(*jy)-(1)] = (*y)[(*jy)-(1)] + (*temp1)*(*a)[(*kplus1)-(1)][(*j)-(1)] + (*alpha)*(*temp2)
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				if (*j) > (*k) {
					(*kx) = (*kx) + (*incx)
					(*ky) = (*ky) + (*incy)
				}
			//Label80:
			}
		}
	} else {
		//*
		//*        Form  y  when lower triangle of A is stored.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*j)-(1)]
				(*temp2) = (*zero)
				(*y)[(*j)-(1)] = (*y)[(*j)-(1)] + (*temp1)*(*a)[0][(*j)-(1)]
				(*l) = 1 - (*j)
				for (*i) = (*j) + 1; (*i) <= (MIN((*n), (*j)+(*k))); (*i)++ {
					(*y)[(*i)-(1)] = (*y)[(*i)-(1)] + (*temp1)*(*a)[(*l)+(*i)-(1)][(*j)-(1)]
					(*temp2) = (*temp2) + (*a)[(*l)+(*i)-(1)][(*j)-(1)]*(*x)[(*i)-(1)]
				//Label90:
				}
				(*y)[(*j)-(1)] = (*y)[(*j)-(1)] + (*alpha)*(*temp2)
			//Label100:
			}
		} else {
			(*jx) = (*kx)
			(*jy) = (*ky)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*jx)-(1)]
				(*temp2) = (*zero)
				(*y)[(*jy)-(1)] = (*y)[(*jy)-(1)] + (*temp1)*(*a)[0][(*j)-(1)]
				(*l) = 1 - (*j)
				(*ix) = (*jx)
				(*iy) = (*jy)
				for (*i) = (*j) + 1; (*i) <= (MIN((*n), (*j)+(*k))); (*i)++ {
					(*ix) = (*ix) + (*incx)
					(*iy) = (*iy) + (*incy)
					(*y)[(*iy)-(1)] = (*y)[(*iy)-(1)] + (*temp1)*(*a)[(*l)+(*i)-(1)][(*j)-(1)]
					(*temp2) = (*temp2) + (*a)[(*l)+(*i)-(1)][(*j)-(1)]*(*x)[(*ix)-(1)]
				//Label110:
				}
				(*y)[(*jy)-(1)] = (*y)[(*jy)-(1)] + (*alpha)*(*temp2)
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
			//Label120:
			}
		}
	}
	//*
	return
	//*
	//*     End of Ssbmv .
	//*
}