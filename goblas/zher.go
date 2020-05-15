package goblas

import 

// \brief \b Zher
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zher(UPLO,N,ALPHA,X,incx,A,LDA)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION ALPHA
//       INTEGER incx,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 A(LDA,//),X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zher   performs the hermitian rank 1 operation
//
//    A := alpha//x//x////H + A,
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
//          UPLO is CHARACTER//1
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
//          X is COMPLEX//16 array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)).
//           Before entry, the incremented array X must contain the n
//           element vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           X. incx must not be zero.
// \endverbatim
//
// \param[in,out] A
// \verbatim
//          A is COMPLEX//16 array, dimension ( LDA, N)
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
//           max( 1, n).
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
func Zher(uplo *byte, n *int, alpha *float64, x *[]complex128, incx *int, a *[][]complex128, lda *int) {
	zero := new(complex128)
	temp := new(complex128)
	i := new(int)
	info := new(int)
	ix := new(int)
	j := new(int)
	jx := new(int)
	kx := new(int)
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
	(*zero) = (0.0e+0 + (0.0e+0)*1i)
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
	if !Lsame((*uplo), "u") && !Lsame((*uplo), "l") {
		(*info) = 1
	} else if (*n) < 0 {
		(*info) = 2
	} else if (*incx) == 0 {
		(*info) = 5
	} else if (*lda) < max(func() *int {y := 1; return &y}(), n) {
		(*info) = 7
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("zher  "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*n) == 0) || ((*alpha) == dble((*zero))) {
		return
	}
	//*
	//*     Set the start point in X if the increment is not unity.
	//*
	if (*incx) <= 0 {
		(*kx) = 1 - ((*n)-1)*(*incx)
	} else if (*incx) != 1 {
		(*kx) = 1
	}
	//*
	//*     Start the operations. In this version the elements of A are
	//*     accessed sequentially with one pass through the triangular part
	//*     of A.
	//*
	if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
		//*
		//*        Form  A  when A is stored in upper triangle.
		//*
		if (*incx) == 1 {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*x)[(*j)-1] != (*zero) {
					(*temp) = (*alpha) * DCONJG(((*x)[(*j)-1]))
					for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*i)-1]*(*temp)
						//Label10:
					}
					(*a)[(*j)-1][(*j)-1] = DBLE(((*a)[(*j)-1][(*j)-1])) + DBLE((*x)[(*j)-1]*(*temp))
				} else {
					(*a)[(*j)-1][(*j)-1] = (DBLE(((*a)[(*j)-1][(*j)-1])))
				}
				//Label20:
			}
		} else {
			(*jx) = (*kx)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*x)[(*jx)-1] != (*zero) {
					(*temp) = (*alpha) * DCONJG(((*x)[(*jx)-1]))
					(*ix) = (*kx)
					for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*ix)-1]*(*temp)
						(*ix) = (*ix) + (*incx)
						//Label30:
					}
					(*a)[(*j)-1][(*j)-1] = DBLE(((*a)[(*j)-1][(*j)-1])) + DBLE((*x)[(*jx)-1]*(*temp))
				} else {
					(*a)[(*j)-1][(*j)-1] = (DBLE(((*a)[(*j)-1][(*j)-1])))
				}
				(*jx) = (*jx) + (*incx)
				//Label40:
			}
		}
	} else {
		//*
		//*        Form  A  when A is stored in lower triangle.
		//*
		if (*incx) == 1 {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*x)[(*j)-1] != (*zero) {
					(*temp) = (*alpha) * DCONJG(((*x)[(*j)-1]))
					(*a)[(*j)-1][(*j)-1] = DBLE(((*a)[(*j)-1][(*j)-1])) + DBLE((*temp)*(*x)[(*j)-1])
					for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*i)-1]*(*temp)
						//Label50:
					}
				} else {
					(*a)[(*j)-1][(*j)-1] = (DBLE(((*a)[(*j)-1][(*j)-1])))
				}
				//Label60:
			}
		} else {
			(*jx) = (*kx)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*x)[(*jx)-1] != (*zero) {
					(*temp) = (*alpha) * DCONJG(((*x)[(*jx)-1]))
					(*a)[(*j)-1][(*j)-1] = DBLE(((*a)[(*j)-1][(*j)-1])) + DBLE((*temp)*(*x)[(*jx)-1])
					(*ix) = (*jx)
					for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
						(*ix) = (*ix) + (*incx)
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*ix)-1]*(*temp)
						//Label70:
					}
				} else {
					(*a)[(*j)-1][(*j)-1] = (DBLE(((*a)[(*j)-1][(*j)-1])))
				}
				(*jx) = (*jx) + (*incx)
				//Label80:
			}
		}
	}
	//*
	return
	//*
	//*     End of Zher  .
	//*
}
