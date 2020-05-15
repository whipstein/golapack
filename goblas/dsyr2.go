package goblas

// \brief \b Dsyr2
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Dsyr2(UPLO,N,ALPHA,X,incx,Y,incy,A,LDA)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION ALPHA
//       INTEGER incx,incy,LDA,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION A(LDA,//),X(//),Y(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dsyr2  performs the symmetric rank 2 operation
//
//    A := alpha//x//y////T + alpha//y//x////T + A,
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
//          X is DOUBLE PRECISION array, dimension at least
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
// \param[in] Y
// \verbatim
//          Y is DOUBLE PRECISION array, dimension at least
//           ( 1 + ( n - 1)//abs( incy)).
//           Before entry, the incremented array Y must contain the n
//           element vector y.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//           On entry, incy specifies the increment for the elements of
//           Y. incy must not be zero.
// \endverbatim
//
// \param[in,out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension ( LDA, N)
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
// \ingroup double_blas_level2
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
func Dsyr2(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[][]float64, lda *int) {
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
	kx := new(int)
	ky := new(int)
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
	if !Lsame((*uplo), "u") && !Lsame((*uplo), "l") {
		(*info) = 1
	} else if (*n) < 0 {
		(*info) = 2
	} else if (*incx) == 0 {
		(*info) = 5
	} else if (*incy) == 0 {
		(*info) = 7
	} else if (*lda) < max(func() *int { y := 1; return &y }(), n) {
		(*info) = 9
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte { y := []byte("dsyr2 "); return &y }(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*n) == 0) || ((*alpha) == (*zero)) {
		return
	}
	//*
	//*     Set up the start points in X and Y if the increments are not both
	//*     unity.
	//*
	if ((*incx) != 1) || ((*incy) != 1) {
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
		(*jx) = (*kx)
		(*jy) = (*ky)
	}
	//*
	//*     Start the operations. In this version the elements of A are
	//*     accessed sequentially with one pass through the triangular part
	//*     of A.
	//*
	if Lsame(uplo, func() *byte { y := byte('u'); return &y }()) {
		//*
		//*        Form  A  when A is stored in the upper triangle.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*j)-1] != (*zero)) || ((*y)[(*j)-1] != (*zero)) {
					(*temp1) = (*alpha) * (*y)[(*j)-1]
					(*temp2) = (*alpha) * (*x)[(*j)-1]
					for (*i) = 1; (*i) <= (*j); (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*i)-1]*(*temp1) + (*y)[(*i)-1]*(*temp2)
						//Label10:
					}
				}
				//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*jx)-1] != (*zero)) || ((*y)[(*jy)-1] != (*zero)) {
					(*temp1) = (*alpha) * (*y)[(*jy)-1]
					(*temp2) = (*alpha) * (*x)[(*jx)-1]
					(*ix) = (*kx)
					(*iy) = (*ky)
					for (*i) = 1; (*i) <= (*j); (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*ix)-1]*(*temp1) + (*y)[(*iy)-1]*(*temp2)
						(*ix) = (*ix) + (*incx)
						(*iy) = (*iy) + (*incy)
						//Label30:
					}
				}
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				//Label40:
			}
		}
	} else {
		//*
		//*        Form  A  when A is stored in the lower triangle.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*j)-1] != (*zero)) || ((*y)[(*j)-1] != (*zero)) {
					(*temp1) = (*alpha) * (*y)[(*j)-1]
					(*temp2) = (*alpha) * (*x)[(*j)-1]
					for (*i) = (*j); (*i) <= (*n); (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*i)-1]*(*temp1) + (*y)[(*i)-1]*(*temp2)
						//Label50:
					}
				}
				//Label60:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*jx)-1] != (*zero)) || ((*y)[(*jy)-1] != (*zero)) {
					(*temp1) = (*alpha) * (*y)[(*jy)-1]
					(*temp2) = (*alpha) * (*x)[(*jx)-1]
					(*ix) = (*jx)
					(*iy) = (*jy)
					for (*i) = (*j); (*i) <= (*n); (*i)++ {
						(*a)[(*i)-1][(*j)-1] = (*a)[(*i)-1][(*j)-1] + (*x)[(*ix)-1]*(*temp1) + (*y)[(*iy)-1]*(*temp2)
						(*ix) = (*ix) + (*incx)
						(*iy) = (*iy) + (*incy)
						//Label70:
					}
				}
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				//Label80:
			}
		}
	}
	//*
	return
	//*
	//*     End of Dsyr2 .
	//*
}
