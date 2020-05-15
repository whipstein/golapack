package goblas

import 
// \brief \b Zhpr2
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zhpr2(UPLO,N,ALPHA,X,incx,Y,incy,AP)
//
//       .. Scalar Arguments ..
//       COMPLEX//16 ALPHA
//       INTEGER incx,incy,N
//       CHARACTER UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 AP(//),X(//),Y(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zhpr2  performs the hermitian rank 2 operation
//
//    A := alpha//x//y////H + conjg( alpha)//y//x////H + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an
// n by n hermitian matrix, supplied in packed form.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
//           On entry, UPLO specifies whether the upper or lower
//           triangular part of the matrix A is supplied in the packed
//           array AP as follows:
//
//              UPLO = 'U' or 'u'   The upper triangular part of A is
//                                  supplied in AP.
//
//              UPLO = 'L' or 'l'   The lower triangular part of A is
//                                  supplied in AP.
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
//          ALPHA is COMPLEX//16
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
// \param[in] Y
// \verbatim
//          Y is COMPLEX//16 array, dimension at least
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
// \param[in,out] AP
// \verbatim
//          AP is COMPLEX//16 array, dimension at least
//           ( ( n//( n + 1))/2).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular part of the hermitian matrix
//           packed sequentially, column by column, so that AP(1)
//           contains a( 1, 1), AP( 2) and AP( 3) contain a( 1, 2)
//           and a( 2, 2) respectively, and so on. On exit, the array
//           AP is overwritten by the upper triangular part of the
//           updated matrix.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular part of the hermitian matrix
//           packed sequentially, column by column, so that AP(1)
//           contains a( 1, 1), AP( 2) and AP( 3) contain a( 2, 1)
//           and a( 3, 1) respectively, and so on. On exit, the array
//           AP is overwritten by the lower triangular part of the
//           updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set, they are assumed to be zero, and on exit they
//           are set to zero.
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
func Zhpr2(uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, ap *[]complex128) {
	zero := new(complex128)
	temp1 := new(complex128)
	temp2 := new(complex128)
	i := new(int)
	info := new(int)
	ix := new(int)
	iy := new(int)
	j := new(int)
	jx := new(int)
	jy := new(int)
	k := new(int)
	kk := new(int)
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
	if !Lsame((*uplo), "u") && . !Lsame((*uplo), "l") {
		(*info) = 1
	} else if (*n) < 0 {
		(*info) = 2
	} else if (*incx) == 0 {
		(*info) = 5
	} else if (*incy) == 0 {
		(*info) = 7
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte{y := []byte("zhpr2 "); return &y }(), info)
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
	//*     Start the operations. In this version the elements of the array AP
	//*     are accessed sequentially with one pass through AP.
	//*
	(*kk) = 1
	if Lsame(uplo, func() *byte{y := byte('u'); return &y }()) {
		//*
		//*        Form  A  when upper triangle is stored in AP.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*j) - (1)] != (*zero)) || ((*y)[(*j) - (1)] != (*zero)) {
					(*temp1) = (*alpha) * DCONJG(((*y)[(*j)-(1)]))
					(*temp2) = (DCONJG((*alpha) * (*x)[(*j)-(1)]))
					(*k) = (*kk)
					for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
						(*ap)[(*k)-(1)] = (*ap)[(*k)-(1)] + (*x)[(*i)-(1)]*(*temp1) + (*y)[(*i)-(1)]*(*temp2)
						(*k) = (*k) + 1
					//Label10:
					}
					(*ap)[(*kk)+(*j)-0] = DBLE(((*ap)[(*kk)+(*j)-0])) + DBLE((*x)[(*j)-(1)]*(*temp1)+(*y)[(*j)-(1)]*(*temp2))
				} else {
					(*ap)[(*kk)+(*j)-0] = (DBLE(((*ap)[(*kk)+(*j)-0])))
				}
				(*kk) = (*kk) + (*j)
			//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*jx) - (1)] != (*zero)) || ((*y)[(*jy) - (1)] != (*zero)) {
					(*temp1) = (*alpha) * DCONJG(((*y)[(*jy)-(1)]))
					(*temp2) = (DCONJG((*alpha) * (*x)[(*jx)-(1)]))
					(*ix) = (*kx)
					(*iy) = (*ky)
					for (*k) = (*kk); (*k) <= (*kk)+(*j)-2; (*k)++ {
						(*ap)[(*k)-(1)] = (*ap)[(*k)-(1)] + (*x)[(*ix)-(1)]*(*temp1) + (*y)[(*iy)-(1)]*(*temp2)
						(*ix) = (*ix) + (*incx)
						(*iy) = (*iy) + (*incy)
					//Label30:
					}
					(*ap)[(*kk)+(*j)-0] = DBLE(((*ap)[(*kk)+(*j)-0])) + DBLE((*x)[(*jx)-(1)]*(*temp1)+(*y)[(*jy)-(1)]*(*temp2))
				} else {
					(*ap)[(*kk)+(*j)-0] = (DBLE(((*ap)[(*kk)+(*j)-0])))
				}
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				(*kk) = (*kk) + (*j)
			//Label40:
			}
		}
	} else {
		//*
		//*        Form  A  when lower triangle is stored in AP.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*j) - (1)] != (*zero)) || ((*y)[(*j) - (1)] != (*zero)) {
					(*temp1) = (*alpha) * DCONJG(((*y)[(*j)-(1)]))
					(*temp2) = (DCONJG((*alpha) * (*x)[(*j)-(1)]))
					(*ap)[(*kk)-(1)] = DBLE(((*ap)[(*kk)-(1)])) + DBLE((*x)[(*j)-(1)]*(*temp1)+(*y)[(*j)-(1)]*(*temp2))
					(*k) = (*kk) + 1
					for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
						(*ap)[(*k)-(1)] = (*ap)[(*k)-(1)] + (*x)[(*i)-(1)]*(*temp1) + (*y)[(*i)-(1)]*(*temp2)
						(*k) = (*k) + 1
					//Label50:
					}
				} else {
					(*ap)[(*kk)-(1)] = (DBLE(((*ap)[(*kk)-(1)])))
				}
				(*kk) = (*kk) + (*n) - (*j) + 1
			//Label60:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if ((*x)[(*jx) - (1)] != (*zero)) || ((*y)[(*jy) - (1)] != (*zero)) {
					(*temp1) = (*alpha) * DCONJG(((*y)[(*jy)-(1)]))
					(*temp2) = (DCONJG((*alpha) * (*x)[(*jx)-(1)]))
					(*ap)[(*kk)-(1)] = DBLE(((*ap)[(*kk)-(1)])) + DBLE((*x)[(*jx)-(1)]*(*temp1)+(*y)[(*jy)-(1)]*(*temp2))
					(*ix) = (*jx)
					(*iy) = (*jy)
					for (*k) = (*kk) + 1; (*k) <= (*kk)+(*n)-(*j); (*k)++ {
						(*ix) = (*ix) + (*incx)
						(*iy) = (*iy) + (*incy)
						(*ap)[(*k)-(1)] = (*ap)[(*k)-(1)] + (*x)[(*ix)-(1)]*(*temp1) + (*y)[(*iy)-(1)]*(*temp2)
					//Label70:
					}
				} else {
					(*ap)[(*kk)-(1)] = (DBLE(((*ap)[(*kk)-(1)])))
				}
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				(*kk) = (*kk) + (*n) - (*j) + 1
			//Label80:
			}
		}
	}
	//*
	return
	//*
	//*     End of Zhpr2 .
	//*
}
