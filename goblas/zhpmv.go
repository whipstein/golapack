package goblas

import 

// \brief \b Zhpmv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zhpmv(UPLO,N,ALPHA,AP,X,incx,BETA,Y,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX//16 ALPHA,BETA
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
// Zhpmv  performs the matrix-vector operation
//
//    y := alpha//A//x + beta//y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian matrix, supplied in packed form.
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
// \param[in] AP
// \verbatim
//          AP is COMPLEX//16 array, dimension at least
//           ( ( n//( n + 1))/2).
//           Before entry with UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular part of the hermitian matrix
//           packed sequentially, column by column, so that AP1
//           contains a( 1, 1), AP( 2) and AP( 3) contain a( 1, 2)
//           and a( 2, 2) respectively, and so on.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular part of the hermitian matrix
//           packed sequentially, column by column, so that AP1
//           contains a( 1, 1), AP( 2) and AP( 3) contain a( 2, 1)
//           and a( 3, 1) respectively, and so on.
//           Note that the imaginary parts of the diagonal elements need
//           not be set and are assumed to be zero.
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
// \param[in] BETA
// \verbatim
//          BETA is COMPLEX//16
//           On entry, BETA specifies the scalar beta. When BETA is
//           supplied as zero then Y need not be set on input.
// \endverbatim
//
// \param[in,out] Y
// \verbatim
//          Y is COMPLEX//16 array, dimension at least
//           ( 1 + ( n - 1)//abs( incy)).
//           Before entry, the incremented array Y must contain the n
//           element vector y. On exit, Y is overwritten by the updated
//           vector y.
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
// \ingroup complex16_blas_level2
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
func Zhpmv(uplo *byte, n *int, alpha *complex128, ap *[]complex128, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	one := new(complex128)
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
	(*one) = (1.0e+0 + (0.0e+0)*1i)
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
		(*info) = 6
	} else if (*incy) == 0 {
		(*info) = 9
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("zhpmv "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*n) == 0) || (((*alpha) == (*zero)) && ((*beta) == (*one))) {
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
	//*     Start the operations. In this version the elements of the array AP
	//*     are accessed sequentially with one pass through AP.
	//*
	//*     First form  y := beta*y.
	//*
	if (*beta) != (*one) {
		if (*incy) == 1 {
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*i)-1] = (*zero)
					//Label10:
				}
			} else {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*i)-1] = (*beta) * (*y)[(*i)-1]
					//Label20:
				}
			}
		} else {
			(*iy) = (*ky)
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*y)[(*iy)-1] = (*zero)
					(*iy) = (*iy) + (*incy)
					//Label30:
				}
			} else {
				for (*i) = 1; (*i) <= (*n); (*i)++ {
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
	(*kk) = 1
	if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
		//*
		//*        Form  y  when AP contains the upper triangle.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*j)-1]
				(*temp2) = (*zero)
				(*k) = (*kk)
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*y)[(*i)-1] = (*y)[(*i)-1] + (*temp1)*(*ap)[(*k)-1]
					(*temp2) = (*temp2) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*i)-1]
					(*k) = (*k) + 1
					//Label50:
				}
				(*y)[(*j)-1] = (*y)[(*j)-1] + (*temp1)*DBLE(((*ap)[(*kk)+(*j)-0])) + (*alpha)*(*temp2)
				(*kk) = (*kk) + (*j)
				//Label60:
			}
		} else {
			(*jx) = (*kx)
			(*jy) = (*ky)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*jx)-1]
				(*temp2) = (*zero)
				(*ix) = (*kx)
				(*iy) = (*ky)
				for (*k) = (*kk); (*k) <= (*kk)+(*j)-2; (*k)++ {
					(*y)[(*iy)-1] = (*y)[(*iy)-1] + (*temp1)*(*ap)[(*k)-1]
					(*temp2) = (*temp2) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*ix)-1]
					(*ix) = (*ix) + (*incx)
					(*iy) = (*iy) + (*incy)
					//Label70:
				}
				(*y)[(*jy)-1] = (*y)[(*jy)-1] + (*temp1)*DBLE(((*ap)[(*kk)+(*j)-0])) + (*alpha)*(*temp2)
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				(*kk) = (*kk) + (*j)
				//Label80:
			}
		}
	} else {
		//*
		//*        Form  y  when AP contains the lower triangle.
		//*
		if ((*incx) == 1) && ((*incy) == 1) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*j)-1]
				(*temp2) = (*zero)
				(*y)[(*j)-1] = (*y)[(*j)-1] + (*temp1)*DBLE(((*ap)[(*kk)-1]))
				(*k) = (*kk) + 1
				for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
					(*y)[(*i)-1] = (*y)[(*i)-1] + (*temp1)*(*ap)[(*k)-1]
					(*temp2) = (*temp2) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*i)-1]
					(*k) = (*k) + 1
					//Label90:
				}
				(*y)[(*j)-1] = (*y)[(*j)-1] + (*alpha)*(*temp2)
				(*kk) = (*kk) + ((*n) - (*j) + 1)
				//Label100:
			}
		} else {
			(*jx) = (*kx)
			(*jy) = (*ky)
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				(*temp1) = (*alpha) * (*x)[(*jx)-1]
				(*temp2) = (*zero)
				(*y)[(*jy)-1] = (*y)[(*jy)-1] + (*temp1)*DBLE(((*ap)[(*kk)-1]))
				(*ix) = (*jx)
				(*iy) = (*jy)
				for (*k) = (*kk) + 1; (*k) <= (*kk)+(*n)-(*j); (*k)++ {
					(*ix) = (*ix) + (*incx)
					(*iy) = (*iy) + (*incy)
					(*y)[(*iy)-1] = (*y)[(*iy)-1] + (*temp1)*(*ap)[(*k)-1]
					(*temp2) = (*temp2) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*ix)-1]
					//Label110:
				}
				(*y)[(*jy)-1] = (*y)[(*jy)-1] + (*alpha)*(*temp2)
				(*jx) = (*jx) + (*incx)
				(*jy) = (*jy) + (*incy)
				(*kk) = (*kk) + ((*n) - (*j) + 1)
				//Label120:
			}
		}
	}
	//*
	return
	//*
	//*     End of Zhpmv .
	//*
}
