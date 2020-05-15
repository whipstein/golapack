package goblas

import 

// \brief \b Ztpmv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ztpmv(UPLO,TRANS,DIAG,N,AP,X,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 AP(//),X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ztpmv  performs one of the matrix-vector operations
//
//    x := A//x,   or   x := A////T//x,   or   x := A////H//x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix, supplied in packed form.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
//           On entry, UPLO specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//              UPLO = 'L' or 'l'   A is a lower triangular matrix.
// \endverbatim
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER//1
//           On entry, TRANS specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   x := A//x.
//
//              TRANS = 'T' or 't'   x := A////T//x.
//
//              TRANS = 'C' or 'c'   x := A////H//x.
// \endverbatim
//
// \param[in] DIAG
// \verbatim
//          DIAG is CHARACTER//1
//           On entry, DIAG specifies whether or not A is unit
//           triangular as follows:
//
//              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
//
//              DIAG = 'N' or 'n'   A is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] AP
// \verbatim
//          AP is COMPLEX//16 array, dimension at least
//           ( ( n//( n + 1))/2).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular matrix packed sequentially,
//           column by column, so that AP1 contains a( 1, 1),
//           AP( 2) and AP( 3) contain a( 1, 2) and a( 2, 2)
//           respectively, and so on.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular matrix packed sequentially,
//           column by column, so that AP1 contains a( 1, 1),
//           AP( 2) and AP( 3) contain a( 2, 1) and a( 3, 1)
//           respectively, and so on.
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced, but are assumed to be unity.
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is COMPLEX//16 array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)).
//           Before entry, the incremented array X must contain the n
//           element vector x. On exit, X is overwritten with the
//           transformed vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           X. incx must not be zero.
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
func Ztpmv(uplo *byte, trans *byte, diag *byte, n *int, ap *[]complex128, x *[]complex128, incx *int) {
	zero := new(complex128)
	temp := new(complex128)
	i := new(int)
	info := new(int)
	ix := new(int)
	j := new(int)
	jx := new(int)
	k := new(int)
	kk := new(int)
	kx := new(int)
	noconj := new(bool)
	nounit := new(bool)
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
	} else if !Lsame((*trans), "n") && !Lsame((*trans), "t") && !Lsame((*trans), "c") {
		(*info) = 2
	} else if !Lsame((*diag), "u") && !Lsame((*diag), "n") {
		(*info) = 3
	} else if (*n) < 0 {
		(*info) = 4
	} else if (*incx) == 0 {
		(*info) = 7
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("ztpmv "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if (*n) == 0 {
		return
	}
	//*
	(*noconj) = (*Lsame(trans, func() *byte {y := byte('t'); return &y}()))
	(*nounit) = (*Lsame(diag, func() *byte {y := byte('n'); return &y}()))
	//*
	//*     Set up the start point in X if the increment is not unity. This
	//*     will be  ( N - 1)*incx  too small for descending loops.
	//*
	if (*incx) <= 0 {
		(*kx) = 1 - ((*n)-1)*(*incx)
	} else if (*incx) != 1 {
		(*kx) = 1
	}
	//*
	//*     Start the operations. In this version the elements of AP are
	//*     accessed sequentially with one pass through AP.
	//*
	if Lsame(trans, func() *byte {y := byte('n'); return &y}()) {
		//*
		//*        Form  x:= A*x.
		//*
		if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
			(*kk) = 1
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*x)[(*j)-1] != (*zero) {
						(*temp) = (*x)[(*j)-1]
						(*k) = (*kk)
						for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
							(*x)[(*i)-1] = (*x)[(*i)-1] + (*temp)*(*ap)[(*k)-1]
							(*k) = (*k) + 1
							//Label10:
						}
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] * (*ap)[(*kk)+(*j)-0]
						}
					}
					(*kk) = (*kk) + (*j)
					//Label20:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*x)[(*jx)-1] != (*zero) {
						(*temp) = (*x)[(*jx)-1]
						(*ix) = (*kx)
						for (*k) = (*kk); (*k) <= (*kk)+(*j)-2; (*k)++ {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] + (*temp)*(*ap)[(*k)-1]
							(*ix) = (*ix) + (*incx)
							//Label30:
						}
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] * (*ap)[(*kk)+(*j)-0]
						}
					}
					(*jx) = (*jx) + (*incx)
					(*kk) = (*kk) + (*j)
					//Label40:
				}
			}
		} else {
			(*kk) = ((*n) * ((*n) + 1)) / 2
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*x)[(*j)-1] != (*zero) {
						(*temp) = (*x)[(*j)-1]
						(*k) = (*kk)
						for (*i) = (*n); (*i) <= (*j)+1; (*i) += -1 {
							(*x)[(*i)-1] = (*x)[(*i)-1] + (*temp)*(*ap)[(*k)-1]
							(*k) = (*k) - 1
							//Label50:
						}
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] * (*ap)[(*kk)-(*n)+(*j)-1]
						}
					}
					(*kk) = (*kk) - ((*n) - (*j) + 1)
					//Label60:
				}
			} else {
				(*kx) = (*kx) + ((*n)-1)*(*incx)
				(*jx) = (*kx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*x)[(*jx)-1] != (*zero) {
						(*temp) = (*x)[(*jx)-1]
						(*ix) = (*kx)
						for (*k) = (*kk); (*k) <= (*kk)-((*n)-((*j)+1)); (*k) += -1 {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] + (*temp)*(*ap)[(*k)-1]
							(*ix) = (*ix) - (*incx)
							//Label70:
						}
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] * (*ap)[(*kk)-(*n)+(*j)-1]
						}
					}
					(*jx) = (*jx) - (*incx)
					(*kk) = (*kk) - ((*n) - (*j) + 1)
					//Label80:
				}
			}
		}
	} else {
		//*
		//*        Form  x := A**T*x  or  x := A**H*x.
		//*
		if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
			(*kk) = ((*n) * ((*n) + 1)) / 2
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*j)-1]
					(*k) = (*kk) - 1
					if *noconj {
						if *nounit {
							(*temp) = (*temp) * (*ap)[(*kk)-1]
						}
						for (*i) = (*j) - 1; (*i) <= 1; (*i) += -1 {
							(*temp) = (*temp) + (*ap)[(*k)-1]*(*x)[(*i)-1]
							(*k) = (*k) - 1
							//Label90:
						}
					} else {
						if *nounit {
							(*temp) = (*temp) * DCONJG(((*ap)[(*kk)-1]))
						}
						for (*i) = (*j) - 1; (*i) <= 1; (*i) += -1 {
							(*temp) = (*temp) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*i)-1]
							(*k) = (*k) - 1
							//Label100:
						}
					}
					(*x)[(*j)-1] = (*temp)
					(*kk) = (*kk) - (*j)
					//Label110:
				}
			} else {
				(*jx) = (*kx) + ((*n)-1)*(*incx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*jx)
					if *noconj {
						if *nounit {
							(*temp) = (*temp) * (*ap)[(*kk)-1]
						}
						for (*k) = (*kk) - 1; (*k) <= (*kk)-(*j)+1; (*k) += -1 {
							(*ix) = (*ix) - (*incx)
							(*temp) = (*temp) + (*ap)[(*k)-1]*(*x)[(*ix)-1]
							//Label120:
						}
					} else {
						if *nounit {
							(*temp) = (*temp) * DCONJG(((*ap)[(*kk)-1]))
						}
						for (*k) = (*kk) - 1; (*k) <= (*kk)-(*j)+1; (*k) += -1 {
							(*ix) = (*ix) - (*incx)
							(*temp) = (*temp) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*ix)-1]
							//Label130:
						}
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) - (*incx)
					(*kk) = (*kk) - (*j)
					//Label140:
				}
			}
		} else {
			(*kk) = 1
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*j)-1]
					(*k) = (*kk) + 1
					if *noconj {
						if *nounit {
							(*temp) = (*temp) * (*ap)[(*kk)-1]
						}
						for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
							(*temp) = (*temp) + (*ap)[(*k)-1]*(*x)[(*i)-1]
							(*k) = (*k) + 1
							//Label150:
						}
					} else {
						if *nounit {
							(*temp) = (*temp) * DCONJG(((*ap)[(*kk)-1]))
						}
						for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
							(*temp) = (*temp) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*i)-1]
							(*k) = (*k) + 1
							//Label160:
						}
					}
					(*x)[(*j)-1] = (*temp)
					(*kk) = (*kk) + ((*n) - (*j) + 1)
					//Label170:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*jx)
					if *noconj {
						if *nounit {
							(*temp) = (*temp) * (*ap)[(*kk)-1]
						}
						for (*k) = (*kk) + 1; (*k) <= (*kk)+(*n)-(*j); (*k)++ {
							(*ix) = (*ix) + (*incx)
							(*temp) = (*temp) + (*ap)[(*k)-1]*(*x)[(*ix)-1]
							//Label180:
						}
					} else {
						if *nounit {
							(*temp) = (*temp) * DCONJG(((*ap)[(*kk)-1]))
						}
						for (*k) = (*kk) + 1; (*k) <= (*kk)+(*n)-(*j); (*k)++ {
							(*ix) = (*ix) + (*incx)
							(*temp) = (*temp) + DCONJG(((*ap)[(*k)-1]))*(*x)[(*ix)-1]
							//Label190:
						}
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) + (*incx)
					(*kk) = (*kk) + ((*n) - (*j) + 1)
					//Label200:
				}
			}
		}
	}
	//*
	return
	//*
	//*     End of Ztpmv .
	//*
}
