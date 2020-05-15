package goblas

import 

// \brief \b Zdotc
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       COMPLEX//16 FUNCTION Zdotc(N,ZX,incx,ZY,incy)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 ZX(//),ZY(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zdotc forms the dot product of two complex vectors
//      Zdotc = X^H // Y
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] ZX
// \verbatim
//          ZX is REAL array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of ZX
// \endverbatim
//
// \param[in] ZY
// \verbatim
//          ZY is REAL array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of ZY
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
// \date November 2017
//
// \ingroup complex16_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, 3/11/78.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Zdotc(n *int, zx *[]complex128, incx *int, zy *[]complex128, incy *int) (zdotcReturn *complex128) {
	zdotcreturn := new(complex128)
	ztemp := new(complex128)
	i := new(int)
	ix := new(int)
	iy := new(int)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	(*ztemp) = (0.0 + (0.0)*1i)
	(*zdotc) = (0.0 + (0.0)*1i)
	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*        code for both increments equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ztemp) = (*ztemp) + DCONJG(((*zx)[(*i)-1]))*(*zy)[(*i)-1]
		}
	} else {
		//*
		//*        code for unequal increments or equal increments
		//*          not equal to 1
		//*
		(*ix) = 1
		(*iy) = 1
		if (*incx) < 0 {
			(*ix) = (-(*n)+1)*(*incx) + 1
		}
		if (*incy) < 0 {
			(*iy) = (-(*n)+1)*(*incy) + 1
		}
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ztemp) = (*ztemp) + DCONJG(((*zx)[(*ix)-1]))*(*zy)[(*iy)-1]
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	(*zdotc) = (*ztemp)
	return
}
