package goblas

// \brief \b Caxpy
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Caxpy(N,CA,CX,incx,CY,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX CA
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       COMPLEX CX(//),CY(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Caxpy constant times a vector plus a vector.
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
// \param[in] CA
// \verbatim
//          CA is COMPLEX
//           On entry, CA specifies the scalar alpha.
// \endverbatim
//
// \param[in] CX
// \verbatim
//          CX is COMPLEX array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of CX
// \endverbatim
//
// \param[in,out] CY
// \verbatim
//          CY is COMPLEX array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of CY
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
// \ingroup complex_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Caxpy(n *int, ca *complex128, cx *[]complex128, incx *int, cy *[]complex128, incy *int) {
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
	//*     .. External Functions ..
	//*     ..
	if (*n) <= 0 {
		return
	}
	if Scabs1((*ca)) == 0.0e+0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*        code for both increments equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*cy)[(*i)-1] = (*cy)[(*i)-1] + (*ca)*(*cx)[(*i)-1]
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
			(*cy)[(*iy)-1] = (*cy)[(*iy)-1] + (*ca)*(*cx)[(*ix)-1]
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	//*
	return
}
