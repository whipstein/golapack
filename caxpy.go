package golapack

// Caxpy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
//
//       .. Scalar Arguments ..
//       COMPLEX CA
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       COMPLEX CX(*),CY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    CAXPY constant times a vector plus a vector.
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
//          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of CX
// \endverbatim
//
// \param[in,out] CY
// \verbatim
//          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
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
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Caxpy(n *int, ca *complex64, cx *[]complex64, cxoff, incx *int, cy *[]complex64, cyoff, incy *int) {
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if absc64(*ca) == 0.0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//        code for both increments equal to 1
		//
		for i = 1; i <= (*n); i++ {
			(*cy)[i-1+(*cyoff)] = (*cy)[i-1+(*cyoff)] + (*ca)*(*cx)[i-1+(*cxoff)]
		}
	} else {
		//
		//        code for unequal increments or equal increments
		//          not equal to 1
		//
		ix = 1
		iy = 1
		if (*incx) < 0 {
			ix = (-(*n)+1)*(*incx) + 1
		}
		if (*incy) < 0 {
			iy = (-(*n)+1)*(*incy) + 1
		}
		for i = 1; i <= (*n); i++ {
			(*cy)[iy-1+(*cyoff)] = (*cy)[iy-1+(*cyoff)] + (*ca)*(*cx)[ix-1+(*cxoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
