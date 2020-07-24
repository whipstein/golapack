package golapack

// Cswap ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
//
//       .. Scalar Arguments ..
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
//   CSWAP interchanges two vectors.
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
// \param[in,out] CX
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
func Cswap(n *int, cx *[]complex64, cxoff, incx *int, cy *[]complex64, cyoff, incy *int) {
	var ctemp complex64
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//       code for both increments equal to 1
		for i = 1; i <= (*n); i++ {
			ctemp = (*cx)[i-1+(*cxoff)]
			(*cx)[i-1+(*cxoff)] = (*cy)[i-1+(*cyoff)]
			(*cy)[i-1+(*cyoff)] = ctemp
		}
	} else {
		//
		//       code for unequal increments or equal increments not equal
		//         to 1
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
			ctemp = (*cx)[ix-1+(*cxoff)]
			(*cx)[ix-1+(*cxoff)] = (*cy)[iy-1+(*cyoff)]
			(*cy)[iy-1+(*cyoff)] = ctemp
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
