package goblas

import "math/cmplx"

// Cdotc ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       COMPLEX FUNCTION Cdotc(n,cx,incx,cy,incy)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       COMPLEX cx(*),cy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Cdotc forms the dot product of two complex vectors
//      Cdotc = x^H * y
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] cx
// \verbatim
//          cx is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of cx
// \endverbatim
//
// \param[in] cy
// \verbatim
//          cy is REAL array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of cy
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
//     jack dongarra, linpack,  3/11/78.
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Cdotc(n *int, cx *[]complex64, incx *int, cy *[]complex64, incy *int) (ctemp complex64) {
	var i, ix, iy int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//        code for both increments equal to 1
		//
		for i = 1; i <= *n; i++ {
			ctemp += complex64(cmplx.Conj(complex128((*cx)[i-1]))) * (*cy)[i-1]
		}
	} else {
		//
		//        code for unequal increments or equal increments
		//          not equal to 1
		//
		ix = 1
		iy = 1
		if *incx < 0 {
			ix = (-(*n)+1)*(*incx) + 1
		}
		if *incy < 0 {
			iy = (-(*n)+1)*(*incy) + 1
		}
		for i = 1; i <= *n; i++ {
			ctemp += complex64(cmplx.Conj(complex128((*cx)[ix-1]))) * (*cy)[iy-1]
			ix += *incx
			iy += *incy
		}
	}
	return
}
