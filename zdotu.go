package golapack

// Zdotu ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 ZX(*),ZY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// ZDOTU forms the dot product of two complex vectors
//      ZDOTU = X^T * Y
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
//          ZX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of ZX
// \endverbatim
//
// \param[in] ZY
// \verbatim
//          ZY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
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
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Zdotu(n *int, zx *[]complex128, zxoff, incx *int, zy *[]complex128, zyoff, incy *int) (zdotuReturn complex128) {
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//        code for both increments equal to 1
		//
		for i = 1; i <= (*n); i++ {
			zdotuReturn = zdotuReturn + (*zx)[i-1+(*zxoff)]*(*zy)[i-1+(*zyoff)]
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
			zdotuReturn = zdotuReturn + (*zx)[ix-1+(*zxoff)]*(*zy)[iy-1+(*zyoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
	return
}
