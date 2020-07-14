package golapack

// Srot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
//
//       .. Scalar Arguments ..
//       REAL C,S
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       REAL SX(*),SY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    applies a plane rotation.
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
// \param[in,out] SX
// \verbatim
//          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of SX
// \endverbatim
//
// \param[in,out] SY
// \verbatim
//          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of SY
// \endverbatim
//
// \param[in] C
// \verbatim
//          C is REAL
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is REAL
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
// \ingroup single_blas_level1
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
func Srot(n *int, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int, c *float32, s *float32) {
	var stemp float32
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//       code for both increments equal to 1
		//
		for i = 1; i <= (*n); i++ {
			stemp = (*c)*(*sx)[i-1+(*sxoff)] + (*s)*(*sy)[i-1+(*syoff)]
			(*sy)[i-1+(*syoff)] = (*c)*(*sy)[i-1+(*syoff)] - (*s)*(*sx)[i-1+(*sxoff)]
			(*sx)[i-1+(*sxoff)] = stemp
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
			stemp = (*c)*(*sx)[ix-1+(*sxoff)] + (*s)*(*sy)[iy-1+(*syoff)]
			(*sy)[iy-1+(*syoff)] = (*c)*(*sy)[iy-1+(*syoff)] - (*s)*(*sx)[ix-1+(*sxoff)]
			(*sx)[ix-1+(*sxoff)] = stemp
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
