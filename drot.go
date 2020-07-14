package golapack

// Drot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION C,S
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DX(*),DY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    DROT applies a plane rotation.
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
// \param[in,out] DX
// \verbatim
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of DX
// \endverbatim
//
// \param[in,out] DY
// \verbatim
//          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of DY
// \endverbatim
//
// \param[in] C
// \verbatim
//          C is DOUBLE PRECISION
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is DOUBLE PRECISION
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
// \ingroup double_blas_level1
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
func Drot(n *int, dx *[]float64, dxoff, incx *int, dy *[]float64, dyoff, incy *int, c *float64, s *float64) {
	var dtemp float64
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//       code for both increments equal to 1
		//
		for i = 1; i <= (*n); i++ {
			dtemp = (*c)*(*dx)[i-1+(*dxoff)] + (*s)*(*dy)[i-1+(*dyoff)]
			(*dy)[i-1+(*dyoff)] = (*c)*(*dy)[i-1+(*dyoff)] - (*s)*(*dx)[i-1+(*dxoff)]
			(*dx)[i-1+(*dxoff)] = dtemp
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
			dtemp = (*c)*(*dx)[ix-1+(*dxoff)] + (*s)*(*dy)[iy-1+(*dyoff)]
			(*dy)[iy-1+(*dyoff)] = (*c)*(*dy)[iy-1+(*dyoff)] - (*s)*(*dx)[ix-1+(*dxoff)]
			(*dx)[ix-1+(*dxoff)] = dtemp
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
