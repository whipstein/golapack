package golapack

// Izamax ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION IZAMAX(N,ZX,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 ZX(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
//          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of ZX
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
// \ingroup aux_blas
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, 1/15/85.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Izamax(n *int, zx *[]complex128, zxoff, incx *int) (izamaxReturn int) {
	var dmax float64
	var i, ix int

	izamaxReturn = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	izamaxReturn = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		dmax = abssumc128((*zx)[0+(*zxoff)])
		for i = 2; i <= (*n); i++ {
			if abssumc128((*zx)[i-1+(*zxoff)]) > dmax {
				izamaxReturn = i
				dmax = abssumc128((*zx)[i-1+(*zxoff)])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		dmax = abssumc128((*zx)[0+(*zxoff)])
		ix = ix + (*incx)
		for i = 2; i <= (*n); i++ {
			if abssumc128((*zx)[ix-1+(*zxoff)]) > dmax {
				izamaxReturn = i
				dmax = abssumc128((*zx)[ix-1+(*zxoff)])
			}
			ix = ix + (*incx)
		}
	}
	return
}
