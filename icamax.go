package golapack

// Icamax ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION ICAMAX(N,CX,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       COMPLEX CX(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Icamax(n *int, cx *[]complex64, cxoff, incx *int) (icamaxReturn int) {
	var smax float32
	var i, ix int

	icamaxReturn = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	icamaxReturn = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		smax = absc64((*cx)[0+(*cxoff)])
		for i = 2; i <= (*n); i++ {
			if absc64((*cx)[i-1+(*cxoff)]) > smax {
				icamaxReturn = i
				smax = absc64((*cx)[i-1+(*cxoff)])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		smax = absc64((*cx)[0+(*cxoff)])
		ix = ix + (*incx)
		for i = 2; i <= (*n); i++ {
			if absc64((*cx)[ix-1+(*cxoff)]) > smax {
				icamaxReturn = i
				smax = absc64((*cx)[ix-1+(*cxoff)])
			}
			ix = ix + (*incx)
		}
	}
	return
}
