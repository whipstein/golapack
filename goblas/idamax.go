package goblas

import "math"

// Idamax ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION Idamax(n,dx,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION dx(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Idamax finds the index of the first element having maximum absolute value.
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
// \param[in] dx
// \verbatim
//          dx is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of dx
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
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Idamax(major *byte, n *int, dx *[]float64, incx *int) (idamaxReturn int) {
	var dmax float64
	var i, ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	idamaxReturn = 0
	if *n < 1 || *incx <= 0 {
		return
	}
	idamaxReturn = 1
	if *n == 1 {
		return
	}
	if *incx == 1 {
		//
		//        code for increment equal to 1
		//
		dmax = math.Abs((*dx)[1-1])
		for i = 2; i <= *n; i++ {
			if math.Abs((*dx)[i-1]) > dmax {
				idamaxReturn = i
				dmax = math.Abs((*dx)[i-1])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		dmax = math.Abs((*dx)[1-1])
		ix += *incx
		for i = 2; i <= *n; i++ {
			if math.Abs((*dx)[ix-1]) > dmax {
				idamaxReturn = i
				dmax = (math.Abs(((*dx)[ix-1])))
			}
			ix += *incx
		}
	}
	return
}
