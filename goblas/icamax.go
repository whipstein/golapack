package goblas

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
//       INTEGER FUNCTION Icamax(n,cx,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX cx(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Icamax finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
//          cx is COMPLEX array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of cx
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
func Icamax(n *int, cx *[]complex64, incx *int) (icamaxReturn int) {
	var smax float32
	var i, ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n < 1 || *incx <= 0 {
		return
	}
	icamaxReturn = 1
	if *n == 1 {
		return
	}
	if *incx == 1 {
		//
		//        code for increment equal to 1
		//
		smax = scabs1(&(*cx)[0])
		for i = 2; i <= *n; i++ {
			if scabs1(&(*cx)[i-1]) > smax {
				icamaxReturn = i
				smax = scabs1(&(*cx)[i-1])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		smax = scabs1(&(*cx)[0])
		ix += *incx
		for i = 2; i <= *n; i++ {
			if scabs1(&(*cx)[ix-1]) > smax {
				icamaxReturn = i
				smax = scabs1(&(*cx)[ix-1])
			}
			ix += *incx
		}
	}
	return
}
