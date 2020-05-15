package goblas

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
//       INTEGER FUNCTION Izamax(n,zx,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 zx(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Izamax finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
// \param[in] zx
// \verbatim
//          zx is COMPLEX*16 array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of zx
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
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Izamax(n *int, zx *[]complex128, incx *int) (izamaxReturn int) {
	var dmax float64
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
	izamaxReturn = 1
	if *n == 1 {
		return
	}
	if *incx == 1 {
		//
		//        code for increment equal to 1
		//
		dmax = dcabs1(&(*zx)[0])
		for i = 2; i <= *n; i++ {
			if dcabs1(&(*zx)[i-1]) > dmax {
				izamaxReturn = i
				dmax = dcabs1(&(*zx)[i-1])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		dmax = dcabs1(&(*zx)[0])
		ix += *incx
		for i = 2; i <= *n; i++ {
			if dcabs1(&(*zx)[ix-1]) > dmax {
				izamaxReturn = i
				dmax = dcabs1(&(*zx)[ix-1])
			}
			ix += *incx
		}
	}
	return
}
