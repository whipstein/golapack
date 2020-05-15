package goblas

import "math"

// Dasum ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dasum(n,dx,incx)
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
//    Dasum takes the sum of the absolute values.
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
// \ingroup double_blas_level1
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
func Dasum(major *byte, n *int, dx *[]float64, incx *int) (dasumReturn float64) {
	var i, m, mp1, nincx int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 || *incx <= 0 {
		return
	}
	if *incx == 1 {
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 6
		if m != 0 {
			for i = 0; i < m; i++ {
				dasumReturn += math.Abs((*dx)[i])
			}
			if *n < 6 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 6 {
			dasumReturn += math.Abs((*dx)[i-1]) + math.Abs((*dx)[i]) + math.Abs((*dx)[i+1]) + math.Abs((*dx)[i+2]) + math.Abs((*dx)[i+3]) + math.Abs((*dx)[i+4])
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 0; i < nincx; i += *incx {
			dasumReturn += math.Abs(((*dx)[i]))
		}
	}
	return
}
