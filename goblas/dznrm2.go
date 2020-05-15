package goblas

import (
	"math"
)

// Dznrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dznrm2(n,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dznrm2 returns the euclidean dznrm2Return of a vector via the function
// name, so that
//
//    Dznrm2 := sqrt( x**H*x )
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
// \param[in] x
// \verbatim
//          x is COMPLEX*16 array, dimension n
//         complex vector with n elements
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of x
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
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to ZLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Dznrm2(n *int, x *[]complex128, incx *int) (dznrm2Return float64) {
	var scale, ssq, temp float64
	var ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n < 1 || *incx < 1 {
		dznrm2Return = 0.0
	} else {
		scale = 0.0
		ssq = 1.0
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL ZLASSQ( n, x, incx, scale, ssq )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += *incx {
			if real((*x)[ix-1]) != 0.0 {
				temp = math.Abs(real((*x)[ix-1]))
				if scale < temp {
					ssq = 1.0 + ssq*math.Pow(scale/temp, 2)
					scale = temp
				} else {
					ssq += math.Pow(temp/scale, 2)
				}
			}
			if imag((*x)[ix-1]) != 0.0 {
				temp = math.Abs(imag((*x)[ix-1]))
				if scale < temp {
					ssq = 1.0 + ssq*math.Pow(scale/temp, 2)
					scale = temp
				} else {
					ssq += math.Pow(temp/scale, 2)
				}
			}
		}
		dznrm2Return = scale * math.Sqrt(ssq)
	}

	return
}
