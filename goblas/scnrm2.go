package goblas

import (
	"math"
)

// Scnrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION Scnrm2(n,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Scnrm2 returns the euclidean norm of a vector via the function
// name, so that
//
//    Scnrm2 := sqrt( x**H*x )
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
//          x is COMPLEX array, dimension n
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
// \ingroup single_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to CLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Scnrm2(n *int, x *[]complex64, incx *int) (norm float32) {
	var scale, ssq, temp float32
	var ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n < 1 || *incx < 1 {
		norm = 0.0
	} else {
		scale = 0.0
		ssq = 1.0
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL CLASSQ( n, x, incx, scale, ssq )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += *incx {
			if real((*x)[ix-1]) != 0.0 {
				temp = float32(math.Abs(float64(real((*x)[ix-1]))))
				if scale < temp {
					ssq = 1.0 + ssq*float32(math.Pow(float64(scale/temp), 2))
					scale = temp
				} else {
					ssq += float32(math.Pow(float64(temp/scale), 2))
				}
			}
			if imag((*x)[ix-1]) != 0.0 {
				temp = float32(math.Abs(float64(imag((*x)[ix-1]))))
				if scale < temp {
					ssq = 1.0 + ssq*float32(math.Pow(float64(scale/temp), 2))
					scale = temp
				} else {
					ssq += float32(math.Pow(float64(temp/scale), 2))
				}
			}
		}
		norm = scale * float32(math.Sqrt(float64(ssq)))
	}

	return
}
