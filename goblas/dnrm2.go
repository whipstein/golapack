package goblas

import (
	"math"
)

// Dnrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dnrm2(n,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dnrm2 returns the euclidean dnrm2Return of a vector via the function
// name, so that
//
//    Dnrm2 := sqrt( x'*x )
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
//          x is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incx ) )
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
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to DLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Dnrm2(major *byte, n *int, x *[]float64, incx *int) (dnrm2Return float64) {
	var absxi, scale, ssq float64
	var ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n < 1 || *incx < 1 {
		dnrm2Return = 0.0
	} else if *n == 1 {
		dnrm2Return = math.Abs((*x)[0])
	} else {
		scale = 0.0
		ssq = 1.0
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL DLASSQ( n, x, incx, scale, ssq )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += *incx {
			if (*x)[ix-1] != 0.0 {
				absxi = math.Abs((*x)[ix-1])
				if scale < absxi {
					ssq = 1.0 + ssq*math.Pow(scale/absxi, 2)
					scale = absxi
				} else {
					ssq += math.Pow(absxi/scale, 2)
				}
			}
		}
		dnrm2Return = scale * math.Sqrt(ssq)
	}

	return
}
