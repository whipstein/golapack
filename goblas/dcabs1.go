package goblas

import "math"

// dcabs1 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION dcabs1(Z)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 Z
//       ..
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// dcabs1 computes |Re(.)| + |Im(.)| of a double complex number
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] Z
// \verbatim
//          Z is COMPLEX*16
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
//  =====================================================================
func dcabs1(z *complex128) float64 {
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	return math.Abs(real(*z)) + math.Abs(imag(*z))
}
