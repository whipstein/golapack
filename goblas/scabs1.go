package goblas

import "math"

// scabs1 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION scabs1(Z)
//
//       .. Scalar Arguments ..
//       COMPLEX Z
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// scabs1 computes |Re(.)| + |Im(.)| of a complex number
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] Z
// \verbatim
//          Z is COMPLEX
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
//  =====================================================================
func scabs1(z *complex64) float32 {
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	return float32(math.Abs(float64(real(*z)) + math.Abs(float64(imag(*z)))))
}
