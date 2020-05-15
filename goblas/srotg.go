package goblas

// Srotg ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SROTG(sa,sb,c,s)
//
//       .. Scalar Arguments ..
//       REAL c,s,sa,sb
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    SROTG construct givens plane rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] sa
// \verbatim
//          sa is REAL
// \endverbatim
//
// \param[in] sb
// \verbatim
//          sb is REAL
// \endverbatim
//
// \param[out] c
// \verbatim
//          c is REAL
// \endverbatim
//
// \param[out] s
// \verbatim
//          s is REAL
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
//     jack dongarra, linpack, 3/11/78.
// \endverbatim
//
//  =====================================================================
func Srotg(major *byte, sa, sb, c, s *float32) {
	var r, roe, scale, z float32
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	roe = *sb
	if absf32(*sa) > absf32(*sb) {
		roe = *sa
	}
	scale = absf32(*sa) + absf32(*sb)
	if scale == 0.0 {
		*c = 1.0
		*s = 0.0
		r = 0.0
		z = 0.0
	} else {
		r = scale * sqrtf32(powf32((*sa)/scale, 2)+powf32((*sb)/scale, 2))
		r = signf32(1, roe) * r
		*c = (*sa) / r
		*s = (*sb) / r
		z = 1.0
		if absf32(*sa) > absf32(*sb) {
			z = (*s)
		}
		if absf32(*sb) >= absf32(*sa) && *c != 0.0 {
			z = 1.0 / (*c)
		}
	}
	*sa = r
	*sb = z
	return
}
