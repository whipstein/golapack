package goblas

import (
	"math"
)

// Drotg ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drotg(da,db,c,s)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION c,da,db,s
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Drotg construct givens plane rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] da
// \verbatim
//          da is DOUBLE PRECISION
// \endverbatim
//
// \param[in] db
// \verbatim
//          db is DOUBLE PRECISION
// \endverbatim
//
// \param[out] c
// \verbatim
//          c is DOUBLE PRECISION
// \endverbatim
//
// \param[out] s
// \verbatim
//          s is DOUBLE PRECISION
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
// \endverbatim
//
//  =====================================================================
func Drotg(major *byte, da *float64, db *float64, c *float64, s *float64) {
	var r, roe, scale, z float64
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	roe = *db
	if math.Abs(*da) > math.Abs(*db) {
		roe = *da
	}
	scale = math.Abs(*da) + math.Abs(*db)
	if scale == 0.0 {
		*c = 1.0
		*s = 0.0
	} else {
		r = scale * math.Sqrt(math.Pow((*da)/scale, 2)+math.Pow((*db)/scale, 2))
		if roe < 0.0 {
			r = -r
		}
		*c = (*da) / r
		*s = (*db) / r
		z = 1.0
		if math.Abs(*da) > math.Abs(*db) {
			z = *s
		}
		if math.Abs(*db) >= math.Abs(*da) && *c != 0.0 {
			z = 1.0 / (*c)
		}
	}
	*da = r
	*db = z
}
