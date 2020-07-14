package golapack

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
//       SUBROUTINE DROTG(DA,DB,C,S)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION C,DA,DB,S
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    DROTG construct givens plane rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] DA
// \verbatim
//          DA is DOUBLE PRECISION
// \endverbatim
//
// \param[in] DB
// \verbatim
//          DB is DOUBLE PRECISION
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION
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
func Drotg(da *float64, db *float64, c *float64, s *float64) {
	var r, roe, scale, z float64

	roe = (*db)
	if absf64(*da) > absf64(*db) {
		roe = (*da)
	}
	scale = absf64(*da) + absf64(*db)
	if scale == 0.0 {
		(*c) = 1.0
		(*s) = 0.0
		r = 0.0
		z = 0.0
	} else {
		r = scale * sqrtf64(powf64((*da)/scale, 2)+powf64((*db)/scale, 2))
		r = signf64(1.0, roe) * r
		(*c) = (*da) / r
		(*s) = (*db) / r
		z = 1.0
		if absf64(*da) > absf64(*db) {
			z = (*s)
		}
		if absf64(*db) >= absf64(*da) && (*c) != 0.0 {
			z = 1.0 / (*c)
		}
	}
	(*da) = r
	(*db) = z
}
