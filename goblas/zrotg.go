package goblas

import "math"
import 
// \brief \b Zrotg
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zrotg(CA,CB,C,S)
//
//       .. Scalar Arguments ..
//       COMPLEX//16 CA,CB,S
//       DOUBLE PRECISION C
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Zrotg determines a double complex Givens rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] CA
// \verbatim
//          CA is COMPLEX//16
// \endverbatim
//
// \param[in] CB
// \verbatim
//          CB is COMPLEX//16
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is COMPLEX//16
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
// \ingroup complex16_blas_level1
//
//  =====================================================================
func Zrotg(ca *complex128, cb *complex128, c *float64, s *complex128) {
	alpha := new(complex128)
	norm := new(float64)
	scale := new(float64)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	if cdabs ((*ca)) == 0.0 {
		(*c) = 0.0
		(*s) = (1.0 + (0.0)*1i)
		(*ca) = (*cb)
	} else {
		(*scale) = CDABS(ca) + CDABS(cb)
		(*norm) = (*scale) * DSQRT(math.pow((CDABS((*ca)/DCMPLX(scale, func() *float64{y := 0.0; return &y}()))), 2)+math.pow((CDABS((*cb)/DCMPLX(scale, func() *float64{y := 0.0; return &y}()))), 2))
		(*alpha) = (*ca) / CDABS(ca)
		(*c) = CDABS(ca) / (*norm)
		(*s) = (*alpha) * DCONJG((*cb)) / (*norm)
		(*ca) = (*alpha) * complex((*norm), (0))
	}
	return
}
