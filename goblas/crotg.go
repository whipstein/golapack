package goblas

import "math"
import 
// \brief \b Crotg
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Crotg(CA,CB,C,S)
//
//       .. Scalar Arguments ..
//       COMPLEX CA,CB,S
//       REAL C
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Crotg determines a complex Givens rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] CA
// \verbatim
//          CA is COMPLEX
// \endverbatim
//
// \param[in] CB
// \verbatim
//          CB is COMPLEX
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is REAL
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is COMPLEX
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
// \ingroup complex_blas_level1
//
//  =====================================================================
func Crotg(ca *complex128, cb *complex128, c *float64, s *complex128) {
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
	if cabs ((*ca)) == 0. {
		(*c) = 0.
		(*s) = (1. + (0.)*1i)
		(*ca) = (*cb)
	} else {
		(*scale) = CABS((*ca)) + CABS((*cb))
		(*norm) = (*scale) * SQRT(math.pow((CABS((*ca)/(*scale))), 2)+math.pow((CABS((*cb)/(*scale))), 2))
		(*alpha) = (*ca) / CABS((*ca))
		(*c) = CABS((*ca)) / (*norm)
		(*s) = (*alpha) * CONJG((*cb)) / (*norm)
		(*ca) = (*alpha) * complex((*norm), (0))
	}
	return
}
