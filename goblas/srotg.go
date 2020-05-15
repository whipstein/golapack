package goblas

import "math"
import 
// \brief \b Srotg
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Srotg(SA,SB,C,S)
//
//       .. Scalar Arguments ..
//       REAL C,S,SA,SB
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Srotg construct givens plane rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SA
// \verbatim
//          SA is REAL
// \endverbatim
//
// \param[in] SB
// \verbatim
//          SB is REAL
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is REAL
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is REAL
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
func Srotg(sa *float64, sb *float64, c *float64, s *float64) {
	r := new(float64)
	roe := new(float64)
	scale := new(float64)
	z := new(float64)
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
	(*roe) = (*sb)
	if ABS((*sa)) > ABS(sb) {
		(*roe) = (*sa)
	}
	(*scale) = ABS((*sa)) + ABS((*sb))
	if (*scale) == 0.0 {
		(*c) = 1.0
		(*s) = 0.0
		(*r) = 0.0
		(*z) = 0.0
	} else {
		(*r) = (*scale) * SQRT(math.pow(((*sa)/(*scale)), 2)+math.pow(((*sb)/(*scale)), 2))
		(*r) = sign(func() *float64{y := 1.0; return &y }(), roe) * (*r)
		(*c) = (*sa) / (*r)
		(*s) = (*sb) / (*r)
		(*z) = 1.0
		if ABS((*sa)) > ABS(sb) {
			(*z) = (*s)
		}
		if abs ((*sb)) >= abs ((*sa)) && (*c) != 0.0 {
			(*z) = 1.0 / (*c)
		}
	}
	(*sa) = (*r)
	(*sb) = (*z)
	return
}
