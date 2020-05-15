package goblas

import "math"
import 

// \brief \b Drotg
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drotg(DA,DB,C,S)
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
//    Drotg construct givens plane rotation.
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
	(*roe) = (*db)
	if ABS((*da)) > dabs(db) {
		(*roe) = (*da)
	}
	(*scale) = ABS((*da)) + ABS((*db))
	if (*scale) == 0.0 {
		(*c) = 1.0
		(*s) = 0.0
		(*r) = 0.0
		(*z) = 0.0
	} else {
		(*r) = (*scale) * DSQRT(math.pow(((*da)/(*scale)), 2)+math.pow(((*db)/(*scale)), 2))
		(*r) = dsign(func() *float64 {y := 1.0; return &y}(), roe) * (*r)
		(*c) = (*da) / (*r)
		(*s) = (*db) / (*r)
		(*z) = 1.0
		if ABS((*da)) > dabs(db) {
			(*z) = (*s)
		}
		if dabs((*db)) >= dabs((*da)) && (*c) != 0.0 {
			(*z) = 1.0 / (*c)
		}
	}
	(*da) = (*r)
	(*db) = (*z)
	return
}
