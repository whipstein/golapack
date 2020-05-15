package goblas

import 

// \brief \b Dcabs1
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dcabs1(Z)
//
//       .. Scalar Arguments ..
//       COMPLEX//16 Z
//       ..
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dcabs1 computes |Re(.)| + |Im(.)| of a double complex number
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] Z
// \verbatim
//          Z is COMPLEX//16
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
func Dcabs1(z *complex128) (dcabs1Return *float64) {
	dcabs1return := new(float64)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     ..
	//*  =====================================================================
	//*
	//*     .. Intrinsic Functions ..
	//*
	(*dcabs1) = ABS(DBLE((*z))) + ABS(DIMAG(z))
	return
}
