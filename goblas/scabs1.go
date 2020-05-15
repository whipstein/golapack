package goblas

import 

// \brief \b Scabs1
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION Scabs1(Z)
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
// Scabs1 computes |Re(.)| + |Im(.)| of a complex number
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
func Scabs1(z *complex128) (scabs1Return *float64) {
	scabs1return := new(float64)
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
	//*     .. Intrinsic Functions ..
	//*     ..
	(*scabs1) = ABS(real((*z))) + ABS(imag((*z)))
	return
}
