package goblas

import 
import "math"

// \brief \b Scnrm2
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION Scnrm2(N,X,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       COMPLEX X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Scnrm2 returns the euclidean norm of a vector via the function
// name, so that
//
//    Scnrm2 := sqrt( x////H//x)
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is COMPLEX array, dimension (N)
//         complex vector with N elements
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of X
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
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to CLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Scnrm2(n *int, x *[]complex128, incx *int) (scnrm2Return *float64) {
	scnrm2return := new(float64)
	one := new(float64)
	zero := new(float64)
	norm := new(float64)
	scale := new(float64)
	ssq := new(float64)
	temp := new(float64)
	ix := new(int)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. Parameters ..
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
	//*     ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	if (*n) < 1 || (*incx) < 1 {
		(*norm) = (*zero)
	} else {
		(*scale) = (*zero)
		(*ssq) = (*one)
		//*        The following loop is equivalent to this call to the LAPACK
		//*        auxiliary routine:
		//*        CALL CLASSQ( N, X, incx, SCALE, SSQ)
		//*
		for (*ix) = 1; (*ix) <= 1+((*n)-1)*(*incx); (*ix) += (*incx) {
			if real(((*x)[(*ix)-1])) != (*zero) {
				(*temp) = (ABS(real(((*x)[(*ix)-1]))))
				if (*scale) < (*temp) {
					(*ssq) = (*one) + (*ssq)*math.pow(((*scale)/(*temp)), 2)
					(*scale) = (*temp)
				} else {
					(*ssq) = (*ssq) + math.pow(((*temp)/(*scale)), 2)
				}
			}
			if imag(((*x)[(*ix)-1])) != (*zero) {
				(*temp) = (ABS(imag(((*x)[(*ix)-1]))))
				if (*scale) < (*temp) {
					(*ssq) = (*one) + (*ssq)*math.pow(((*scale)/(*temp)), 2)
					(*scale) = (*temp)
				} else {
					(*ssq) = (*ssq) + math.pow(((*temp)/(*scale)), 2)
				}
			}
			//Label10:
		}
		(*norm) = (*scale) * SQRT((*ssq))
	}
	//*
	(*scnrm2) = (*norm)
	return
	//*
	//*     End of Scnrm2.
	//*
}
