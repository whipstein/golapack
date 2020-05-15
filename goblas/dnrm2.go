package goblas

import "math"
import 

// \brief \b Dnrm2
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dnrm2(N,X,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dnrm2 returns the euclidean norm of a vector via the function
// name, so that
//
//    Dnrm2 := sqrt( x'//x)
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
//          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of DX
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
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to DLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Dnrm2(n *int, x *[]float64, incx *int) (dnrm2Return *float64) {
	dnrm2return := new(float64)
	one := new(float64)
	zero := new(float64)
	absxi := new(float64)
	norm := new(float64)
	scale := new(float64)
	ssq := new(float64)
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
	} else if (*n) == 1 {
		(*norm) = (ABS(((*x)[0])))
	} else {
		(*scale) = (*zero)
		(*ssq) = (*one)
		//*        The following loop is equivalent to this call to the LAPACK
		//*        auxiliary routine:
		//*        CALL DLASSQ( N, X, incx, SCALE, SSQ)
		//*
		for (*ix) = 1; (*ix) <= 1+((*n)-1)*(*incx); (*ix) += (*incx) {
			if (*x)[(*ix)-1] != (*zero) {
				(*absxi) = (ABS(((*x)[(*ix)-1])))
				if (*scale) < (*absxi) {
					(*ssq) = (*one) + (*ssq)*math.pow(((*scale)/(*absxi)), 2)
					(*scale) = (*absxi)
				} else {
					(*ssq) = (*ssq) + math.pow(((*absxi)/(*scale)), 2)
				}
			}
			//Label10:
		}
		(*norm) = (*scale) * SQRT((*ssq))
	}
	//*
	(*dnrm2) = (*norm)
	return
	//*
	//*     End of Dnrm2.
	//*
}
