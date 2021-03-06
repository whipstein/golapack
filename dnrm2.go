package golapack

// Dnrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// DNRM2 returns the euclidean norm of a vector via the function
// name, so that
//
//    DNRM2 := sqrt( x'*x )
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
//          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*absf64( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
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
func Dnrm2(n *int, x *[]float64, xoff, incx *int) (dnrm2Return float64) {
	var absxi, norm, scale, ssq float64
	var ix int
	var one float64 = 1.0
	var zero float64 = 0.0

	if (*n) < 1 || (*incx) < 1 {
		norm = zero
	} else if (*n) == 1 {
		norm = absf64((*x)[0+(*xoff)])
	} else {
		scale = zero
		ssq = one
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += (*incx) {
			if (*x)[ix-1+(*xoff)] != zero {
				absxi = absf64((*x)[ix-1+(*xoff)])
				if scale < absxi {
					ssq = one + ssq*powf64(scale/absxi, 2)
					scale = absxi
				} else {
					ssq = ssq + powf64(absxi/scale, 2)
				}
			}
		}
		norm = scale * sqrtf64(ssq)
	}
	dnrm2Return = norm
	return
}
