package golapack

// Snrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION SNRM2(N,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       REAL X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SNRM2 returns the euclidean norm of a vector via the function
// name, so that
//
//    SNRM2 := sqrt( x'*x ).
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
//          X is REAL array, dimension ( 1 + ( N - 1 )*absf32( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of SX
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
//     Modified on 14-October-1993 to inline the call to SLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Snrm2(n *int, x *[]float32, xoff, incx *int) (snrm2Return float32) {
	var absxi, norm, scale, ssq float32
	var ix int
	var one float32 = 1.0
	var zero float32 = 0.0

	if (*n) < 1 || (*incx) < 1 {
		norm = zero
	} else if (*n) == 1 {
		norm = absf32((*x)[0+(*xoff)])
	} else {
		scale = zero
		ssq = one
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += (*incx) {
			if (*x)[ix-1+(*xoff)] != zero {
				absxi = absf32((*x)[ix-1+(*xoff)])
				if scale < absxi {
					ssq = one + ssq*powf32(scale/absxi, 2)
					scale = absxi
				} else {
					ssq = ssq + powf32(absxi/scale, 2)
				}
			}
		}
		norm = scale * sqrtf32(ssq)
	}

	snrm2Return = norm
	return
}
