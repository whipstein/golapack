package golapack

// Dznrm2 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// DZNRM2 returns the euclidean norm of a vector via the function
// name, so that
//
//    DZNRM2 := sqrt( x**H*x )
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
//          X is COMPLEX*16 array, dimension (N)
//         complex vector with N elements
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
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
// \ingroup double_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to ZLASSQ.
//     Sven Hammarling, Nag Ltd.
// \endverbatim
//
//  =====================================================================
func Dznrm2(n *int, x *[]complex128, xoff, incx *int) (dznrm2Return float64) {
	var norm, scale, ssq, temp float64
	var ix int
	var one float64 = 1.0
	var zero float64 = 0.0

	if (*n) < 1 || (*incx) < 1 {
		norm = zero
	} else {
		scale = zero
		ssq = one
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += (*incx) {
			if real((*x)[ix-1+(*xoff)]) != zero {
				temp = absf64(real((*x)[ix-1+(*xoff)]))
				if scale < temp {
					ssq = one + ssq*powf64(scale/temp, 2)
					scale = temp
				} else {
					ssq = ssq + powf64(temp/scale, 2)
				}
			}
			if imag((*x)[ix-1+(*xoff)]) != zero {
				temp = absf64(imag((*x)[ix-1+(*xoff)]))
				if scale < temp {
					ssq = one + ssq*powf64(scale/temp, 2)
					scale = temp
				} else {
					ssq = ssq + powf64(temp/scale, 2)
				}
			}
		}
		norm = scale * sqrtf64(ssq)
	}

	return norm
}
