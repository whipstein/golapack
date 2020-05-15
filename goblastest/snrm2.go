package main

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
//          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
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
func Snrm2(n int, x *[]float32, incx int) float32 {
	var absxi, norm, scale, ssq float32
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	if n < 1 || incx < 1 {
		norm = 0
	} else if n == 1 {
		norm = absf32((*x)[0])
	} else {
		scale = 0
		ssq = 1
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
		//
		for ix := 0; ix < 1+(n-1)*incx; ix += incx {
			if (*x)[ix] != 0 {
				absxi = absf32((*x)[ix])
				if scale < absxi {
					ssq = 1 + ssq*powf32(scale/absxi, 2)
					scale = absxi
				} else {
					ssq += powf32(absxi/scale, 2)
				}
			}
			//Label10:
		}
		norm = scale * sqrtf32(ssq)
	}
	//
	return norm
	//
	//     End of SNRM2.
	//
}
