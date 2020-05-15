package goblas

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
//       REAL FUNCTION SNRM2(n,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       REAL x(*)
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
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] x
// \verbatim
//          x is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
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
func Snrm2(major *byte, n *int, x *[]float32, incx *int) float32 {
	var absxi, norm, scale, ssq float32
	var ix int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n < 1 || *incx < 1 {
		norm = 0.0
	} else if *n == 1 {
		norm = absf32((*x)[0])
	} else {
		scale = 0.0
		ssq = 1.0
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL SLASSQ( n, x, incx, scale, ssq )
		//
		for ix = 1; ix <= 1+((*n)-1)*(*incx); ix += *incx {
			if (*x)[ix-1] != 0.0 {
				absxi = absf32((*x)[ix-1])
				if scale < absxi {
					ssq = 1.0 + ssq*powf32((scale/absxi), 2)
					scale = absxi
				} else {
					ssq = ssq + powf32((absxi/scale), 2)
				}
			}
		}
		norm = scale * sqrtf32(ssq)
	}

	return norm
	//
	//     End of SNRM2.
	//
}
