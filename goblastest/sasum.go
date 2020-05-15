package main

// Sasum ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION SASUM(N,SX,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       REAL SX(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    SASUM takes the sum of the absolute values.
//    uses unrolled loops for increment equal to one.
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
// \param[in] SX
// \verbatim
//          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
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
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Sasum(n int, sx *[]float32, incx int) float32 {
	var stemp float32
	var m, mp1, nincx int
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
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	stemp = 0.0
	if n <= 0 || incx <= 0 {
		return 0
	}
	if incx == 1 {
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = mod(n, 6)
		if m != 0 {
			for i := 0; i < m; i++ {
				stemp += absf32((*sx)[i])
			}
			if n < 6 {
				return stemp
			}
		}
		mp1 = m + 1
		for i := mp1 - 1; i < n; i += 6 {
			stemp += absf32((*sx)[i]) + absf32((*sx)[i+1]) + absf32((*sx)[i+2]) + absf32((*sx)[i+3]) + absf32((*sx)[i+4]) + absf32((*sx)[i+5])
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = n * incx
		for i := 0; i < nincx; i += incx {
			stemp += absf32((*sx)[i])
		}
	}
	return stemp
}
