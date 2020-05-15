package goblas

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
//       REAL FUNCTION SASUM(n,sx,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       REAL sx(*)
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
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] sx
// \verbatim
//          sx is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of sx
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
func Sasum(major *byte, n *int, sx *[]float32, incx *int) float32 {
	var stemp float32
	var i, m, mp1, nincx int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	stemp = 0.0
	if *n <= 0 || *incx <= 0 {
		return stemp
	}
	if *incx == 1 {
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 6
		if m != 0 {
			for i = 1; i <= m; i++ {
				stemp += absf32((*sx)[i-1])
			}
			if *n < 6 {
				return stemp
			}
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 6 {
			stemp += absf32((*sx)[i-1]) + absf32((*sx)[i]) + absf32((*sx)[i+1]) + absf32((*sx)[i+2]) + absf32((*sx)[i+3]) + absf32((*sx)[i+4])
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += *incx {
			stemp += absf32((*sx)[i-1])
		}
	}
	return stemp
}
