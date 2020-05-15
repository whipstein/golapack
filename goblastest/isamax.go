package main

// Isamax ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION ISAMAX(N,SX,INCX)
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
//    ISAMAX finds the index of the first element having maximum absolute value.
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
// \ingroup aux_blas
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
func Isamax(n int, sx *[]float32, incx int) int {
	var smax float32
	var isamaxReturn, ix int
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
	if n < 1 || incx <= 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	if incx == 1 {
		//
		//        code for increment equal to 1
		//
		smax = absf32((*sx)[0])
		for i := 1; i < n; i++ {
			if absf32((*sx)[i]) > smax {
				isamaxReturn = i + 1
				smax = absf32((*sx)[i])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		smax = absf32((*sx)[0])
		ix += incx
		for i := 1; i < n; i++ {
			if absf32((*sx)[ix-1]) > smax {
				isamaxReturn = i + 1
				smax = absf32((*sx)[ix-1])
			}
			ix += incx
		}
	}
	return isamaxReturn
}
