package main

// Sswap ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       REAL SX(*),SY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    SSWAP interchanges two vectors.
//    uses unrolled loops for increments equal to 1.
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
// \param[in,out] SX
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
// \param[in,out] SY
// \verbatim
//          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of SY
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
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Sswap(n int, sx *[]float32, incx int, sy *[]float32, incy int) {
	var stemp float32
	var ix, iy, m, mp1 int
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
	if n <= 0 {
		return
	}
	if incx == 1 && incy == 1 {
		//
		//       code for both increments equal to 1
		//
		//
		//       clean-up loop
		//
		m = mod(n, 3)
		if m != 0 {
			for i := 0; i < m; i++ {
				stemp = (*sx)[i]
				(*sx)[i] = (*sy)[i]
				(*sy)[i] = stemp
			}
			if n < 3 {
				return
			}
		}
		mp1 = m + 1
		for i := mp1 - 1; i < n; i += 3 {
			stemp = (*sx)[i]
			(*sx)[i] = (*sy)[i]
			(*sy)[i] = stemp
			stemp = (*sx)[i+1]
			(*sx)[i+1] = (*sy)[i+1]
			(*sy)[i+1] = stemp
			stemp = (*sx)[i+2]
			(*sx)[i+2] = (*sy)[i+2]
			(*sy)[i+2] = stemp
		}
	} else {
		//
		//       code for unequal increments or equal increments not equal
		//         to 1
		//
		ix = 1
		iy = 1
		if incx < 0 {
			ix = (-n+1)*incx + 1
		}
		if incy < 0 {
			iy = (-n+1)*incy + 1
		}
		for i := 0; i < n; i++ {
			stemp = (*sx)[ix-1]
			(*sx)[ix-1] = (*sy)[iy-1]
			(*sy)[iy-1] = stemp
			ix += incx
			iy += incy
		}
	}
	return
}
