package goblas

// Srot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SROT(n,sx,incx,sy,incy,c,s)
//
//       .. Scalar Arguments ..
//       REAL c,s
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       REAL sx(*),sy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    applies a plane rotation.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vectors
// \endverbatim
//
// \param[in,out] sx
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
// \param[in,out] sy
// \verbatim
//          sy is REAL array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of sy
// \endverbatim
//
// \param[in] c
// \verbatim
//          c is REAL
// \endverbatim
//
// \param[in] s
// \verbatim
//          s is REAL
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
func Srot(major *byte, n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, c, s *float32) {
	var stemp float32
	var i, ix, iy int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//       code for both increments equal to 1
		//
		for i = 1; i <= *n; i++ {
			stemp = (*c)*(*sx)[i-1] + (*s)*(*sy)[i-1]
			(*sy)[i-1] = (*c)*(*sy)[i-1] - (*s)*(*sx)[i-1]
			(*sx)[i-1] = stemp
		}
	} else {
		//
		//       code for unequal increments or equal increments not equal
		//         to 1
		//
		ix = 1
		iy = 1
		if *incx < 0 {
			ix = (-(*n)+1)*(*incx) + 1
		}
		if *incy < 0 {
			iy = (-(*n)+1)*(*incy) + 1
		}
		for i = 1; i <= *n; i++ {
			stemp = (*c)*(*sx)[ix-1] + (*s)*(*sy)[iy-1]
			(*sy)[iy-1] = (*c)*(*sy)[iy-1] - (*s)*(*sx)[ix-1]
			(*sx)[ix-1] = stemp
			ix += *incx
			iy += *incy
		}
	}
	return
}
