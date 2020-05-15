package goblas

// Saxpy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SAXPY(n,sa,sx,incx,sy,incy)
//
//       .. Scalar Arguments ..
//       REAL sa
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
//    SAXPY constant times a vector plus a vector.
//    uses unrolled loops for increments equal to one.
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
// \param[in] sa
// \verbatim
//          sa is REAL
//           On entry, sa specifies the scalar alpha.
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
func Saxpy(major *byte, n *int, sa *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	var i, ix, iy, m, mp1 int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 {
		return
	}
	if *sa == 0.0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//        code for both increments equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 4
		if m != 0 {
			for i = 1; i <= m; i++ {
				(*sy)[i-1] += (*sa) * (*sx)[i-1]
			}
		}
		if *n < 4 {
			return
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 4 {
			(*sy)[i-1] += (*sa) * (*sx)[i-1]
			(*sy)[i] += (*sa) * (*sx)[i]
			(*sy)[i+1] += (*sa) * (*sx)[i+1]
			(*sy)[i+2] += (*sa) * (*sx)[i+2]
		}
	} else {
		//
		//        code for unequal increments or equal increments
		//          not equal to 1
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
			(*sy)[iy-1] += (*sa) * (*sx)[ix-1]
			ix += *incx
			iy += *incy
		}
	}
	return
}
