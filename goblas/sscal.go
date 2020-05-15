package goblas

// Sscal ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSCAL(n,sa,sx,incx)
//
//       .. Scalar Arguments ..
//       REAL sa
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
//    SSCAL scales a vector by a constant.
//    uses unrolled loops for increment equal to 1.
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
func Sscal(major *byte, n *int, sa *float32, sx *[]float32, incx *int) {
	var i, m, mp1, nincx int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 || *incx <= 0 {
		return
	}
	if *incx == 1 {
		//
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 5
		if m != 0 {
			for i = 1; i <= m; i++ {
				(*sx)[i-1] *= *sa
			}
			if *n < 5 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 5 {
			(*sx)[i-1] *= *sa
			(*sx)[i] *= *sa
			(*sx)[i+1] *= *sa
			(*sx)[i+2] *= *sa
			(*sx)[i+3] *= *sa
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += *incx {
			(*sx)[i-1] *= *sa
		}
	}
	return
}
