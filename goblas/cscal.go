package goblas

// Cscal ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cscal(n,ca,cx,incx)
//
//       .. Scalar Arguments ..
//       COMPLEX ca
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX cx(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Cscal scales a vector by a constant.
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
// \param[in] ca
// \verbatim
//          ca is COMPLEX
//           On entry, ca specifies the scalar alpha.
// \endverbatim
//
// \param[in,out] cx
// \verbatim
//          cx is COMPLEX array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of cx
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
// \ingroup complex_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack,  3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Cscal(n *int, ca *complex64, cx *[]complex64, incx *int) {
	var i, nincx int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
	//     November 2017
	//
	if *n <= 0 || *incx <= 0 {
		return
	}
	if *incx == 1 {
		//
		//        code for increment equal to 1
		//
		for i = 1; i <= *n; i++ {
			(*cx)[i-1] = (*ca) * (*cx)[i-1]
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += *incx {
			(*cx)[i-1] = (*ca) * (*cx)[i-1]
		}
	}
	return
}
