package goblas

// Zdscal ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zdscal(n,da,zx,incx)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION da
//       INTEGER incx,n
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 zx(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Zdscal scales a vector by a constant.
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
// \param[in] da
// \verbatim
//          da is DOUBLE PRECISION
//           On entry, da specifies the scalar alpha.
// \endverbatim
//
// \param[in,out] zx
// \verbatim
//          zx is COMPLEX*16 array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of zx
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
// \ingroup complex16_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Zdscal(n *int, da *float64, zx *[]complex128, incx *int) {
	var i, nincx int
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
		for i = 1; i <= *n; i++ {
			(*zx)[i-1] = complex(*da, 0.0) * (*zx)[i-1]
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = *n * (*incx)
		for i = 1; i <= nincx; i += *incx {
			(*zx)[i-1] = complex(*da, 0.0) * (*zx)[i-1]
		}
	}
}
