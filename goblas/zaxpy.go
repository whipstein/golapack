package goblas

// Zaxpy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zaxpy(n,za,zx,incx,zy,incy)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 za
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 zx(*),zy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Zaxpy constant times a vector plus a vector.
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
// \param[in] za
// \verbatim
//          za is COMPLEX*16
//           On entry, za specifies the scalar alpha.
// \endverbatim
//
// \param[in] zx
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
// \param[in,out] zy
// \verbatim
//          zy is COMPLEX*16 array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of zy
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
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Zaxpy(n *int, za *complex128, zx *[]complex128, incx *int, zy *[]complex128, incy *int) {
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
	if dcabs1(za) == 0.0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//        code for both increments equal to 1
		//
		for i = 1; i <= *n; i++ {
			(*zy)[i-1] += (*za) * (*zx)[i-1]
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
			(*zy)[iy-1] += (*za) * (*zx)[ix-1]
			ix += *incx
			iy += *incy
		}
	}
}
