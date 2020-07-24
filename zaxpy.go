package golapack

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
//       SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 ZA
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 ZX(*),ZY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    ZAXPY constant times a vector plus a vector.
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
// \param[in] ZA
// \verbatim
//          ZA is COMPLEX*16
//           On entry, ZA specifies the scalar alpha.
// \endverbatim
//
// \param[in] ZX
// \verbatim
//          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of ZX
// \endverbatim
//
// \param[in,out] ZY
// \verbatim
//          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of ZY
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
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Zaxpy(n *int, za *complex128, zx *[]complex128, zxoff, incx *int, zy *[]complex128, zyoff, incy *int) {
	var i, ix, iy int

	if (*n) <= 0 {
		return
	}
	if abssumc128(*za) == 0.0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//        code for both increments equal to 1
		//
		for i = 1; i <= (*n); i++ {
			(*zy)[i-1+(*zyoff)] = (*zy)[i-1+(*zyoff)] + (*za)*(*zx)[i-1+(*zxoff)]
		}
	} else {
		//
		//        code for unequal increments or equal increments
		//          not equal to 1
		//
		ix = 1
		iy = 1
		if (*incx) < 0 {
			ix = (-(*n)+1)*(*incx) + 1
		}
		if (*incy) < 0 {
			iy = (-(*n)+1)*(*incy) + 1
		}
		for i = 1; i <= (*n); i++ {
			(*zy)[iy-1+(*zyoff)] = (*zy)[iy-1+(*zyoff)] + (*za)*(*zx)[ix-1+(*zxoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}

}
