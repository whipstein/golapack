package golapack

// Dswap ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DX(*),DY(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    DSWAP interchanges two vectors.
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
// \param[in,out] DX
// \verbatim
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of DX
// \endverbatim
//
// \param[in,out] DY
// \verbatim
//          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of DY
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
// \ingroup double_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Dswap(n *int, dx *[]float64, dxoff, incx *int, dy *[]float64, dyoff, incy *int) {
	var dtemp float64
	var i, ix, iy, m, mp1 int

	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//       code for both increments equal to 1
		//
		//
		//       clean-up loop
		//
		m = modint(*n, int(3))
		if m != 0 {
			for i = 1; i <= m; i++ {
				dtemp = (*dx)[i-1+(*dxoff)]
				(*dx)[i-1+(*dxoff)] = (*dy)[i-1+(*dyoff)]
				(*dy)[i-1+(*dyoff)] = dtemp
			}
			if (*n) < 3 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 3 {
			dtemp = (*dx)[i-1+(*dxoff)]
			(*dx)[i-1+(*dxoff)] = (*dy)[i-1+(*dyoff)]
			(*dy)[i-1+(*dyoff)] = dtemp
			dtemp = (*dx)[i+1-1+(*dxoff)]
			(*dx)[i+1-1+(*dxoff)] = (*dy)[i+1-1+(*dyoff)]
			(*dy)[i+1-1+(*dyoff)] = dtemp
			dtemp = (*dx)[i+2-1+(*dxoff)]
			(*dx)[i+2-1+(*dxoff)] = (*dy)[i+2-1+(*dyoff)]
			(*dy)[i+2-1+(*dyoff)] = dtemp
		}
	} else {
		//
		//       code for unequal increments or equal increments not equal
		//         to 1
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
			dtemp = (*dx)[ix-1+(*dxoff)]
			(*dx)[ix-1+(*dxoff)] = (*dy)[iy-1+(*dyoff)]
			(*dy)[iy-1+(*dyoff)] = dtemp
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
