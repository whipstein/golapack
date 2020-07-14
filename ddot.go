package golapack

// Ddot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
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
//    DDOT forms the dot product of two vectors.
//    uses unrolled loops for increments equal to one.
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
// \param[in] DX
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
// \param[in] DY
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
func Ddot(n *int, dx *[]float64, dxoff, incx *int, dy *[]float64, dyoff, incy *int) (ddotReturn float64) {
	var dtemp float64
	var i, ix, iy, m, mp1 int

	ddotReturn = 0.0
	dtemp = 0.0
	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//        code for both increments equal to 1
		//
		//
		//        clean-up loop
		//
		m = modint(*n, int(5))
		if m != 0 {
			for i = 1; i <= m; i++ {
				dtemp = dtemp + (*dx)[i-1+(*dxoff)]*(*dy)[i-1+(*dyoff)]
			}
			if (*n) < 5 {
				ddotReturn = dtemp
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 5 {
			dtemp = dtemp + (*dx)[i-1+(*dxoff)]*(*dy)[i-1+(*dyoff)] + (*dx)[i+1-1+(*dxoff)]*(*dy)[i+1-1+(*dyoff)] + (*dx)[i+2-1+(*dxoff)]*(*dy)[i+2-1+(*dyoff)] + (*dx)[i+3-1+(*dxoff)]*(*dy)[i+3-1+(*dyoff)] + (*dx)[i+4-1+(*dxoff)]*(*dy)[i+4-1+(*dyoff)]
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
			dtemp = dtemp + (*dx)[ix-1+(*dxoff)]*(*dy)[iy-1+(*dyoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
	ddotReturn = dtemp
	return
}
