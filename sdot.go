package golapack

// Sdot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
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
//    SDOT forms the dot product of two vectors.
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
// \param[in] SX
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
// \param[in] SY
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
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Sdot(n *int, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int) (sdotReturn float32) {
	var stemp float32
	var i, ix, iy, m, mp1 int

	stemp = 0.0
	sdotReturn = 0.0
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
				stemp = stemp + (*sx)[i-1+(*sxoff)]*(*sy)[i-1+(*syoff)]
			}
			if (*n) < 5 {
				sdotReturn = stemp
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 5 {
			stemp = stemp + (*sx)[i-1+(*sxoff)]*(*sy)[i-1+(*syoff)] + (*sx)[i+1-1+(*sxoff)]*(*sy)[i+1-1+(*syoff)] + (*sx)[i+2-1+(*sxoff)]*(*sy)[i+2-1+(*syoff)] + (*sx)[i+3-1+(*sxoff)]*(*sy)[i+3-1+(*syoff)] + (*sx)[i+4-1+(*sxoff)]*(*sy)[i+4-1+(*syoff)]
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
			stemp = stemp + (*sx)[ix-1+(*sxoff)]*(*sy)[iy-1+(*syoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
	sdotReturn = stemp
	return
}
