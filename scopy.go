package golapack

// Scopy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
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
//    SCOPY copies a vector, x, to a vector, y.
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
// \param[out] SY
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
func Scopy(n *int, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int) {
	var i, ix, iy, m, mp1 int

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
		m = modint(*n, int(7))
		if m != 0 {
			for i = 1; i <= m; i++ {
				(*sy)[i-1+(*syoff)] = (*sx)[i-1+(*sxoff)]
			}
			if (*n) < 7 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 7 {
			(*sy)[i-1+(*syoff)] = (*sx)[i-1+(*sxoff)]
			(*sy)[i+(*syoff)] = (*sx)[i+(*sxoff)]
			(*sy)[i+1+(*syoff)] = (*sx)[i+1+(*sxoff)]
			(*sy)[i+2+(*syoff)] = (*sx)[i+2+(*sxoff)]
			(*sy)[i+3+(*syoff)] = (*sx)[i+3+(*sxoff)]
			(*sy)[i+4+(*syoff)] = (*sx)[i+4+(*sxoff)]
			(*sy)[i+5+(*syoff)] = (*sx)[i+5+(*sxoff)]
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
			(*sy)[iy-1+(*syoff)] = (*sx)[ix-1+(*sxoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
