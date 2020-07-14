package golapack

// Sswap ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
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
//    SSWAP interchanges two vectors.
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
// \param[in,out] SX
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
// \param[in,out] SY
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
func Sswap(n *int, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int) {
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
				(*sx)[i-1+(*sxoff)], (*sy)[i-1+(*syoff)] = (*sy)[i-1+(*syoff)], (*sx)[i-1+(*sxoff)]
			}
			if (*n) < 3 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 3 {
			(*sx)[i-1+(*sxoff)], (*sy)[i-1+(*syoff)] = (*sy)[i-1+(*syoff)], (*sx)[i-1+(*sxoff)]
			(*sx)[i+(*sxoff)], (*sy)[i+(*syoff)] = (*sy)[i+(*syoff)], (*sx)[i+(*sxoff)]
			(*sx)[i+1+(*sxoff)], (*sy)[i+1+(*syoff)] = (*sy)[i+1+(*syoff)], (*sx)[i+1+(*sxoff)]
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
			(*sx)[ix-1+(*sxoff)], (*sy)[iy-1+(*syoff)] = (*sy)[iy-1+(*syoff)], (*sx)[ix-1+(*sxoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
