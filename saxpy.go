package golapack

// Saxpy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
//
//       .. Scalar Arguments ..
//       REAL SA
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
//    SAXPY constant times a vector plus a vector.
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
// \param[in] SA
// \verbatim
//          SA is REAL
//           On entry, SA specifies the scalar alpha.
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
func Saxpy(n *int, sa *float32, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int) {
	var i, ix, iy, m, mp1 int

	if (*n) <= 0 {
		return
	}
	if (*sa) == 0.0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//
		//        code for both increments equal to 1
		//
		//
		//        clean-up loop
		//
		m = modint(*n, int(4))
		if m != 0 {
			for i = 1; i <= m; i++ {
				(*sy)[i-1+(*syoff)] = (*sy)[i-1+(*syoff)] + (*sa)*(*sx)[i-1+(*sxoff)]
			}
		}
		if (*n) < 4 {
			return
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 4 {
			(*sy)[i-1+(*syoff)] = (*sy)[i-1+(*syoff)] + (*sa)*(*sx)[i-1+(*sxoff)]
			(*sy)[i+1-1+(*syoff)] = (*sy)[i+1-1+(*syoff)] + (*sa)*(*sx)[i+1-1+(*sxoff)]
			(*sy)[i+2-1+(*syoff)] = (*sy)[i+2-1+(*syoff)] + (*sa)*(*sx)[i+2-1+(*sxoff)]
			(*sy)[i+3-1+(*syoff)] = (*sy)[i+3-1+(*syoff)] + (*sa)*(*sx)[i+3-1+(*sxoff)]
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
			(*sy)[iy-1+(*syoff)] = (*sy)[iy-1+(*syoff)] + (*sa)*(*sx)[ix-1+(*sxoff)]
			ix = ix + (*incx)
			iy = iy + (*incy)
		}
	}
}
