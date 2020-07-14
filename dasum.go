package golapack

// Dasum ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DX(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    DASUM takes the sum of the absolute values.
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
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*absf64( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of DX
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
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Dasum(n *int, dx *[]float64, dxoff, incx *int) (dasumReturn float64) {
	var dtemp float64
	var i, m, mp1, nincx int

	dasumReturn = 0.0
	dtemp = 0.0
	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = modint(*n, int(6))
		if m != 0 {
			for i = 1; i <= m; i++ {
				dtemp = dtemp + absf64((*dx)[i-1+(*dxoff)])
			}
			if (*n) < 6 {
				dasumReturn = dtemp
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 6 {
			dtemp = dtemp + absf64((*dx)[i-1+(*dxoff)]) + absf64((*dx)[i+1-1+(*dxoff)]) + absf64((*dx)[i+2-1+(*dxoff)]) + absf64((*dx)[i+3-1+(*dxoff)]) + absf64((*dx)[i+4-1+(*dxoff)]) + absf64((*dx)[i+5-1+(*dxoff)])
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += (*incx) {
			dtemp = dtemp + absf64((*dx)[i-1+(*dxoff)])
		}
	}
	dasumReturn = dtemp
	return
}
