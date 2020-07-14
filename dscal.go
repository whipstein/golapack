package golapack

// Dscal ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DSCAL(N,DA,DX,INCX)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION DA
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
//    DSCAL scales a vector by a constant.
//    uses unrolled loops for increment equal to 1.
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
// \param[in] DA
// \verbatim
//          DA is DOUBLE PRECISION
//           On entry, DA specifies the scalar alpha.
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
func Dscal(n *int, da *float64, dx *[]float64, dxoff, incx *int) {
	var i, m, mp1, nincx int

	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		//
		//        clean-up loop
		//
		m = modint(*n, int(5))
		if m != 0 {
			for i = 1; i <= m; i++ {
				(*dx)[i-1+(*dxoff)] = (*da) * (*dx)[i-1+(*dxoff)]
			}
			if (*n) < 5 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= (*n); i += 5 {
			(*dx)[i-1+(*dxoff)] = (*da) * (*dx)[i-1+(*dxoff)]
			(*dx)[i+1-1+(*dxoff)] = (*da) * (*dx)[i+1-1+(*dxoff)]
			(*dx)[i+2-1+(*dxoff)] = (*da) * (*dx)[i+2-1+(*dxoff)]
			(*dx)[i+3-1+(*dxoff)] = (*da) * (*dx)[i+3-1+(*dxoff)]
			(*dx)[i+4-1+(*dxoff)] = (*da) * (*dx)[i+4-1+(*dxoff)]
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += (*incx) {
			(*dx)[i-1+(*dxoff)] = (*da) * (*dx)[i-1+(*dxoff)]
		}
	}
}
