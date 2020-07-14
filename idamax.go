package golapack

// Idamax ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION IDAMAX(N,DX,INCX)
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
//    IDAMAX finds the index of the first element having maximum absolute value.
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
// \ingroup aux_blas
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
func Idamax(n *int, dx *[]float64, dxoff, incx *int) (idamaxReturn int) {
	var dmax float64
	var i, ix int

	idamaxReturn = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	idamaxReturn = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		dmax = absf64((*dx)[0+(*dxoff)])
		for i = 2; i <= (*n); i++ {
			if absf64((*dx)[i-1+(*dxoff)]) > dmax {
				idamaxReturn = i
				dmax = absf64((*dx)[i-1+(*dxoff)])
			}
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		ix = 1
		dmax = absf64((*dx)[0+(*dxoff)])
		ix = ix + (*incx)
		for i = 2; i <= (*n); i++ {
			if absf64((*dx)[ix-1+(*dxoff)]) > dmax {
				idamaxReturn = i
				dmax = absf64((*dx)[ix-1+(*dxoff)])
			}
			ix = ix + (*incx)
		}
	}
	return
}
