package golapack

// Scasum ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION SCASUM(N,CX,INCX)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,N
//       ..
//       .. Array Arguments ..
//       COMPLEX CX(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    SCASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
//    returns a single precision result.
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
// \param[in,out] CX
// \verbatim
//          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of SX
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
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Scasum(n *int, cx *[]complex64, cxoff, incx *int) (scasumReturn float32) {
	var stemp float32
	var i, nincx int

	scasumReturn = 0.0
	stemp = 0.0
	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		for i = 1; i <= (*n); i++ {
			stemp = stemp + absf32(real((*cx)[i-1+(*cxoff)])) + absf32(imag((*cx)[i-1+(*cxoff)]))
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += (*incx) {
			stemp = stemp + absf32(real((*cx)[i-1+(*cxoff)])) + absf32(imag((*cx)[i-1+(*cxoff)]))
		}
	}
	scasumReturn = stemp
	return
}
