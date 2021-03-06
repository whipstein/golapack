package golapack

// Cscal ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE CSCAL(N,CA,CX,INCX)
//
//       .. Scalar Arguments ..
//       COMPLEX CA
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
//    CSCAL scales a vector by a constant.
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
// \param[in] CA
// \verbatim
//          CA is COMPLEX
//           On entry, CA specifies the scalar alpha.
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
//         storage spacing between elements of CX
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
// \ingroup complex_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack,  3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Cscal(n *int, ca *complex64, cx *[]complex64, cxoff, incx *int) {
	var i, nincx int

	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//
		//        code for increment equal to 1
		//
		for i = 1; i <= (*n); i++ {
			(*cx)[i-1+(*cxoff)] = (*ca) * (*cx)[i-1+(*cxoff)]
		}
	} else {
		//
		//        code for increment not equal to 1
		//
		nincx = (*n) * (*incx)
		for i = 1; i <= nincx; i += (*incx) {
			(*cx)[i-1+(*cxoff)] = (*ca) * (*cx)[i-1+(*cxoff)]
		}
	}
}
