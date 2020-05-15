package goblas

// Dcopy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Dcopy(n,dx,incx,dy,incy)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION dx(*),dy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Dcopy copies a vector, x, to a vector, y.
//    uses unrolled loops for increments equal to 1.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] dx
// \verbatim
//          dx is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of dx
// \endverbatim
//
// \param[out] dy
// \verbatim
//          dy is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of dy
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
//     modified 12/3/93, array1 declarations changed to array(*)
// \endverbatim
//
//  =====================================================================
func Dcopy(major *byte, n *int, dx *[]float64, incx *int, dy *[]float64, incy *int) {
	var i, ix, iy, m, mp1 int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	if *n <= 0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//        code for both increments equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 7
		if m != 0 {
			for i = 0; i < m; i++ {
				(*dy)[i] = (*dx)[i]
			}
			if *n < 7 {
				return
			}
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 7 {
			(*dy)[i-1] = (*dx)[i-1]
			(*dy)[i] = (*dx)[i]
			(*dy)[i+1] = (*dx)[i+1]
			(*dy)[i+2] = (*dx)[i+2]
			(*dy)[i+3] = (*dx)[i+3]
			(*dy)[i+4] = (*dx)[i+4]
			(*dy)[i+5] = (*dx)[i+5]
		}
	} else {
		//
		//        code for unequal increments or equal increments
		//          not equal to 1
		//
		ix = 1
		iy = 1
		if *incx < 0 {
			ix = (-(*n)+1)*(*incx) + 1
		}
		if *incy < 0 {
			iy = (-(*n)+1)*(*incy) + 1
		}
		for i = 1; i <= *n; i++ {
			(*dy)[iy-1] = (*dx)[ix-1]
			ix += *incx
			iy += *incy
		}
	}
	return
}
