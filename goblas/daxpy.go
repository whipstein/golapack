package goblas

// Daxpy ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Daxpy(n,da,dx,incx,dy,incy)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION da
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
//    Daxpy constant times a vector plus a vector.
//    uses unrolled loops for increments equal to one.
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
// \param[in] da
// \verbatim
//          da is DOUBLE PRECISION
//           On entry, da specifies the scalar alpha.
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
// \param[in,out] dy
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
func Daxpy(major *byte, n *int, da *float64, dx *[]float64, incx *int, dy *[]float64, incy *int) {
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
	if *da == 0.0e0 {
		return
	}
	if *incx == 1 && *incy == 1 {
		//
		//        code for both increments equal to 1
		//
		//
		//        clean-up loop
		//
		m = (*n) % 4
		if m != 0 {
			for i = 0; i < m; i++ {
				(*dy)[i] += (*da) * (*dx)[i]
			}
		}
		if *n < 4 {
			return
		}
		mp1 = m + 1
		for i = mp1; i <= *n; i += 4 {
			(*dy)[i-1] += (*da) * (*dx)[i-1]
			(*dy)[i] += (*da) * (*dx)[i]
			(*dy)[i+1] += (*da) * (*dx)[i+1]
			(*dy)[i+2] += (*da) * (*dx)[i+2]
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
			(*dy)[iy-1] += (*da) * (*dx)[ix-1]
			ix += (*incx)
			iy += (*incy)
		}
	}
	return
}
