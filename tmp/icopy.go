package goblas

import 

// icopy copies an integer vector x to an integer vector y.
// Uses unrolled loops for increments equal to 1.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Icopy( n, sx, incx, sy, incy)
//
//       .. Scalar Arguments ..
//       intEGER            incx, incy, N
//       ..
//       .. Array Arguments ..
//       intEGER            sx(*), sy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// icopy copies an integer vector x to an integer vector y.
// Uses unrolled loops for increments equal to 1.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The length of the vectors sx and sy.
// \endverbatim
//
// \param[in] sx
// \verbatim
//          sx is intEGER array, dimension (1+(N-1)*abs(incx))
//          The vector X.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is intEGER
//          The spacing between consecutive elements of sx.
// \endverbatim
//
// \param[out] sy
// \verbatim
//          sy is intEGER array, dimension (1+(N-1)*abs(incy))
//          The vector Y.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is intEGER
//          The spacing between consecutive elements of sy.
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
// \date December 2016
//
// \ingroup aux_lin
//
//  =====================================================================
func Icopy(n *int, sx *[]int, incx *int, sy *[]int, incy *int) {
	i := new(int)
	ix := new(int)
	iy := new(int)
	m := new(int)
	mp1 := new(int)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*n) <= 0 {
		return
	}
	if (*(incx)) == 1 && (*(incy)) == 1 {
		goto Label20
	}
	//
	//     Code for unequal increments or equal increments not equal to 1
	//
	(*ix) = 1
	(*iy) = 1
	if (*(incx)) < 0 {
		(*ix) = (-(*n)+1)*(*(incx)) + 1
	}
	if (*(incy)) < 0 {
		(*iy) = (-(*n)+1)*(*(incy)) + 1
	}
	for (*i) = 1; (*i) <= (*n); (*i)++ {
		(*(sy))[(*iy)-(1)] = (*(sx))[(*ix)-(1)]
		(*ix) = (*ix) + (*(incx))
		(*iy) = (*iy) + (*(incy))
		//Label10:
	}
	return
	//
	//     Code for both increments equal to 1
	//
	//     Clean-up loop
	//
Label20:
	;
	(*m) = (MOD((*n), int(7)))
	if (*m) == 0 {
		goto Label40
	}
	for (*i) = 1; (*i) <= (*m); (*i)++ {
		(*(sy))[(*i)-(1)] = (*(sx))[(*i)-(1)]
		//Label30:
	}
	if (*n) < 7 {
		return
	}
Label40:
	;
	(*mp1) = (*m) + 1
	for (*i) = (*mp1); (*i) <= (*n); (*i) += 7 {
		(*(sy))[(*i)-(1)] = (*(sx))[(*i)-(1)]
		(*(sy))[(*i)+0] = (*(sx))[(*i)+0]
		(*(sy))[(*i)+1] = (*(sx))[(*i)+1]
		(*(sy))[(*i)+2] = (*(sx))[(*i)+2]
		(*(sy))[(*i)+3] = (*(sx))[(*i)+3]
		(*(sy))[(*i)+4] = (*(sx))[(*i)+4]
		(*(sy))[(*i)+5] = (*(sx))[(*i)+5]
		//Label50:
	}
	return
	//
	//     End of icopy
	//
}
