package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dlaord sorts the elements of a vector x in increasing or decreasing
// order.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlaord( job, n, x, incx)
//
//       .. Scalar Arguments ..
//       CHARACTER          job
//       intEGER            incx, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   X(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlaord sorts the elements of a vector x in increasing or decreasing
// order.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] job
// \verbatim
//          job is CHARACTER
//          = 'I':  Sort in increasing order
//          = 'D':  Sort in decreasing order
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The length of the vector X.
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension
//                         (1+(N-1)*incX)
//          On entry, the vector of length n to be sorted.
//          On exit, the vector x is sorted in the prescribed order.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is intEGER
//          The spacing between successive elements of X.  incx >= 0.
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
// \ingroup double_lin
//
//  =====================================================================
func Dlaord(job *byte, n *int, x *[]float64, incx *int) {
	i := new(int)
	inc := new(int)
	ix := new(int)
	ixnext := new(int)
	temp := new(float64)
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
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	(*inc) = (ABS((*(incx))))
	if blas.Lsame((job), func() *byte {y := byte('I'); return &y }()) {
		//
		//        Sort in increasing order
		//
		for (*i) = 2; (*i) <= (*(n)); (*i)++ {
			(*ix) = 1 + ((*i)-1)*(*inc)
		Label10:
			;
			if (*ix) == 1 {
				goto Label20
			}
			(*ixnext) = (*ix) - (*inc)
			if (*(x))[(*ix)-1] > (*(x))[(*ixnext)-1] {
				goto Label20
			} else {
				(*temp) = (*(x))[(*ix)-1]
				(*(x))[(*ix)-1] = (*(x))[(*ixnext)-1]
				(*(x))[(*ixnext)-1] = (*temp)
			}
			(*ix) = (*ixnext)
			goto Label10
		Label20:
		}
		//
	} else if blas.Lsame((job), func() *byte {y := byte('D'); return &y }()) {
		//
		//        Sort in decreasing order
		//
		for (*i) = 2; (*i) <= (*(n)); (*i)++ {
			(*ix) = 1 + ((*i)-1)*(*inc)
		Label30:
			;
			if (*ix) == 1 {
				goto Label40
			}
			(*ixnext) = (*ix) - (*inc)
			if (*(x))[(*ix)-1] < (*(x))[(*ixnext)-1] {
				goto Label40
			} else {
				(*temp) = (*(x))[(*ix)-1]
				(*(x))[(*ix)-1] = (*(x))[(*ixnext)-1]
				(*(x))[(*ixnext)-1] = (*temp)
			}
			(*ix) = (*ixnext)
			goto Label30
		Label40:
		}
	}
	return
	//
	//     End of Dlaord
	//
}
