package goblas

import 

// Dget06 computes a test ratio to compare two values for rcond.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dget06( rcond, rcondc)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   rcond, rcondc
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dget06 computes a test ratio to compare two values for rcond.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The estimate of the reciprocal of the condition number of a,
//          as computed by DGECON.
// \endverbatim
//
// \param[in] rcondc
// \verbatim
//          rcondc is DOUBLE PRECISION
//          The reciprocal of the condition number of a, computed as
//          ( 1/norm(a)) / norm(inv(a)).
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
func Dget06(rcond *float64, rcondc *float64) (dget06Return *float64) {
	dget06Return = new(float64)
	zero := new(float64)
	one := new(float64)
	eps := new(float64)
	rat := new(float64)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	if (*(rcond)) > (*zero) {
		if (*(rcondc)) > (*zero) {
			(*rat) = MAX((*(rcond)), (*(rcondc)))/Min((*(rcond)), (*(rcondc))) - ((*one) - (*eps))
		} else {
			(*rat) = (*(rcond)) / (*eps)
		}
	} else {
		if (*(rcondc)) > (*zero) {
			(*rat) = (*(rcondc)) / (*eps)
		} else {
			(*rat) = (*zero)
		}
	}
	(*(dget06Return)) = (*rat)
	return
	//
	//     End of Dget06
	//
}
