package goblas

import 

// Dgennd tests that its argument has a non-negative diagonal.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       LOGICAL FUNCTION Dgennd (m, n, a, lda)
//
//       .. Scalar Arguments ..
//       inTEGER m, n, lda
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION a( lda, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Dgennd tests that its argument has a non-negative diagonal.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows in A.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns in A.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda, N)
//          The matrix.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          Leading dimension of A.
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
func Dgennd(m *int, n *int, a *[][]float64, lda *int) (dgenndReturn *bool) {
	dgenndReturn = new(bool)
	zero := new(float64)
	i := new(int)
	k := new(int)
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
	//     .. Parameters ..
	(*zero) = 0.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsics ..
	//     ..
	//     .. Executable Statements ..
	(*k) = (Min((*(m)), (*(n))))
	for (*i) = 1; (*i) <= (*k); (*i)++ {
		if (*(a))[(*i)-1][(*i)-1] < (*zero) {
			(*(dgenndReturn)) = false
			return
		}
	}
	(*(dgenndReturn)) = true
	return
}
