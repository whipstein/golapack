package goblas

import 

// Chkxer ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Chkxer( srnamt, infot, nout, lerr, ok)
//
//       .. Scalar Arguments ..
//       LOGICAL            lerr, ok
//       CHARACTER*(*)      srnamt
//       intEGER            infot, nout
//
//
// \par Purpose:
//  =============
//
// \verbatim
// \endverbatim
//
//  Arguments:
//  ==========
//
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG Ltd.
//
// \date June 2017
//
// \ingroup complex_lin
//
//  =====================================================================
func Chkxer(srnamt *[]byte, infot *int, nout *int, lerr *bool, ok *bool) {
	//
	//  -- lapACK test routine (version 3.7.1) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2017
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	if !(*(lerr)) {
		WRITE((*(nout)), *func() *[]byte {
			y :=[]byte(" *** illegal value of parameter number %2d not detected by %6s ***\n")
			return &y
		}(), (*(infot)), (*(srnamt))[1:lenTrim((srnamt))-1])
		(*(ok)) = false
	}
	(*(lerr)) = false
	return
	//
	//
	//     End of Chkxer.
	//
}
