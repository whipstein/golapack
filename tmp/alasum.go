package goblas

import 

// Alasum prints a summary of results from one of the -CHK- routines.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Alasum( type_, nout, nfail, nrun, nerrs)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        type_
//       intEGER            nfail, nout, nrun, nerrs
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Alasum prints a summary of results from one of the -CHK- routines.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] type_
// \verbatim
//          type_ is CHARACTER*3
//          The lapACK path name.
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number on which results are to be printed.
//          nout >= 0.
// \endverbatim
//
// \param[in] nfail
// \verbatim
//          nfail is intEGER
//          The number of tests which did not pass the threshold ratio.
// \endverbatim
//
// \param[in] nrun
// \verbatim
//          nrun is intEGER
//          The total number of tests.
// \endverbatim
//
// \param[in] nerrs
// \verbatim
//          nerrs is intEGER
//          The number of error messages recorded.
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
func Alasum(_type *[]byte, nout *int, nfail *int, nrun *int, nerrs *int) {
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
	//     .. Executable Statements ..
	//
	if (*(nfail)) > 0 {
		WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %3s: %6d out of %6d tests failed to pass the threshold\n"); return &y }(), (*(_type)), (*(nfail)), (*(nrun)))
	} else {
		WRITE((*(nout)), *func() *[]byte {
			y :=[]byte("\n All tests for %3s routines passed the threshold ( %6d tests run)\n")
			return &y
		}(), (*(_type)), (*(nrun)))
	}
	if (*(nerrs)) > 0 {
		WRITE((*(nout)), *func() *[]byte {y :=[]byte("      %6d error messages recorded\n"); return &y }(), (*(nerrs)))
	}
	//
	return
	//
	//     End of Alasum
	//
}