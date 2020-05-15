package goblas

import 

// Alaesm prints a summary of results from one of the -err- routines.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Alaesm( path, ok, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            ok
//       CHARACTER*3        path
//       intEGER            nout
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Alaesm prints a summary of results from one of the -err- routines.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The lapACK path name.
// \endverbatim
//
// \param[in] ok
// \verbatim
//          ok is LOGICAL
//          The flag from CHKXER that indicates whether or not the tests
//          of error exits passed.
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number on which results are to be printed.
//          nout >= 0.
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
func Alaesm(path *[]byte, ok *bool, nout *int) {
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
	if *(ok) {
		WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %3s routines passed the tests of the error exits\n"); return &y}(), (*(path)))
	} else {
		WRITE((*(nout)), *func() *[]byte {y :=[]byte(" *** %3s routines failed the tests of the error exits ***\n"); return &y}(), (*(path)))
	}
	//
	return
	//
	//     End of Alaesm
	//
}
