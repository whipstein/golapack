package goblas

import (
	"bytes"
)

// Xerbla is a special version of Xerbla to be used only as part of
// the test program for testing error exits from the lapACK routines.
// Error messages are printed if info.NE.infot or if srname.NE.srnamt,
// where infot and srnamt are values stored in common.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Xerbla( srname, info)
//
//       .. Scalar Arguments ..
//       CHARACTER*(*)      srname
//       intEGER            info
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// This is a special version of Xerbla to be used only as part of
// the test program for testing error exits from the lapACK routines.
// Error messages are printed if info.NE.infot or if srname.NE.srnamt,
// where infot and srnamt are values stored in common.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] srname
// \verbatim
//          srname is CHARACTER*(*)
//          The name of the subroutine calling Xerbla.  This name should
//          match the common variable srnamt.
// \endverbatim
//
// \param[in] info
// \verbatim
//          info is intEGER
//          The error return code from the calling subroutine.  info
//          should equal the common variable infot.
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
// \ingroup aux_eig
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  The following variables are passed via the common blocks infoc and
//  srnamc:
//
//  infot   intEGER      Expected integer return code
//  nout    intEGER      Unit number for printing error messages
//  ok      LOGICAL      Set to .TRUE. if info = infot and
//                       srname = srnamt, otherwise set to .FALSE.
//  lerr    LOGICAL      Set to .TRUE., indicating that Xerbla was called
//  srnamt  CHARACTER*(*) Expected name of calling subroutine
// \endverbatim
//
//  =====================================================================
func Xerbla(srname *[]byte, info *int) {
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	infot := new(int)
	nout := new(int)
	common.infoc.lerr = new(bool)
	common.infoc.ok = new(bool)
	common.infoc.nout = new(int)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new([]byte)
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
	//     .. Scalars in common ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. common blocks ..
	infot = common.infoc.infot
	nout = common.infoc.nout
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Executable Statements ..
	//
	(*lerr) = true
	if (*(info)) != (*infot) {
		if (*infot) != 0 {
			WRITE((*nout), *func() *[]byte {
				y := []byte(" *** Xerbla was called from %s with info = %6d instead of %2d ***\n")
				return &y
			}(), (*srnamt)[1:lenTrim(srnamt)-(1)], (*(info)), (*infot))
		} else {
			WRITE((*nout), *func() *[]byte {
				y := []byte(" *** On entry to %s parameter number %6d had an illegal value ***\n")
				return &y
			}(), (*(srname))[1:lenTrim((srname))-(1)], (*(info)))
		}
		(*ok) = false
	}
	if bytes.Compare(*srname, *srnamt) != 0 {
		WRITE((*nout), *func() *[]byte {y := []byte(" *** Xerbla was called with srname = %s instead of %9s ***\n"); return &y }(), (*(srname))[1:lenTrim((srname))-(1)], (*srnamt)[1:lenTrim(srnamt)-(1)])
		(*ok) = false
	}
	return
	//
	//
	//     End of Xerbla
	//
}
