package goblas

import 

// Derrtr tests the error exits for the DOUBLE PRECISION triangular
// routines.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrtr( path, nunit)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        path
//       intEGER            nunit
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Derrtr tests the error exits for the DOUBLE PRECISION triangular
// routines.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The lapACK path name for the routines to be tested.
// \endverbatim
//
// \param[in] nunit
// \verbatim
//          nunit is intEGER
//          The unit number for output.
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
func Derrtr(path *[]byte, nunit *int) {
	nmax := new(int)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	info := new(int)
	rcond := new(float64)
	scale := new(float64)
	iw := func() *[]int {
		arr := make([]int, 2)
		return &arr
	}()
	a := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	b := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	r1 := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	r2 := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	w := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	x := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	infot := new(int)
	nout := new(int)
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.nout = new(float64)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new(int)
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
	(*nmax) = 2
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Scalars in common ..
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
	(*nout) = (*(nunit))
	WRITE((*nout), *func() *[]byte {y :=[]byte(" %v\n"); return &y }())
	(*c2) = (*(path))[1]
	(*a)[0][0] = 1.
	(*a)[0][1] = 2.
	(*a)[1][1] = 3.
	(*a)[1][0] = 4.
	(*ok) = true
	//
	if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TR"); return &y }()) {
		//
		//        Test error exits for the general triangular routines.
		//
		//        Dtrtri
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtrtri"); return &y }()
		(*infot) = 1
		Dtrtri(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtri"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtrtri(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtri"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtrtri(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtri"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtrtri(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtri"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtrti2
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtrti2"); return &y }()
		(*infot) = 1
		Dtrti2(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrti2"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtrti2(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrti2"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtrti2(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrti2"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtrti2(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrti2"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtrtrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtrtrs"); return &y }()
		(*infot) = 1
		Dtrtrs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 9
		Dtrtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrtrs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtrrfs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtrrfs"); return &y }()
		(*infot) = 1
		Dtrrfs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), x, func() *int {y := 2; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 9
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 2; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 11
		Dtrrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrrfs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtrcon
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtrcon"); return &y }()
		(*infot) = 1
		Dtrcon(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtrcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtrcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtrcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 6
		Dtrcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtrcon"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dlatrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dlatrs"); return &y }()
		(*infot) = 1
		Dlatrs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dlatrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dlatrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dlatrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dlatrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dlatrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatrs"); return &y }(), infot, nout, lerr, ok)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TP"); return &y }()) {
		//
		//        Test error exits for the packed triangular routines.
		//
		//        Dtptri
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtptri"); return &y }()
		(*infot) = 1
		Dtptri(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptri"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtptri(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptri"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtptri(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptri"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtptrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtptrs"); return &y }()
		(*infot) = 1
		Dtptrs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtptrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtptrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtptrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtptrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dtptrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtptrs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtprfs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtprfs"); return &y }()
		(*infot) = 1
		Dtprfs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, b, func() *int {y := 1; return &y }(), x, func() *int {y := 2; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 10
		Dtprfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, b, func() *int {y := 2; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtprfs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtpcon
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtpcon"); return &y }()
		(*infot) = 1
		Dtpcon(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtpcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtpcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtpcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtpcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtpcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtpcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtpcon"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dlatps
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dlatps"); return &y }()
		(*infot) = 1
		Dlatps(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatps"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dlatps(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatps"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dlatps(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), a, x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatps"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dlatps(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatps"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dlatps(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, a, x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatps"); return &y }(), infot, nout, lerr, ok)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TB"); return &y }()) {
		//
		//        Test error exits for the banded triangular routines.
		//
		//        Dtbtrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtbtrs"); return &y }()
		(*infot) = 1
		Dtbtrs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 6
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 10
		Dtbtrs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbtrs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtbrfs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtbrfs"); return &y }()
		(*infot) = 1
		Dtbrfs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 6
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), x, func() *int {y := 2; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 10
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), x, func() *int {y := 2; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 12
		Dtbrfs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), x, func() *int {y := 1; return &y }(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbrfs"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dtbcon
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtbcon"); return &y }()
		(*infot) = 1
		Dtbcon(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtbcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dtbcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtbcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dtbcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dtbcon(func() *byte {y := byte('1'); return &y }(), func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dtbcon"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dlatbs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dlatbs"); return &y }()
		(*infot) = 1
		Dlatbs(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 6
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dlatbs(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, scale, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dlatbs"); return &y }(), infot, nout, lerr, ok)
	}
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrtr
	//
}
