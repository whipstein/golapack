package goblas

import 

// Derrls tests the error exits for the DOUBLE PRECISION least squares
// driver routines (Dgels, SGELSS, SGELSY, SGELSD).
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrls( path, nunit)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        path
//       inTEGER            nunit
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Derrls tests the error exits for the DOUBLE PRECISION least squares
// driver routines (Dgels, SGELSS, SGELSY, SGELSD).
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
//          nunit is inTEGER
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
func Derrls(path *[]byte, nunit *int) {
	nmax := new(int)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	info := new(int)
	irnk := new(int)
	rcond := new(float64)
	ip := func() *[]int {
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
	b := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	s := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	w := func() *[]float64 {
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
	//     .. Scalars in Common ..
	//     ..
	//     .. Common blocks ..
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
	(*a)[0][0] = 1.0e+0
	(*a)[0][1] = 2.0e+0
	(*a)[1][1] = 3.0e+0
	(*a)[1][0] = 4.0e+0
	(*ok) = true
	//
	if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("LS"); return &y }()) {
		//
		//        Test error exits for the least squares driver routines.
		//
		//        Dgels
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgels "); return &y }()
		(*infot) = 1
		Dgels(func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgels(func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgels(func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dgels(func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 6
		Dgels(func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), w, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 8
		Dgels(func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		(*infot) = 10
		Dgels(func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgels "); return &y }(), infot, nout, lerr, ok)
		//
		//        Dgelss
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgelss"); return &y }()
		(*infot) = 1
		Dgelss(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelss"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgelss(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelss"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgelss(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelss"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dgelss(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), S, rcond, IRNK, w, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelss"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dgelss(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 2; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelss"); return &y }(), infot, nout, lerr, ok)
		//
		//        Dgelsy
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgelsy"); return &y }()
		(*infot) = 1
		Dgelsy(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), ip, rcond, IRNK, w, func() *int {y := 10; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgelsy(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), ip, rcond, IRNK, w, func() *int {y := 10; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgelsy(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), ip, rcond, IRNK, w, func() *int {y := 10; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dgelsy(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), ip, rcond, IRNK, w, func() *int {y := 10; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dgelsy(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), ip, rcond, IRNK, w, func() *int {y := 10; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 12
		Dgelsy(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), ip, rcond, IRNK, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), infot, nout, lerr, ok)
		//
		//        DgelsD
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("DgelsD"); return &y }()
		(*infot) = 1
		Dgelsd(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 10; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgelsd(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 10; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgelsd(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 10; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 5
		Dgelsd(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), S, rcond, IRNK, w, func() *int {y := 10; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dgelsd(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), S, rcond, IRNK, w, func() *int {y := 10; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 12
		Dgelsd(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), S, rcond, IRNK, w, func() *int {y := 1; return &y }(), ip, info)
		Chkxer(func() *[]byte {y :=[]byte("DgelsD"); return &y }(), infot, nout, lerr, ok)
	}
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrls
	//
}
