package goblas

import 

// Derrlq tests the error exits for the DOUBLE PRECISION routines
// that use the LQ decomposition of a general matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrlq( path, nunit)
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
// Derrlq tests the error exits for the DOUBLE PRECISION routines
// that use the LQ decomposition of a general matrix.
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
func Derrlq(path *[]byte, nunit *int) {
	nmax := new(int)
	i := new(int)
	info := new(int)
	j := new(int)
	a := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	af := func() *[][]float64 {
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
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	(*nout) = (*(nunit))
	WRITE((*nout), *func() *[]byte {y :=[]byte(" %v\n"); return &y }())
	//
	//     Set the variables to innocuous values.
	//
	for (*j) = 1; (*j) <= (*nmax); (*j)++ {
		for (*i) = 1; (*i) <= (*nmax); (*i)++ {
			(*a)[(*i)-(1)][(*j)-(1)] = 1. / DBLE((*i)+(*j))
			(*af)[(*i)-(1)][(*j)-(1)] = 1. / DBLE((*i)+(*j))
			//Label10:
		}
		(*b)[(*j)-(1)] = 0.
		(*w)[(*j)-(1)] = 0.
		(*x)[(*j)-(1)] = 0.
		//Label20:
	}
	(*ok) = true
	//
	//     Error exits for LQ factorization
	//
	//     Dgelqf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgelqf"); return &y }()
	(*infot) = 1
	Dgelqf(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgelqf(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgelqf(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgelqf(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqf"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgelq2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgelq2"); return &y }()
	(*infot) = 1
	Dgelq2(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelq2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgelq2(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelq2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgelq2(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelq2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgelqs
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgelqs"); return &y }()
	(*infot) = 1
	Dgelqs(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgelqs(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgelqs(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgelqs(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgelqs(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgelqs(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dgelqs(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), infot, nout, lerr, ok)
	//
	//     DORGLQ
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORGLQ"); return &y }()
	(*infot) = 1
	Dorglq(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorglq(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorglq(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorglq(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorglq(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorglq(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dorglq(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORGLQ"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorgl2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgl2"); return &y }()
	(*infot) = 1
	Dorgl2(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgl2(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgl2(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgl2(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgl2(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorgl2(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgl2"); return &y }(), infot, nout, lerr, ok)
	//
	//     DORMLQ
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORMLQ"); return &y }()
	(*infot) = 1
	Dormlq(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMLQ"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorml2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorml2"); return &y }()
	(*infot) = 1
	Dorml2(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorml2(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorml2(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dorml2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorml2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrlq
	//
}
