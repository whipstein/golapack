package goblas

import 

// Derrrq tests the error exits for the DOUBLE PRECISION routines
// that use the RQ decomposition of a general matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrrq( path, nunit)
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
// Derrrq tests the error exits for the DOUBLE PRECISION routines
// that use the RQ decomposition of a general matrix.
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
func Derrrq(path *[]byte, nunit *int) {
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
			(*a)[(*i)-1][(*j)-1] = 1. / DBLE((*i)+(*j))
			(*af)[(*i)-1][(*j)-1] = 1. / DBLE((*i)+(*j))
			//Label10:
		}
		(*b)[(*j)-1] = 0.
		(*w)[(*j)-1] = 0.
		(*x)[(*j)-1] = 0.
		//Label20:
	}
	(*ok) = true
	//
	//     Error exits for RQ factorization
	//
	//     Dgerqf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgerqf"); return &y }()
	(*infot) = 1
	Dgerqf(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgerqf(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgerqf(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgerqf(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 2; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqf"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgerq2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgerq2"); return &y }()
	(*infot) = 1
	Dgerq2(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerq2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgerq2(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerq2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgerq2(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerq2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgerqs
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgerqs"); return &y }()
	(*infot) = 1
	Dgerqs(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgerqs(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgerqs(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgerqs(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgerqs(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgerqs(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dgerqs(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorgrq
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgrq"); return &y }()
	(*infot) = 1
	Dorgrq(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgrq(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgrq(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgrq(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgrq(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorgrq(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dorgrq(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgrq"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorgr2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgr2"); return &y }()
	(*infot) = 1
	Dorgr2(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgr2(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgr2(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgr2(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgr2(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 2; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorgr2(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgr2"); return &y }(), infot, nout, lerr, ok)
	//
	//     DORMRQ
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORMRQ"); return &y }()
	(*infot) = 1
	Dormrq(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormrq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormrq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormrq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormrq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("DORMRQ"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dormr2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dormr2"); return &y }()
	(*infot) = 1
	Dormr2(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormr2(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormr2(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dormr2(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dormr2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrrq
	//
}
