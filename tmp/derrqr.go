package goblas

import 

// Derrqr tests the error exits for the DOUBLE PRECISION routines
// that use the QR decomposition of a general matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrqr( path, nunit)
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
// Derrqr tests the error exits for the DOUBLE PRECISION routines
// that use the QR decomposition of a general matrix.
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
func Derrqr(path *[]byte, nunit *int) {
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
	//     Error exits for QR factorization
	//
	//     Dgeqrf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrf"); return &y }()
	(*infot) = 1
	Dgeqrf(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrf(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqrf(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrf"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgeqrf(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrf"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgeqrfp
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrfp"); return &y }()
	(*infot) = 1
	Dgeqrfp(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrfp"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrfp(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrfp"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqrfp(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrfp"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgeqrfp(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), b, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrfp"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgeqr2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqr2"); return &y }()
	(*infot) = 1
	Dgeqr2(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqr2(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqr2(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgeqr2p
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqr2p"); return &y }()
	(*infot) = 1
	Dgeqr2p(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2p"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqr2p(func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2p"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqr2p(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqr2p"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dgeqrs
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrs"); return &y }()
	(*infot) = 1
	Dgeqrs(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrs(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrs(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, b, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgeqrs(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgeqrs(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgeqrs(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dgeqrs(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, b, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorgqr
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgqr"); return &y }()
	(*infot) = 1
	Dorgqr(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgqr(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgqr(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgqr(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgqr(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorgqr(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dorgqr(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgqr"); return &y }(), infot, nout, lerr, ok)
	//
	//     DORG2R
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORG2R"); return &y }()
	(*infot) = 1
	Dorg2r(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorg2r(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorg2r(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorg2r(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorg2r(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 2; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorg2r(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("DORG2R"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dormqr
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dormqr"); return &y }()
	(*infot) = 1
	Dormqr(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormqr"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dorm2r
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorm2r"); return &y }()
	(*infot) = 1
	Dorm2r(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2r(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 2; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorm2r(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dorm2r(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), x, af, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2r"); return &y }(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrqr
	//
}
