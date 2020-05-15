package goblas

import 

// Derrql tests the error exits for the DOUBLE PRECISION routines
// that use the QL decomposition of a general matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrql( path, nunit)
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
// Derrql tests the error exits for the DOUBLE PRECISION routines
// that use the QL decomposition of a general matrix.
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
func Derrql(path *[]byte, nunit *int) {
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
	//     .. Scalars in Common ..
	//     ..
	//     .. Common blocks ..
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
	WRITE((*nout), *func() *[]byte {y :=[]byte(" %v\n"); return &y}())
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
	//     Error exits for QL factorization
	//
	//     Dgeqlf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqlf"); return &y}()
	(*infot) = 1
	Dgeqlf(-1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqlf"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqlf(func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), b, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqlf"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqlf(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqlf"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgeqlf(func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 1; return &y}(), b, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqlf"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dgeql2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeql2"); return &y}()
	(*infot) = 1
	Dgeql2(-1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeql2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeql2(func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeql2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeql2(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeql2"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dgeqls
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqls"); return &y}()
	(*infot) = 1
	Dgeqls(-1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqls(func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqls(func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgeqls(func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgeqls(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 2; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgeqls(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 10
	Dgeqls(func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 1; return &y}(), x, b, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqls"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dorgql
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgql"); return &y}()
	(*infot) = 1
	Dorgql(-1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgql(func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorgql(func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 2; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgql(func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorgql(func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorgql(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 8
	Dorgql(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), x, w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dorgql"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dorg2l
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorg2l"); return &y}()
	(*infot) = 1
	Dorg2l(-1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorg2l(func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorg2l(func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorg2l(func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorg2l(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorg2l(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorg2l"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dormql
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dormql"); return &y}()
	(*infot) = 1
	Dormql(func() *byte {y := byte('/'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dormql(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 2; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dormql(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 10
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormql(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 12
	Dormql(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 2; return &y}(), w, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dormql"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dorm2l
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorm2l"); return &y}()
	(*infot) = 1
	Dorm2l(func() *byte {y := byte('/'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dorm2l(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 2; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dorm2l(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 10
	Dorm2l(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), x, af, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dorm2l"); return &y}(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrql
	//
}
