package goblas

import 

// Derrgt tests the error exits for the DOUBLE PRECISION tridiagonal
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
//       SUBROUTinE DerrGt( path, nunit)
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
// Derrgt tests the error exits for the DOUBLE PRECISION tridiagonal
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
func DerrGt(path *[]byte, nunit *int) {
	nmax := new(int)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	info := new(int)
	anorm := new(float64)
	rcond := new(float64)
	ip := func() *[]int {
		arr := make([]int, 2)
		return &arr
	}()
	iw := func() *[]int {
		arr := make([]int, 2)
		return &arr
	}()
	b := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	c := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	cf := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	d := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	df := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	e := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	ef := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	f := func() *[]float64 {
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
	WRITE((*nout), *func() *[]byte {y :=[]byte(" %v\n"); return &y}())
	(*c2) = (*(path))[1]
	(*d)[0] = 1.
	(*d)[1] = 2.
	(*df)[0] = 1.
	(*df)[1] = 2.
	(*e)[0] = 3.
	(*e)[1] = 4.
	(*ef)[0] = 3.
	(*ef)[1] = 4.
	(*anorm) = 1.0
	(*ok) = true
	//
	if Lsamen(func() *int {y := 2; return &y}(), c2, func() *[]byte {y :=[]byte("GT"); return &y}()) {
		//
		//        Test error exits for the general tridiagonal routines.
		//
		//        Dgttrf
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgttrf"); return &y}()
		(*infot) = 1
		Dgttrf(-1, c, d, e, f, ip, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgttrf"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dgttrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgttrs"); return &y}()
		(*infot) = 1
		Dgttrs(func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), c, d, e, f, ip, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgttrs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgttrs(func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), c, d, e, f, ip, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgttrs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgttrs(func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, c, d, e, f, ip, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgttrs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 10
		Dgttrs(func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), c, d, e, f, ip, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dgttrs"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dgtrfs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgtrfs"); return &y}()
		(*infot) = 1
		Dgtrfs(func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), c, d, e, cf, df, ef, f, ip, b, func() *int {y := 1; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgtrfs(func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), c, d, e, cf, df, ef, f, ip, b, func() *int {y := 1; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 3
		Dgtrfs(func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, c, d, e, cf, df, ef, f, ip, b, func() *int {y := 1; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 13
		Dgtrfs(func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), c, d, e, cf, df, ef, f, ip, b, func() *int {y := 1; return &y}(), x, func() *int {y := 2; return &y}(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 15
		Dgtrfs(func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), c, d, e, cf, df, ef, f, ip, b, func() *int {y := 2; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtrfs"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dgtcon
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgtcon"); return &y}()
		(*infot) = 1
		Dgtcon(func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), c, d, e, f, ip, anorm, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtcon"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgtcon(func() *byte {y := byte('I'); return &y}(), -1, c, d, e, f, ip, anorm, rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtcon"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 8
		Dgtcon(func() *byte {y := byte('I'); return &y}(), func() *int {y := 0; return &y}(), c, d, e, f, ip, -(*anorm), rcond, w, iw, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgtcon"); return &y}(), infot, nout, lerr, ok)
		//
	} else if Lsamen(func() *int {y := 2; return &y}(), c2, func() *[]byte {y :=[]byte("PT"); return &y}()) {
		//
		//        Test error exits for the positive definite tridiagonal
		//        routines.
		//
		//        Dpttrf
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dpttrf"); return &y}()
		(*infot) = 1
		Dpttrf(-1, d, e, info)
		Chkxer(func() *[]byte {y :=[]byte("Dpttrf"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dpttrs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dpttrs"); return &y}()
		(*infot) = 1
		Dpttrs(-1, func() *int {y := 0; return &y}(), d, e, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dpttrs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dpttrs(func() *int {y := 0; return &y}(), -1, d, e, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dpttrs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 6
		Dpttrs(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), d, e, x, func() *int {y := 1; return &y}(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dpttrs"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dptrfs
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dptrfs"); return &y}()
		(*infot) = 1
		Dptrfs(-1, func() *int {y := 0; return &y}(), d, e, df, ef, b, func() *int {y := 1; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dptrfs(func() *int {y := 0; return &y}(), -1, d, e, df, ef, b, func() *int {y := 1; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 8
		Dptrfs(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), d, e, df, ef, b, func() *int {y := 1; return &y}(), x, func() *int {y := 2; return &y}(), r1, r2, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptrfs"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 10
		Dptrfs(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), d, e, df, ef, b, func() *int {y := 2; return &y}(), x, func() *int {y := 1; return &y}(), r1, r2, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptrfs"); return &y}(), infot, nout, lerr, ok)
		//
		//        Dptcon
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dptcon"); return &y}()
		(*infot) = 1
		Dptcon(-1, d, e, anorm, rcond, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptcon"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 4
		Dptcon(func() *int {y := 0; return &y}(), d, e, -(*anorm), rcond, w, info)
		Chkxer(func() *[]byte {y :=[]byte("Dptcon"); return &y}(), infot, nout, lerr, ok)
	}
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrgt
	//
}
