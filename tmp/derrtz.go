package goblas

import 

// Derrtz tests the error exits for STZRZF.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrtz( path, nunit)
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
// Derrtz tests the error exits for STZRZF.
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
func Derrtz(path *[]byte, nunit *int) {
	nmax := new(int)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	info := new(int)
	a := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	tau := func() *[]float64 {
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
	(*a)[0][0] = 1.e+0
	(*a)[0][1] = 2.e+0
	(*a)[1][1] = 3.e+0
	(*a)[1][0] = 4.e+0
	(*w)[0] = 0.0e+0
	(*w)[1] = 0.0e+0
	(*ok) = true
	//
	if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y :=[]byte("TZ"); return &y }()) {
		//
		//        Test error exits for the trapezoidal routines.
		//
		//        Dtzrzf
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dtzrzf"); return &y }()
		(*infot) = 1
		Dtzrzf(-1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), tau, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtzrzf"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 2
		Dtzrzf(func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), tau, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtzrzf"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 4
		Dtzrzf(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), tau, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtzrzf"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dtzrzf(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 2; return &y }(), tau, w, func() *int {y := 0; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtzrzf"); return &y }(), infot, nout, lerr, ok)
		(*infot) = 7
		Dtzrzf(func() *int {y := 2; return &y }(), func() *int {y := 3; return &y }(), a, func() *int {y := 2; return &y }(), tau, w, func() *int {y := 1; return &y }(), info)
		Chkxer(func() *[]byte {y :=[]byte("Dtzrzf"); return &y }(), infot, nout, lerr, ok)
	}
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrtz
	//
}
