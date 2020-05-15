package goblas

import 

// Derrqp tests the error exits for Dgeqp3.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrqp( path, nunit)
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
// Derrqp tests the error exits for Dgeqp3.
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
func Derrqp(path *[]byte, nunit *int) {
	nmax := new(int)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	info := new(int)
	lw := new(int)
	ip := func() *[]int {
		arr := make([]int, 3)
		return &arr
	}()
	a := func() *[][]float64 {
		arr := make([][]float64, 3)
		for u := 0; u < 3; u++ {
			arr[u] = make([]float64, 3)
		}
		return &arr
	}()
	tau := func() *[]float64 {
		arr := make([]float64, 3)
		return &arr
	}()
	w := func() *[]float64 {
		arr := make([]float64, -1)
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
	(*nmax) = 3
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
	(*LW) = 3*(*nmax) + 1
	(*a)[0][0] = 1.0e+0
	(*a)[0][1] = 2.0e+0
	(*a)[1][1] = 3.0e+0
	(*a)[1][0] = 4.0e+0
	(*ok) = true
	//
	if Lsamen(func() *int {y := 2; return &y}(), c2, func() *[]byte {y :=[]byte("QP"); return &y}()) {
		//
		//        Test error exits for QR factorization with pivoting
		//
		//        Dgeqp3
		//
		(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqp3"); return &y}()
		(*infot) = 1
		Dgeqp3(-1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), ip, tau, w, LW, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgeqp3"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 2
		Dgeqp3(func() *int {y := 1; return &y}(), -1, a, func() *int {y := 1; return &y}(), ip, tau, w, LW, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgeqp3"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 4
		Dgeqp3(func() *int {y := 2; return &y}(), func() *int {y := 3; return &y}(), a, func() *int {y := 1; return &y}(), ip, tau, w, LW, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgeqp3"); return &y}(), infot, nout, lerr, ok)
		(*infot) = 8
		Dgeqp3(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), ip, tau, w, (*LW)-10, info)
		Chkxer(func() *[]byte {y :=[]byte("Dgeqp3"); return &y}(), infot, nout, lerr, ok)
	}
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrqp
	//
}
