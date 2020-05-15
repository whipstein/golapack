package goblas

import 

// Derrqrt tests the error exits for the DOUBLE PRECISION routines
// that use the QRT decomposition of a general matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrqrt( path, nunit)
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
// Derrqrt tests the error exits for the DOUBLE PRECISION routines
// that use the QRT decomposition of a general matrix.
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
func Derrqrt(path *[]byte, nunit *int) {
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
	t := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	w := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	c := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
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
			(*c)[(*i)-1][(*j)-1] = 1. / DBLE((*i)+(*j))
			(*t)[(*i)-1][(*j)-1] = 1. / DBLE((*i)+(*j))
		}
		(*w)[(*j)-1] = 0.
	}
	(*ok) = true
	//
	//     Error exits for QRT factorization
	//
	//     Dgeqrt
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrt"); return &y}()
	(*infot) = 1
	Dgeqrt(-1, func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrt(func() *int {y := 0; return &y}(), -1, func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgeqrt(func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgeqrt(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dgeqrt(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dgeqrt2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrt2"); return &y}()
	(*infot) = 1
	Dgeqrt2(-1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrt2(func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqrt2(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 6
	Dgeqrt2(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt2"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dgeqrt3
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqrt3"); return &y}()
	(*infot) = 1
	Dgeqrt3(-1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt3"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgeqrt3(func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt3"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgeqrt3(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt3"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 6
	Dgeqrt3(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dgeqrt3"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dgemqrt
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y}()
	(*infot) = 1
	Dgemqrt(func() *byte {y := byte('/'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dgemqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 6
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgemqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 8
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 10
	Dgemqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 0; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 12
	Dgemqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), c, func() *int {y := 0; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dgemqrt"); return &y}(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrqrt
	//
}
