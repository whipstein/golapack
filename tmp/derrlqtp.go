package goblas

import 

// Derrlqtp tests the error exits for the REAL routines
// that use the LQT decomposition of a triangular-pentagonal matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrlqtp( path, nunit)
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
// Derrlqtp tests the error exits for the REAL routines
// that use the LQT decomposition of a triangular-pentagonal matrix.
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
func Derrlqtp(path *[]byte, nunit *int) {
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
	b := func() *[][]float64 {
		arr := make([][]float64, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]float64, 2)
		}
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
		(*w)[(*j)-1] = 0.0
	}
	(*ok) = true
	//
	//     Error exits for TPLQT factorization
	//
	//     DTPLQT
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DTPLQT"); return &y}()
	(*infot) = 1
	Dtplqt(-1, func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtplqt(func() *int {y := 1; return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtplqt(func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), -1, func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtplqt(func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtplqt(func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtplqt(func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 6
	Dtplqt(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 8
	Dtplqt(func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 2; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 10
	Dtplqt(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 2; return &y}(), a, func() *int {y := 2; return &y}(), b, func() *int {y := 2; return &y}(), t, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPLQT"); return &y}(), infot, nout, lerr, ok)
	//
	//     Dtplqt2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dtplqt2"); return &y}()
	(*infot) = 1
	Dtplqt2(-1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtplqt2(func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtplqt2(func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, a, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dtplqt2(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), b, func() *int {y := 2; return &y}(), t, func() *int {y := 2; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dtplqt2(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), b, func() *int {y := 1; return &y}(), t, func() *int {y := 2; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 9
	Dtplqt2(func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 2; return &y}(), b, func() *int {y := 2; return &y}(), t, func() *int {y := 1; return &y}(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtplqt2"); return &y}(), infot, nout, lerr, ok)
	//
	//     DTPMLQT
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DTPMLQT"); return &y}()
	(*infot) = 1
	Dtpmlqt(func() *byte {y := byte('/'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('/'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 5
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 0; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	(*infot) = 6
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), -1, func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 7
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 9
	Dtpmlqt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 2; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 11
	Dtpmlqt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 0; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 13
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 0; return &y}(), c, func() *int {y := 1; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	(*infot) = 15
	Dtpmlqt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), a, func() *int {y := 1; return &y}(), t, func() *int {y := 1; return &y}(), b, func() *int {y := 1; return &y}(), c, func() *int {y := 0; return &y}(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPMLQT"); return &y}(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrlqt
	//
}
