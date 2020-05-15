package goblas

import 

// Derrqrtp tests the error exits for the REAL routines
// that use the QRT decomposition of a triangular-pentagonal matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrqrtp( path, nunit)
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
// Derrqrtp tests the error exits for the REAL routines
// that use the QRT decomposition of a triangular-pentagonal matrix.
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
func Derrqrtp(path *[]byte, nunit *int) {
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
			(*c)[(*i)-(1)][(*j)-(1)] = 1. / DBLE((*i)+(*j))
			(*t)[(*i)-(1)][(*j)-(1)] = 1. / DBLE((*i)+(*j))
		}
		(*w)[(*j)-(1)] = 0.0
	}
	(*ok) = true
	//
	//     Error exits for TPQRT factorization
	//
	//     DTPQRT
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DTPQRT"); return &y }()
	(*infot) = 1
	Dtpqrt(-1, func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtpqrt(func() *int {y := 1; return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtpqrt(func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), -1, func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtpqrt(func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtpqrt(func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtpqrt(func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 6
	Dtpqrt(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 8
	Dtpqrt(func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 10
	Dtpqrt(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), t, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("DTPQRT"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dtpqrt2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }()
	(*infot) = 1
	Dtpqrt2(-1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtpqrt2(func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtpqrt2(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dtpqrt2(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), b, func() *int {y := 2; return &y }(), t, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dtpqrt2(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 1; return &y }(), t, func() *int {y := 2; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 9
	Dtpqrt2(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 2; return &y }(), b, func() *int {y := 2; return &y }(), t, func() *int {y := 1; return &y }(), info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpqrt2"); return &y }(), infot, nout, lerr, ok)
	//
	//     Dtpmqrt
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }()
	(*infot) = 1
	Dtpmqrt(func() *byte {y := byte('/'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 3
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 5
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	(*infot) = 6
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 7
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 9
	Dtpmqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 9
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 11
	Dtpmqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 0; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 13
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 0; return &y }(), c, func() *int {y := 1; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 15
	Dtpmqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), b, func() *int {y := 1; return &y }(), c, func() *int {y := 0; return &y }(), w, info)
	Chkxer(func() *[]byte {y :=[]byte("Dtpmqrt"); return &y }(), infot, nout, lerr, ok)
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
