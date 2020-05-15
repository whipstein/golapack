package goblas

import 

// DerrorhrCol tests the error exits for DorhrCol that does
// Householder re_construction from the ouput of tall-skinny
// factorization DLAtsQR.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE DerrorhrCol( path, nunit)
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
// DerrorhrCol tests the error exits for DorhrCol that does
// Householder re_construction from the ouput of tall-skinny
// factorization DLAtsQR.
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
// \date November 2019
//
// \ingroup double_lin
//
//  =====================================================================
func DerrorhrCol(path *int, nunit *int) {
	_len := func() *[]byte {
		arr := make([]byte, -1)
		return &arr
	}()
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
	d := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	_len := func() *[]byte {
		arr := make([]byte, -1)
		return &arr
	}()
	infot := new(int)
	nout := new(int)
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.nout = new(float64)
	common.infoc.infot = new(int)
	srnamt := new(int)
	common.srnamc.srnamt = new(int)
	//
	//  -- lapACK test routine (version 3.9.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2019
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
			(*a)[(*i)-(1)][(*j)-(1)] = 1.e+0 / DBLE((*i)+(*j))
			(*t)[(*i)-(1)][(*j)-(1)] = 1.e+0 / DBLE((*i)+(*j))
		}
		(*d)[(*j)-(1)] = 0.e+0
	}
	(*ok) = true
	//
	//     Error exits for Householder re_construction
	//
	//     DorhrCol
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DorhrCol"); return &y }()
	//
	(*infot) = 1
	DorhrCol(-1, func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	(*infot) = 2
	DorhrCol(func() *int {y := 0; return &y }(), -1, func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	DorhrCol(func() *int {y := 1; return &y }(), func() *int {y := 2; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	(*infot) = 3
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), -1, a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	(*infot) = 5
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, -1, t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 0; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	DorhrCol(func() *int {y := 2; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	(*infot) = 7
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, -1, d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	DorhrCol(func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 1; return &y }(), a, func() *int {y := 1; return &y }(), t, func() *int {y := 0; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	DorhrCol(func() *int {y := 4; return &y }(), func() *int {y := 3; return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 4; return &y }(), t, func() *int {y := 1; return &y }(), d, info)
	Chkxer(func() *[]byte {y :=[]byte("DorhrCol"); return &y }(), infot, nout, lerr, ok)
	//
	//     Print a summary line.
	//
	Alaesm(path, ok, nout)
	//
	return
	//
	//     End of DerrorhrCol
	//
}
