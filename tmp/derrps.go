package goblas

import 

// Derrps tests the error exits for the DOUBLE PRECISION routines
// for DPSTRF.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Derrps( path, nunit)
//
//       .. Scalar Arguments ..
//       intEGER            nunit
//       CHARACTER*3        path
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Derrps tests the error exits for the DOUBLE PRECISION routines
// for DPSTRF.
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
func Derrps(path *[]byte, nunit *int) {
	nmax := new(int)
	i := new(int)
	info := new(int)
	j := new(int)
	rank := new(int)
	a := func() *[][]float64 {
		arr := make([][]float64, 4)
		for u := 0; u < 4; u++ {
			arr[u] = make([]float64, 4)
		}
		return &arr
	}()
	work := func() *[]float64 {
		arr := make([]float64, -1)
		return &arr
	}()
	piv := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	infot := new(int)
	nout := new(int)
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
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
	(*nmax) = 4
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
			//
			//Label100:
		}
		(*piv)[(*j)-(1)] = (*j)
		(*work)[(*j)-(1)] = 0.
		(*work)[(*nmax)+(*j)-(1)] = 0.
		//
		//Label110:
	}
	(*ok) = true
	//
	//
	//        Test error exits of the routines that use the Cholesky
	//        decomposition of a symmetric positive semidefinite matrix.
	//
	//        DPSTRF
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DPSTRF"); return &y }()
	(*infot) = 1
	Dpstrf(func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("DPSTRF"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dpstrf(func() *byte {y := byte('U'); return &y }(), -1, a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("DPSTRF"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dpstrf(func() *byte {y := byte('U'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("DPSTRF"); return &y }(), infot, nout, lerr, ok)
	//
	//        Dpstf2
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dpstf2"); return &y }()
	(*infot) = 1
	Dpstf2(func() *byte {y := byte('/'); return &y }(), func() *int {y := 0; return &y }(), a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("Dpstf2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 2
	Dpstf2(func() *byte {y := byte('U'); return &y }(), -1, a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("Dpstf2"); return &y }(), infot, nout, lerr, ok)
	(*infot) = 4
	Dpstf2(func() *byte {y := byte('U'); return &y }(), func() *int {y := 2; return &y }(), a, func() *int {y := 1; return &y }(), piv, rank, -1., work, info)
	Chkxer(func() *[]byte {y :=[]byte("Dpstf2"); return &y }(), infot, nout, lerr, ok)
	//
	//
	//     Print a summary line.
	//
	Alaesm((path), ok, nout)
	//
	return
	//
	//     End of Derrps
	//
}
