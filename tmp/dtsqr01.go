package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtsqr01 tests DGEQR, DGELq, DGEMLQ and DGEMQR.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dtsqr01(tssw, m,N, mb, nb, result)
//
//       .. Scalar Arguments ..
//       intEGER m, n, mb
//       .. Return values ..
//       DOUBLE PRECISION result(6)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtsqr01 tests DGEQR, DGELq, DGEMLQ and DGEMQR.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] tssw
// \verbatim
//          tssw is CHARACTER
//          'ts' for testing tall skinny QR
//               and anything else for testing short wide LQ
// \endverbatim
// \param[in] M
// \verbatim
//          M is intEGER
//          Number of rows in test matrix.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          Number of columns in test matrix.
// \endverbatim
// \param[in] mb
// \verbatim
//          mb is intEGER
//          Number of row in row block in test matrix.
// \endverbatim
//
// \param[in] nb
// \verbatim
//          nb is intEGER
//          Number of columns in column block test matrix.
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (6)
//          Results of each of the six tests below.
//
//          result1 = | A - Q R | or | A - L Q |
//          result(2) = | I - Q^H Q | or | I - Q Q^H |
//          result(3) = | Q C - Q C |
//          result(4) = | Q^H C - Q^H C |
//          result(5) = | C Q - C Q |
//          result(6) = | C Q^H - C Q^H |
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
// \date April 2012
//
// \ingroup double_lin
//
//  =====================================================================
func Dtsqr01(tssw *byte, m *int, n *int, mb *int, nb *int, result *[]float64) {
	allocatable := new(float64)
	q := new(float64)
	r := new(float64)
	rwork := new(float64)
	work := new(float64)
	t := new(float64)
	cf := new(float64)
	df := new(float64)
	a := new(float64)
	c := new(float64)
	d := new(float64)
	lq := new(float64)
	one := new(float64)
	zero := new(float64)
	testzeros := new(bool)
	ts := new(bool)
	info := new(int)
	j := new(int)
	k := new(int)
	l := new(int)
	lwork := new(int)
	tsize := new(int)
	mnb := new(int)
	anorm := new(float64)
	eps := new(float64)
	resid := new(float64)
	cnorm := new(float64)
	dnorm := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	tquery := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	workquery := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	common.srnamc.srnamt = new(int)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     April 2012
	//
	//     .. Scalar Arguments ..
	//     .. Return values ..
	//
	//  =====================================================================
	//
	//     ..
	//     .. Local allocatable arrays
	//
	//     .. Parameters ..
	(*zero) = 0.0
	(*one) = 1.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Data statements ..
	(*iseed)[0], (*iseed)[1], (*iseed)[2], (*iseed)[3] = 1988, 1989, 1990, 1991
	//
	//     TEST TALL SKInnY OR SHORT WIDE
	//
	(*ts) = (*blas.Lsame((tssw), func() *[]byte {y := []byte("ts"); return &y }()))
	//
	//     TEST MAtriCES WITH half OF MAtrix BEinG zeroS
	//
	(*testzeros) = false
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*k) = (Min((*(m)), (*(n))))
	(*l) = (*MAX((m), (n), func() *int {y := 1; return &y }()))
	(*mnb) = (MAX((*(mb)), (*(nb))))
	(*lwork) = MAX(3, (*l)) * (*mnb)
	//
	//     Dynamically allocate local arrays
	//
	ALLOCATE(a((m), (n)), af((m), (n)), q(L, L), r((m), L), rwork(l), c((m), (n)), cf((m), (n)), d((n), (m)), df((n), (m)), Lq(L, (n)))
	//
	//     Put random numbers into A and copy to af
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), a(func() *int {y := 1; return &y }(), j))
	}
	if *testzeros {
		if (*(m)) >= 4 {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), iseed, (*(m))/2, a((*(m))/4, j))
			}
		}
	}
	Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), a, (m), af, (m))
	//
	if *ts {
		//
		//     Factor the matrix A in the array af.
		//
		Dgeqr((m), (n), af, (m), tquery, -1, workquery, -1, info)
		(*tsize) = (*int(&((*tquery)[0])))
		(*lwork) = (*int(&((*workquery)[0])))
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (m), k, af, (m), tquery, tsize, cf, (m), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, af, (m), tquery, tsize, cf, (m), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, af, (m), tquery, tsize, cf, (m), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, af, (m), tquery, tsize, df, (n), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, af, (m), tquery, tsize, df, (n), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		ALLOCATE(t(tsize))
		ALLOCATE(work(lwork))
		(*srnamt) = *func() *[]byte {y := []byte("DGEQR"); return &y }()
		Dgeqr((m), (n), af, (m), t, tsize, work, lwork, info)
		//
		//     Generate the m-by-m matrix Q
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (m), (m), zero, one, q, (m))
		(*srnamt) = *func() *[]byte {y := []byte("DGEMQR"); return &y }()
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (m), k, af, (m), t, tsize, q, (m), work, lwork, info)
		//
		//     Copy R
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), zero, zero, r, (m))
		Dlacpy(func() *[]byte {y := []byte("Upper"); return &y }(), (m), (n), af, (m), r, (m))
		//
		//     Compute |R - Q'*a| / |A| and store in result1
		//
		Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), -(*one), q, (m), a, (m), one, r, (m))
		(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), a, (m), rwork))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), r, (m), rwork))
		if (*anorm) > (*zero) {
			(*(result))[0] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*anorm))
		} else {
			(*(result))[0] = (*zero)
		}
		//
		//     Compute |I - Q'*Q| and store in result(2)
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (m), (m), zero, one, r, (m))
		Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('C'); return &y }(), (m), (m), -(*one), q, (m), one, r, (m))
		(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y := []byte("Upper"); return &y }(), (m), r, (m), rwork))
		(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*(m))))
		//
		//     Generate random m-by-n matrix C and a copy cf
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), c(func() *int {y := 1; return &y }(), j))
		}
		(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), c, (m), rwork))
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
		//
		//     Apply Q to C as Q*C
		//
		(*srnamt) = *func() *[]byte {y := []byte("DGEMQR"); return &y }()
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, af, (m), t, tsize, cf, (m), work, lwork, info)
		//
		//     Compute |Q*C - Q*C| / |C|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), -(*one), q, (m), c, (m), one, cf, (m))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), cf, (m), rwork))
		if (*cnorm) > (*zero) {
			(*(result))[2] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*cnorm))
		} else {
			(*(result))[2] = (*zero)
		}
		//
		//     Copy C into cf again
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
		//
		//     Apply Q to C as QT*C
		//
		(*srnamt) = *func() *[]byte {y := []byte("DGEMQR"); return &y }()
		Dgemqr(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, af, (m), t, tsize, cf, (m), work, lwork, info)
		//
		//     Compute |QT*C - QT*C| / |C|
		//
		Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), -(*one), q, (m), c, (m), one, cf, (m))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), cf, (m), rwork))
		if (*cnorm) > (*zero) {
			(*(result))[3] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*cnorm))
		} else {
			(*(result))[3] = (*zero)
		}
		//
		//     Generate random n-by-m matrix D and a copy df
		//
		for (*j) = 1; (*j) <= (*(m)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (n), d(func() *int {y := 1; return &y }(), j))
		}
		(*dnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), d, (n), rwork))
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
		//
		//     Apply Q to D as D*Q
		//
		(*srnamt) = *func() *[]byte {y := []byte("DGEMQR"); return &y }()
		Dgemqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, af, (m), t, tsize, df, (n), work, lwork, info)
		//
		//     Compute |D*Q - D*Q| / |D|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (m), -(*one), d, (n), q, (m), one, df, (n))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
		if (*dnorm) > (*zero) {
			(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
		} else {
			(*(result))[4] = (*zero)
		}
		//
		//     Copy D into df again
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
		//
		//     Apply Q to D as D*QT
		//
		Dgemqr(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, af, (m), t, tsize, df, (n), work, lwork, info)
		//
		//     Compute |D*QT - D*QT| / |D|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), (m), -(*one), d, (n), q, (m), one, df, (n))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
		if (*cnorm) > (*zero) {
			(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
		} else {
			(*(result))[5] = (*zero)
		}
		//
		//     Short and wide
		//
	} else {
		Dgelq((m), (n), af, (m), tquery, -1, workquery, -1, info)
		(*tsize) = (*int(&((*tquery)[0])))
		(*lwork) = (*int(&((*workquery)[0])))
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (n), k, af, (m), tquery, tsize, q, (n), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, af, (m), tquery, tsize, df, (n), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, af, (m), tquery, tsize, df, (n), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, af, (m), tquery, tsize, cf, (m), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, af, (m), tquery, tsize, cf, (m), workquery, -1, info)
		(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
		ALLOCATE(t(tsize))
		ALLOCATE(work(lwork))
		(*srnamt) = *func() *[]byte {y := []byte("DGELQ"); return &y }()
		Dgelq((m), (n), af, (m), t, tsize, work, lwork, info)
		//
		//
		//     Generate the n-by-n matrix Q
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, one, q, (n))
		(*srnamt) = *func() *[]byte {y := []byte("DGEMLQ"); return &y }()
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (n), k, af, (m), t, tsize, q, (n), work, lwork, info)
		//
		//     Copy R
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), zero, zero, Lq, L)
		Dlacpy(func() *[]byte {y := []byte("Lower"); return &y }(), (m), (n), af, (m), Lq, L)
		//
		//     Compute |L - A*Q'| / |A| and store in result1
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), (n), -(*one), a, (m), q, (n), one, Lq, L)
		(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), a, (m), rwork))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), Lq, l, rwork))
		if (*anorm) > (*zero) {
			(*(result))[0] = (*resid) / ((*eps) * MAX(1, (*(n))) * (*anorm))
		} else {
			(*(result))[0] = (*zero)
		}
		//
		//     Compute |I - Q'*Q| and store in result(2)
		//
		Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, one, Lq, L)
		Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('C'); return &y }(), (n), (n), -(*one), q, (n), one, Lq, L)
		(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y := []byte("Upper"); return &y }(), (n), Lq, l, rwork))
		(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*(n))))
		//
		//     Generate random m-by-n matrix C and a copy cf
		//
		for (*j) = 1; (*j) <= (*(m)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (n), d(func() *int {y := 1; return &y }(), j))
		}
		(*dnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), d, (n), rwork))
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
		//
		//     Apply Q to C as Q*C
		//
		Dgemlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, af, (m), t, tsize, df, (n), work, lwork, info)
		//
		//     Compute |Q*D - Q*D| / |D|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (n), -(*one), q, (n), d, (n), one, df, (n))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
		if (*dnorm) > (*zero) {
			(*(result))[2] = (*resid) / ((*eps) * MAX(1, (*(n))) * (*dnorm))
		} else {
			(*(result))[2] = (*zero)
		}
		//
		//     Copy D into df again
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
		//
		//     Apply Q to D as QT*D
		//
		Dgemlq(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, af, (m), t, tsize, df, (n), work, lwork, info)
		//
		//     Compute |QT*D - QT*D| / |D|
		//
		Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (n), -(*one), q, (n), d, (n), one, df, (n))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
		if (*dnorm) > (*zero) {
			(*(result))[3] = (*resid) / ((*eps) * MAX(1, (*(n))) * (*dnorm))
		} else {
			(*(result))[3] = (*zero)
		}
		//
		//     Generate random n-by-m matrix D and a copy df
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), c(func() *int {y := 1; return &y }(), j))
		}
		(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), c, (m), rwork))
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
		//
		//     Apply Q to C as C*Q
		//
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, af, (m), t, tsize, cf, (m), work, lwork, info)
		//
		//     Compute |C*Q - C*Q| / |C|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (n), -(*one), c, (m), q, (n), one, cf, (m))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
		if (*cnorm) > (*zero) {
			(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*(n))) * (*cnorm))
		} else {
			(*(result))[4] = (*zero)
		}
		//
		//     Copy C into cf again
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
		//
		//     Apply Q to D as D*QT
		//
		Dgemlq(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, af, (m), t, tsize, cf, (m), work, lwork, info)
		//
		//     Compute |C*QT - C*QT| / |C|
		//
		Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), (n), -(*one), c, (m), q, (n), one, cf, (m))
		(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), cf, (m), rwork))
		if (*cnorm) > (*zero) {
			(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*(n))) * (*cnorm))
		} else {
			(*(result))[5] = (*zero)
		}
		//
	}
	//
	//     Deallocate all arrays
	//
	DEALLOCATE(a, af, q, r, rwork, work, t, c, d, cf, df)
	//
	return
}
