package goblas

import 

// Dlqt04 tests DGELQT and DGEMLQT.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlqt04(m,N,nb,result)
//
//       .. Scalar Arguments ..
//       intEGER m, n, nb, ldt
//       .. Return values ..
//       DOUBLE PRECISION result(6)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlqt04 tests DGELQT and DGEMLQT.
// \endverbatim
//
//  Arguments:
//  ==========
//
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
//
// \param[in] nb
// \verbatim
//          nb is intEGER
//          Block size of test matrix.  nb <= Min(m,N).
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (6)
//          Results of each of the six tests below.
//
//          result1 = | A - L Q |
//          result(2) = | I - Q Q^H |
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
func Dlqt04(m *int, n *int, nb *int, result *[]float64) {
	ldt := new(int)
	allocatable := new(float64)
	q := new(float64)
	l := new(float64)
	rwork := new(float64)
	work := new(float64)
	t := new(float64)
	cf := new(float64)
	df := new(float64)
	a := new(float64)
	c := new(float64)
	d := new(float64)
	one := new(float64)
	zero := new(float64)
	info := new(int)
	j := new(int)
	k := new(int)
	Ll := new(int)
	lwork := new(int)
	anorm := new(float64)
	eps := new(float64)
	resid := new(float64)
	cnorm := new(float64)
	dnorm := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
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
	//     ..
	//     .. Data statements ..
	(*iseed)[0], (*iseed)[1], (*iseed)[2], (*iseed)[3] = 1988, 1989, 1990, 1991
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*k) = (Min((*(m)), (*(n))))
	(*lL) = (MAX((*(m)), (*(n))))
	(*lwork) = MAX(2, (*lL)) * MAX(2, (*lL)) * (*(nb))
	//
	//     Dynamically allocate local arrays
	//
	ALLOCATE(a((m), (n)), af((m), (n)), q((n), (n)), L(LL, (n)), rwork(LL), work(lwork), t((nb), (n)), c((m), (n)), cf((m), (n)), d((n), (m)), df((n), (m)))
	//
	//     Put random numbers into A and copy to af
	//
	(*lDT) = (*(nb))
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), a(func() *int {y := 1; return &y }(), j))
	}
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), a, (m), af, (m))
	//
	//     Factor the matrix A in the array af.
	//
	Dgelqt((m), (n), (nb), af, (m), t, ldt, work, info)
	//
	//     Generate the n-by-n matrix Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, q, (n))
	Dgemlqt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (n), k, (nb), af, (m), t, ldt, q, (n), work, info)
	//
	//     Copy R
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, l, LL)
	Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (m), (n), af, (m), l, LL)
	//
	//     Compute |L - A*Q'| / |A| and store in result1
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), (n), -(*one), a, (m), q, (n), one, l, LL)
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), a, (m), rwork))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), l, LL, rwork))
	if (*anorm) > (*zero) {
		(*(result))[0] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*anorm))
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute |I - Q'*Q| and store in result(2)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, l, LL)
	Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('C'); return &y }(), (n), (n), -(*one), q, (n), one, l, LL)
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (n), l, LL, rwork))
	(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*(n))))
	//
	//     Generate random m-by-n matrix C and a copy cf
	//
	for (*j) = 1; (*j) <= (*(m)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (n), d(func() *int {y := 1; return &y }(), j))
	}
	(*dnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), d, (n), rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
	//
	//     Apply Q to C as Q*C
	//
	Dgemlqt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, (nb), af, (m), t, (nb), df, (n), work, info)
	//
	//     Compute |Q*D - Q*D| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (n), -(*one), q, (n), d, (n), one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
	if (*dnorm) > (*zero) {
		(*(result))[2] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
	} else {
		(*(result))[2] = (*zero)
	}
	//
	//     Copy D into df again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
	//
	//     Apply Q to D as QT*D
	//
	Dgemlqt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, (nb), af, (m), t, (nb), df, (n), work, info)
	//
	//     Compute |QT*D - QT*D| / |D|
	//
	Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (n), -(*one), q, (n), d, (n), one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
	if (*dnorm) > (*zero) {
		(*(result))[3] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
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
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
	//
	//     Apply Q to C as C*Q
	//
	Dgemlqt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, (nb), af, (m), t, (nb), cf, (m), work, info)
	//
	//     Compute |C*Q - C*Q| / |C|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (n), -(*one), c, (m), q, (n), one, cf, (m))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
	} else {
		(*(result))[4] = (*zero)
	}
	//
	//     Copy C into cf again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
	//
	//     Apply Q to D as D*QT
	//
	Dgemlqt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, (nb), af, (m), t, (nb), cf, (m), work, info)
	//
	//     Compute |C*QT - C*QT| / |C|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), (n), -(*one), c, (m), q, (n), one, cf, (m))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), cf, (m), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
	} else {
		(*(result))[5] = (*zero)
	}
	//
	//     Deallocate all arrays
	//
	DEALLOCATE(a, af, q, l, rwork, work, t, c, d, cf, df)
	//
	return
}
