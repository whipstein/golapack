package goblas

import 

// Dqrt05 tests DTPLQT and DTPMLQT.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqt05(m,N,L,nb,result)
//
//       .. Scalar Arguments ..
//       intEGER lwork, m, n, l, nb, ldt
//       .. Return values ..
//       DOUBLE PRECISION result(6)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt05 tests DTPLQT and DTPMLQT.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          Number of rows in lower part of the test matrix.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          Number of columns in test matrix.
// \endverbatim
//
// \param[in] L
// \verbatim
//          L is intEGER
//          The number of rows of the upper trapezoidal part the
//          lower test matrix.  0 <= L <= M.
// \endverbatim
//
// \param[in] nb
// \verbatim
//          nb is intEGER
//          Block size of test matrix.  nb <= N.
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (6)
//          Results of each of the six tests below.
//
//          result(1) = | A - Q R |
//          result(2) = | I - Q^H Q |
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
func Dqt05(m *int, n *int, L *int, nb *int, result *[]float64) {
	lwork := new(int)
	ldt := new(int)
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
	one := new(float64)
	zero := new(float64)
	info := new(int)
	j := new(int)
	k := new(int)
	n2 := new(int)
	np1 := new(int)
	i := new(int)
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
	//     .. Data statements ..
	(*iseed)[0], (*iseed)[1], (*iseed)[2], (*iseed)[3] = 1988, 1989, 1990, 1991
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*k) = (*(m))
	(*n2) = (*(m)) + (*(n))
	if (*(n)) > 0 {
		(*nP1) = (*(m)) + 1
	} else {
		(*nP1) = 1
	}
	(*lwork) = (*n2) * (*n2) * (*(nb))
	//
	//     Dynamically allocate all arrays
	//
	ALLOCATE(a((m), n2), af((m), n2), q(n2, n2), r(n2, n2), rwork(n2), work(lwork), t((nb), (m)), c(n2, (m)), cf(n2, (m)), d((m), n2), df((m), n2))
	//
	//     Put random stuff into A
	//
	(*lDT) = (*(nb))
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), n2, zero, zero, a, (m))
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (nb), (m), zero, zero, t, (nb))
	for (*j) = 1; (*j) <= (*(m)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (*(m))-(*j)+1, a(J, j))
	}
	if (*(n)) > 0 {
		for (*j) = 1; (*j) <= (*(n))-(*(l)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), a(func() *int {y := 1; return &y }(), Min((*(n))+(*(m)), (*(m))+1)+(*j)-1))
		}
	}
	if (*(l)) > 0 {
		for (*j) = 1; (*j) <= (*(l)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, (*(m))-(*j)+1, a(J, Min((*(n))+(*(m)), (*(n))+(*(m))-(*(l))+1)+(*j)-1))
		}
	}
	//
	//     Copy the matrix A to the array af.
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), n2, a, (m), af, (m))
	//
	//     Factor the matrix A in the array af.
	//
	Dtplqt((m), (n), (l), (nb), af, (m), af(func() *int {y := 1; return &y }(), np1), (m), t, ldt, work, info)
	//
	//     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, n2, zero, one, q, n2)
	Dgemlqt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), n2, n2, k, (nb), af, (m), t, ldt, q, n2, work, info)
	//
	//     Copy L
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, n2, zero, zero, r, n2)
	Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (m), n2, af, (m), r, n2)
	//
	//     Compute |L - A*Q*T| / |A| and store in result(1)
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), n2, n2, -(*one), a, (m), q, n2, one, r, n2)
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), n2, a, (m), rwork))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), n2, r, n2, rwork))
	if (*anorm) > (*zero) {
		(*(result))[0] = (*resid) / ((*eps) * (*anorm) * MAX(1, (*n2)))
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute |I - Q*Q'| and store in result(2)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, n2, zero, one, r, n2)
	Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('N'); return &y }(), n2, n2, -(*one), q, n2, one, r, n2)
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), n2, r, n2, rwork))
	(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*n2)))
	//
	//     Generate random m-by-n matrix C and a copy cf
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, (m), zero, one, c, n2)
	for (*j) = 1; (*j) <= (*(m)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, n2, c(func() *int {y := 1; return &y }(), j))
	}
	(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), n2, (m), c, n2, rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, (m), c, n2, cf, n2)
	//
	//     Apply Q to C as Q*C
	//
	Dtpmlqt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, (l), (nb), af(func() *int {y := 1; return &y }(), np1), (m), t, ldt, cf, n2, cf(np1, func() *int {y := 1; return &y }()), n2, work, info)
	//
	//     Compute |Q*C - Q*C| / |C|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), n2, (m), n2, -(*one), q, n2, c, n2, one, cf, n2)
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), n2, (m), cf, n2, rwork))
	if (*cnorm) > (*zero) {
		(*(result))[2] = (*resid) / ((*eps) * MAX(1, (*n2)) * (*cnorm))
	} else {
		(*(result))[2] = (*zero)
	}
	//
	//     Copy C into cf again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n2, (m), c, n2, cf, n2)
	//
	//     Apply Q to C as QT*C
	//
	Dtpmlqt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, (l), (nb), af(func() *int {y := 1; return &y }(), np1), (m), t, ldt, cf, n2, cf(np1, func() *int {y := 1; return &y }()), n2, work, info)
	//
	//     Compute |QT*C - QT*C| / |C|
	//
	Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), n2, (m), n2, -(*one), q, n2, c, n2, one, cf, n2)
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), n2, (m), cf, n2, rwork))
	if (*cnorm) > (*zero) {
		(*(result))[3] = (*resid) / ((*eps) * MAX(1, (*n2)) * (*cnorm))
	} else {
		(*(result))[3] = (*zero)
	}
	//
	//     Generate random m-by-n matrix D and a copy df
	//
	for (*j) = 1; (*j) <= (*n2); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), d(func() *int {y := 1; return &y }(), j))
	}
	(*dnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), n2, d, (m), rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), n2, d, (m), df, (m))
	//
	//     Apply Q to D as D*Q
	//
	Dtpmlqt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, (l), (nb), af(func() *int {y := 1; return &y }(), np1), (m), t, ldt, df, (m), df(func() *int {y := 1; return &y }(), np1), (m), work, info)
	//
	//     Compute |D*Q - D*Q| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), n2, n2, -(*one), d, (m), q, n2, one, df, (m))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), n2, df, (m), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*n2)) * (*dnorm))
	} else {
		(*(result))[4] = (*zero)
	}
	//
	//     Copy D into df again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), n2, d, (m), df, (m))
	//
	//     Apply Q to D as D*QT
	//
	Dtpmlqt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, (l), (nb), af(func() *int {y := 1; return &y }(), np1), (m), t, ldt, df, (m), df(func() *int {y := 1; return &y }(), np1), (m), work, info)
	//
	//     Compute |D*QT - D*QT| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), n2, n2, -(*one), d, (m), q, n2, one, df, (m))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), n2, df, (m), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*n2)) * (*dnorm))
	} else {
		(*(result))[5] = (*zero)
	}
	//
	//     Deallocate all arrays
	//
	DEALLOCATE(a, af, q, r, rwork, work, t, c, d, cf, df)
	return
}
