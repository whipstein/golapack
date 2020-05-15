package goblas

import 

// Dqrt05 tests DTPQRT and Dtpmqrt.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqrt05(m,N,L,nb,result)
//
//       .. Scalar Arguments ..
//       inTEGER lwork, m, n, L, nb, ldt
//       .. Return values ..
//       DOUBLE PRECISION result(6)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt05 tests DTPQRT and Dtpmqrt.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          Number of rows in lower part of the test matrix.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          Number of columns in test matrix.
// \endverbatim
//
// \param[in] L
// \verbatim
//          L is inTEGER
//          The number of rows of the upper trapezoidal part the
//          lower test matrix.  0 <= L <= M.
// \endverbatim
//
// \param[in] nb
// \verbatim
//          nb is inTEGER
//          Block size of test matrix.  nb <= N.
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (6)
//          Results of each of the six tests below.
//
//          result1 = | A - Q R |
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
func Dqrt05(m *int, n *int, L *int, nb *int, result *[]float64) {
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
	m2 := new(int)
	np1 := new(int)
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
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	(*k) = (*(n))
	(*m2) = (*(m)) + (*(n))
	if (*(m)) > 0 {
		(*np1) = (*(n)) + 1
	} else {
		(*np1) = 1
	}
	(*lwork) = (*m2) * (*m2) * (*(nb))
	//
	//     Dynamically allocate all arrays
	//
	ALLOCATE(a(m2, (n)), af(m2, (n)), q(m2, m2), r(m2, m2), rwork(m2), work(lwork), t((nb), (n)), c(m2, (n)), cf(m2, (n)), d((n), m2), df((n), m2))
	//
	//     Put random stuff into A
	//
	(*ldt) = (*(nb))
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, (n), zero, zero, a, m2)
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (nb), (n), zero, zero, t, (nb))
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y}(), iseed, j, a(func() *int {y := 1; return &y}(), j))
	}
	if (*(m)) > 0 {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y}(), iseed, (*(m))-(*(l)), a(Min((*(n))+(*(m)), (*(n))+1), j))
		}
	}
	if (*(l)) > 0 {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y}(), iseed, Min((*j), (*(l))), a(Min((*(n))+(*(m)), (*(n))+(*(m))-(*(l))+1), j))
		}
	}
	//
	//     Copy the matrix A to the array af.
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, (n), a, m2, af, m2)
	//
	//     Factor the matrix A in the array af.
	//
	Dtpqrt((m), (n), (l), (nb), af, m2, af(np1, func() *int {y := 1; return &y}()), m2, t, ldt, work, info)
	//
	//     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, m2, zero, one, q, m2)
	Dgemqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), m2, m2, k, (nb), af, m2, t, ldt, q, m2, work, info)
	//
	//     Copy R
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, (n), zero, zero, r, m2)
	Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y}(), m2, (n), af, m2, r, m2)
	//
	//     Compute |R - Q'*A| / |A| and store in result1
	//
	Dgemm(func() *byte {y := byte('T'); return &y}(), func() *byte {y := byte('N'); return &y}(), m2, (n), m2, -(*one), q, m2, a, m2, one, r, m2)
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), m2, (n), a, m2, rwork))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), m2, (n), r, m2, rwork))
	if (*anorm) > (*zero) {
		(*(result))[0] = (*resid) / ((*eps) * (*anorm) * MAX(1, (*m2)))
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute |I - Q'*Q| and store in result(2)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, m2, zero, one, r, m2)
	Dsyrk(func() *byte {y := byte('U'); return &y}(), func() *byte {y := byte('C'); return &y}(), m2, m2, -(*one), q, m2, one, r, m2)
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y}(), func() *[]byte {y :=[]byte("Upper"); return &y}(), m2, r, m2, rwork))
	(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*m2)))
	//
	//     Generate random m-by-n matrix C and a copy cf
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y}(), iseed, m2, c(func() *int {y := 1; return &y}(), j))
	}
	(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), m2, (n), c, m2, rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, (n), c, m2, cf, m2)
	//
	//     Apply Q to C as Q*C
	//
	Dtpmqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('N'); return &y}(), (m), (n), k, (l), (nb), af(np1, func() *int {y := 1; return &y}()), m2, t, ldt, cf, m2, cf(np1, func() *int {y := 1; return &y}()), m2, work, info)
	//
	//     Compute |Q*C - Q*C| / |C|
	//
	Dgemm(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), m2, (n), m2, -(*one), q, m2, c, m2, one, cf, m2)
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), m2, (n), cf, m2, rwork))
	if (*cnorm) > (*zero) {
		(*(result))[2] = (*resid) / ((*eps) * MAX(1, (*m2)) * (*cnorm))
	} else {
		(*(result))[2] = (*zero)
	}
	//
	//     Copy C into cf again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), m2, (n), c, m2, cf, m2)
	//
	//     Apply Q to C as QT*C
	//
	Dtpmqrt(func() *byte {y := byte('L'); return &y}(), func() *byte {y := byte('T'); return &y}(), (m), (n), k, (l), (nb), af(np1, func() *int {y := 1; return &y}()), m2, t, ldt, cf, m2, cf(np1, func() *int {y := 1; return &y}()), m2, work, info)
	//
	//     Compute |QT*C - QT*C| / |C|
	//
	Dgemm(func() *byte {y := byte('T'); return &y}(), func() *byte {y := byte('N'); return &y}(), m2, (n), m2, -(*one), q, m2, c, m2, one, cf, m2)
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), m2, (n), cf, m2, rwork))
	if (*cnorm) > (*zero) {
		(*(result))[3] = (*resid) / ((*eps) * MAX(1, (*m2)) * (*cnorm))
	} else {
		(*(result))[3] = (*zero)
	}
	//
	//     Generate random n-by-m matrix D and a copy df
	//
	for (*j) = 1; (*j) <= (*m2); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y}(), iseed, (n), d(func() *int {y := 1; return &y}(), j))
	}
	(*dnorm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), m2, d, (n), rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), (n), m2, d, (n), df, (n))
	//
	//     Apply Q to D as D*Q
	//
	Dtpmqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('N'); return &y}(), (n), (m), (n), (l), (nb), af(np1, func() *int {y := 1; return &y}()), m2, t, ldt, df, (n), df(func() *int {y := 1; return &y}(), np1), (n), work, info)
	//
	//     Compute |D*Q - D*Q| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), (n), m2, m2, -(*one), d, (n), q, m2, one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), m2, df, (n), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*m2)) * (*dnorm))
	} else {
		(*(result))[4] = (*zero)
	}
	//
	//     Copy D into df again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), (n), m2, d, (n), df, (n))
	//
	//     Apply Q to D as D*QT
	//
	Dtpmqrt(func() *byte {y := byte('R'); return &y}(), func() *byte {y := byte('T'); return &y}(), (n), (m), (n), (l), (nb), af(np1, func() *int {y := 1; return &y}()), m2, t, ldt, df, (n), df(func() *int {y := 1; return &y}(), np1), (n), work, info)
	//
	//     Compute |D*QT - D*QT| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('T'); return &y}(), (n), m2, m2, -(*one), d, (n), q, m2, one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (n), m2, df, (n), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*m2)) * (*dnorm))
	} else {
		(*(result))[5] = (*zero)
	}
	//
	//     Deallocate all arrays
	//
	DEALLOCATE(A, af, q, r, rwork, work, t, c, d, cf, df)
	return
}
