package goblas

import 

// Dqrt04 tests Dgeqrt and Dgemqrt.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqrt04(m,N,nb,result)
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
// Dqrt04 tests Dgeqrt and Dgemqrt.
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
func Dqrt04(m *int, n *int, nb *int, result *[]float64) {
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
	l := new(int)
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
	(*l) = (MAX((*(m)), (*(n))))
	(*lwork) = MAX(2, (*l)) * MAX(2, (*l)) * (*(nb))
	//
	//     Dynamically allocate local arrays
	//
	ALLOCATE(a((m), (n)), af((m), (n)), q((m), (m)), r((m), L), rwork(l), work(lwork), t((nb), (n)), c((m), (n)), cf((m), (n)), d((n), (m)), df((n), (m)))
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
	Dgeqrt((m), (n), (nb), af, (m), t, ldt, work, info)
	//
	//     Generate the m-by-m matrix Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), zero, one, q, (m))
	Dgemqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (m), k, (nb), af, (m), t, ldt, q, (m), work, info)
	//
	//     Copy R
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, r, (m))
	Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), (n), af, (m), r, (m))
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
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), zero, one, r, (m))
	Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('C'); return &y }(), (m), (m), -(*one), q, (m), one, r, (m))
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), r, (m), rwork))
	(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*(m))))
	//
	//     Generate random m-by-n matrix C and a copy cf
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), c(func() *int {y := 1; return &y }(), j))
	}
	(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), c, (m), rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
	//
	//     Apply Q to C as Q*C
	//
	Dgemqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, (nb), af, (m), t, (nb), cf, (m), work, info)
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
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
	//
	//     Apply Q to C as QT*C
	//
	Dgemqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, (nb), af, (m), t, (nb), cf, (m), work, info)
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
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
	//
	//     Apply Q to D as D*Q
	//
	Dgemqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, (nb), af, (m), t, (nb), df, (n), work, info)
	//
	//     Compute |D*Q - D*Q| / |D|
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), (m), -(*one), d, (n), q, (m), one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
	if (*cnorm) > (*zero) {
		(*(result))[4] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
	} else {
		(*(result))[4] = (*zero)
	}
	//
	//     Copy D into df again
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
	//
	//     Apply Q to D as D*QT
	//
	Dgemqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, (nb), af, (m), t, (nb), df, (n), work, info)
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
	//     Deallocate all arrays
	//
	DEALLOCATE(a, af, q, r, rwork, work, t, c, d, cf, df)
	//
	return
}
