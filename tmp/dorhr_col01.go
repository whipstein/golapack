package goblas

import 

// DorhrCol01 tests DorhrCol using DLAtsQR, Dgemqrt and DORGtsQR.
// Therefore, DLAtsQR (part of DGEQR), Dgemqrt (part DGEMQR), DORGtsQR
// have to be tested before this test.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE DorhrCol01( m, n, mb1, nb1, nb2, result)
//
//       .. Scalar Arguments ..
//       intEGER           m, n, mb1, nb1, nb2
//       .. Return values ..
//       DOUBLE PRECISION  result(6)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// DorhrCol01 tests DorhrCol using DLAtsQR, Dgemqrt and DORGtsQR.
// Therefore, DLAtsQR (part of DGEQR), Dgemqrt (part DGEMQR), DORGtsQR
// have to be tested before this test.
//
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
// \param[in] N
// \verbatim
//          N is intEGER
//          Number of columns in test matrix.
// \endverbatim
// \param[in] mb1
// \verbatim
//          mb1 is intEGER
//          Number of row in row block in an input test matrix.
// \endverbatim
//
// \param[in] nb1
// \verbatim
//          nb1 is intEGER
//          Number of columns in column block an input test matrix.
// \endverbatim
//
// \param[in] nb2
// \verbatim
//          nb2 is intEGER
//          Number of columns in column block in an output test matrix.
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (6)
//          Results of each of the six tests below.
//          ( C is a M-by-N random matrix, D is a N-by-M random matrix)
//
//          result1 = | A - q * R | / (eps * m * |A|)
//          result(2) = | I - (Q**H) * Q | / (eps * m)
//          result(3) = | q * C - q * C | / (eps * m * |C|)
//          result(4) = | (Q**H) * C - (Q**H) * C | / (eps * m * |C|)
//          result(5) = | (d * Q) - d * Q | / (eps * m * |D|)
//          result(6) = | d * (Q**H) - d * (Q**H) | / (eps * m * |D|)
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
// \ingroup single_lin
//
//  =====================================================================
func DorhrCol01(m *int, n *int, mb1 *int, nb1 *int, nb2 *int, result *[]float64) {
	allocatable := new(float64)
	af := new(float64)
	q := new(float64)
	r := new(float64)
	rwork := new(float64)
	work := new(float64)
	t1 := new(float64)
	t2 := new(float64)
	diag := new(float64)
	c := new(float64)
	cf := new(float64)
	d := new(float64)
	df := new(float64)
	one := new(float64)
	zero := new(float64)
	testzeros := new(bool)
	info := new(int)
	i := new(int)
	j := new(int)
	k := new(int)
	l := new(int)
	lwork := new(int)
	nb1Ub := new(int)
	nb2Ub := new(int)
	nrb := new(int)
	anorm := new(float64)
	eps := new(float64)
	resid := new(float64)
	cnorm := new(float64)
	dnorm := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	workquery := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	_len := func() *[]byte {
		arr := make([]byte, -1)
		return &arr
	}()
	srnamt := new(int)
	common.SRmnAMC.srnamt = new(int)
	//
	//  -- lapACK test routine (version 3.9.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2019
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
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
	srnamt = common.SRmnAMC.srnamt
	//     ..
	//     .. Data statements ..
	(*iseed)[0], (*iseed)[1], (*iseed)[2], (*iseed)[3] = 1988, 1989, 1990, 1991
	//
	//     TEST MAtriCES WITH half OF MAtrix BEinG zeroS
	//
	(*testzeros) = false
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*k) = (Min((*(m)), (*(n))))
	(*l) = (*MAX((m), (n), func() *int {y := 1; return &y }()))
	//
	//     Dynamically allocate local arrays
	//
	ALLOCATE(a((m), (n)), af((m), (n)), q(L, L), r((m), L), rwork(l), c((m), (n)), cf((m), (n)), d((n), (m)), df((n), (m)))
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
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), a, (m), af, (m))
	//
	//     Number of row blocks in DLAtsQR
	//
	(*nRB) = (MAX(1, CEILING(DBLE((*(m))-(*(n)))/DBLE((*(mb1))-(*(n))))))
	//
	ALLOCATE(t1((nb1), (*(n))*(*nRB)))
	ALLOCATE(t2((nb2), (n)))
	ALLOCATE(diag((n)))
	//
	//     Begin determine lwork for the array work and allocate memory.
	//
	//     DLAtsQR requires nb1 to be bounded by N.
	//
	(*nb1Ub) = (Min((*(nb1)), (*(n))))
	//
	//     Dgemqrt requires nb2 to be bounded by N.
	//
	(*nb2Ub) = (Min((*(nb2)), (*(n))))
	//
	dlatsqr((m), (n), (mb1), nb1Ub, af, (m), t1, (nb1), workquery, -1, info)
	(*lwork) = (*int(&((*workquery)[0])))
	Dorgtsqr((m), (n), (mb1), (nb1), af, (m), t1, (nb1), workquery, -1, info)
	(*lwork) = (MAX((*lwork), int(&((*workquery)[0]))))
	//
	//     In Dgemqrt, work is N*nb2Ub if side = 'L',
	//                or  M*nb2Ub if side = 'R'.
	//
	(*lwork) = (*MAX(lwork, (*nb2Ub)*(*(n)), (*nb2Ub)*(*(m))))
	//
	ALLOCATE(work(lwork))
	//
	//     End allocate memory for work.
	//
	//
	//     Begin Householder re_construction routines
	//
	//     Factor the matrix A in the array af.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DLAtsQR"); return &y }()
	dlatsqr((m), (n), (mb1), nb1Ub, af, (m), t1, (nb1), work, lwork, info)
	//
	//     Copy the factor R into the array R.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dlacpy"); return &y }()
	Dlacpy(func() *byte {y := byte('U'); return &y }(), (n), (n), af, (m), r, (m))
	//
	//     Re_construct the orthogonal matrix Q.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORGtsQR"); return &y }()
	Dorgtsqr((m), (n), (mb1), (nb1), af, (m), t1, (nb1), work, lwork, info)
	//
	//     Perform the Householder re_construction, the result is stored
	//     the arrays af and t2.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DorhrCol"); return &y }()
	DorhrCol((m), (n), (nb2), af, (m), t2, (nb2), diag, info)
	//
	//     Compute the factor R_hr corresponding to the Householder
	//     re_constructed Q_hr and place it in the upper triangle of af to
	//     match the Q storage format in Dgeqrt. R_hr = R_tsqr * s,
	//     this means changing the sign of I-th row of the matrix R_tsqr
	//     according to sign of of I-th diagonal element diag(i) of the
	//     matrix S.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dlacpy"); return &y }()
	Dlacpy(func() *byte {y := byte('U'); return &y }(), (n), (n), r, (m), af, (m))
	//
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		if (*diag(i)) == -(*one) {
			Dscal((*(n))+1-(*i), -(*one), af(I, I), (m))
		}
	}
	//
	//     End Householder re_construction routines.
	//
	//
	//     Generate the m-by-m matrix Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), zero, one, q, (m))
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y }()
	Dgemqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (m), k, nb2Ub, af, (m), t2, (nb2), q, (m), work, info)
	//
	//     Copy R
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, r, (m))
	//
	Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), (n), af, (m), r, (m))
	//
	//     TEST 1
	//     Compute |R - (Q**t)*a| / ( eps * m * |A|) and store in result1
	//
	Dgemm(func() *byte {y := byte('T'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), -(*one), q, (m), a, (m), one, r, (m))
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), a, (m), rwork))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), r, (m), rwork))
	if (*anorm) > (*zero) {
		(*(result))[0] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*anorm))
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     TEST 2
	//     Compute |I - (Q**t)*Q| / ( eps * m) and store in result(2)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), zero, one, r, (m))
	Dsyrk(func() *byte {y := byte('U'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (m), -(*one), q, (m), one, r, (m))
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), r, (m), rwork))
	(*(result))[1] = (*resid) / ((*eps) * MAX(1, (*(m))))
	//
	//     Generate random m-by-n matrix C
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), iseed, (m), c(func() *int {y := 1; return &y }(), j))
	}
	(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), c, (m), rwork))
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), c, (m), cf, (m))
	//
	//     Apply Q to C as Q*C = cf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y }()
	Dgemqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), k, nb2Ub, af, (m), t2, (nb2), cf, (m), work, info)
	//
	//     TEST 3
	//     Compute |cf - Q*C| / ( eps *  m * |C|)
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
	//     Apply Q to C as (Q**t)*C = cf
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y }()
	Dgemqrt(func() *byte {y := byte('L'); return &y }(), func() *byte {y := byte('T'); return &y }(), (m), (n), k, nb2Ub, af, (m), t2, (nb2), cf, (m), work, info)
	//
	//     TEST 4
	//     Compute |cf - (Q**t)*C| / ( eps * m * |C|)
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
	//     Apply Q to D as D*Q = df
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y }()
	Dgemqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('N'); return &y }(), (n), (m), k, nb2Ub, af, (m), t2, (nb2), df, (n), work, info)
	//
	//     TEST 5
	//     Compute |df - D*Q| / ( eps * m * |D|)
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
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (m), d, (n), df, (n))
	//
	//     Apply Q to D as D*QT = df
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgemqrt"); return &y }()
	Dgemqrt(func() *byte {y := byte('R'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), k, nb2Ub, af, (m), t2, (nb2), df, (n), work, info)
	//
	//     TEST 6
	//     Compute |df - D*(Q**t)| / ( eps * m * |D|)
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('T'); return &y }(), (n), (m), (m), -(*one), d, (n), q, (m), one, df, (n))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (m), df, (n), rwork))
	if (*dnorm) > (*zero) {
		(*(result))[5] = (*resid) / ((*eps) * MAX(1, (*(m))) * (*dnorm))
	} else {
		(*(result))[5] = (*zero)
	}
	//
	//     Deallocate all arrays
	//
	DEALLOCATE(a, af, q, r, rwork, work, t1, t2, diag, c, d, cf, df)
	//
	return
	//
	//     End of DorhrCol01
	//
}
