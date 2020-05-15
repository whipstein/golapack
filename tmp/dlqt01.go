package goblas

import 

// Dlqt01 tests Dgelqf, which computes the LQ factorization of an m-by-n
// matrix a, and partially tests DORGLQ which forms the n-by-n
// orthogonal matrix Q.
//
// Dlqt01 compares L with A*Q', and checks that Q is orthogonal.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlqt01( m, n, a, af, q, l, lda, tau, work, lwork,
//                          rwork, result)
//
//       .. Scalar Arguments ..
//       intEGER            lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), af( lda, *), L( lda, *),
//      $                   q( lda, *), result(*), rwork(*), tau(*),
//      $                   work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlqt01 tests Dgelqf, which computes the LQ factorization of an m-by-n
// matrix a, and partially tests DORGLQ which forms the n-by-n
// orthogonal matrix Q.
//
// Dlqt01 compares L with A*Q', and checks that Q is orthogonal.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix A.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The m-by-n matrix A.
// \endverbatim
//
// \param[out] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the LQ factorization of a, as returned by Dgelqf.
//          See Dgelqf for further details.
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,N)
//          The n-by-n orthogonal matrix Q.
// \endverbatim
//
// \param[out] L
// \verbatim
//          L is DOUBLE PRECISION array, dimension (lda,max(m,N))
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the arrays a, af, Q and L.
//          lda >= max(m,N).
// \endverbatim
//
// \param[out] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (min(m,N))
//          The scalar factors of the elementary reflectors, as returned
//          by Dgelqf.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (lwork)
// \endverbatim
//
// \param[in] lwork
// \verbatim
//          lwork is intEGER
//          The dimension of the array work.
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (max(m,N))
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (2)
//          The test ratios:
//          result1 = norm( L - A*Q') / ( N * norm(a) * eps)
//          result(2) = norm( I - Q*Q') / ( N * eps)
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
func Dlqt01(m *int, n *int, a *[][]float64, af *[][]float64, q *[][]float64, l *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
	zero := new(float64)
	one := new(float64)
	rogue := new(float64)
	info := new(int)
	minmn := new(int)
	anorm := new(float64)
	eps := new(float64)
	resid := new(float64)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	common.srnamc.srnamt = new(int)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	(*rogue) = -1.0e+10
	//     ..
	//     .. Local Scalars ..
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
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Executable Statements ..
	//
	(*minmn) = (Min((*(m)), (*(n))))
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     Copy the matrix A to the array af.
	//
	Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), (a), (lda), (af), (lda))
	//
	//     Factorize the matrix A in the array af.
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgelqf"); return &y }()
	Dgelqf((m), (n), (af), (lda), (tau), (work), (lwork), info)
	//
	//     Copy details of Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), rogue, rogue, (q), (lda))
	if (*(n)) > 1 {
		Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), (*(n))-1, &((*(af))[0][1]), (lda), &((*(q))[0][1]), (lda))
	}
	//
	//     Generate the n-by-n matrix Q
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("DORGLQ"); return &y }()
	Dorglq((n), (n), minmn, (q), (lda), (tau), (work), (lwork), info)
	//
	//     Copy L
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, (l), (lda))
	Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (m), (n), (af), (lda), (l), (lda))
	//
	//     Compute L - A*Q'
	//
	Dgemm(func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (m), (n), (n), -(*one), (a), (lda), (q), (lda), one, (l), (lda))
	//
	//     Compute norm( L - Q'*a) / ( N * norm(a) * eps) .
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (a), (lda), (rwork)))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (l), (lda), (rwork)))
	if (*anorm) > (*zero) {
		(*(result))[0] = (((*resid) / DBLE(MAX(1, (*(n))))) / (*anorm)) / (*eps)
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute I - Q*Q'
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, (l), (lda))
	Dsyrk(func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), (n), (n), -(*one), (q), (lda), one, (l), (lda))
	//
	//     Compute norm( I - Q*Q') / ( N * eps) .
	//
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (n), (l), (lda), (rwork)))
	//
	(*(result))[1] = ((*resid) / DBLE(MAX(1, (*(n))))) / (*eps)
	//
	return
	//
	//     End of Dlqt01
	//
}
