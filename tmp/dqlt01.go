package goblas

import 

// Dqlt01 tests Dgeqlf, which computes the QL factorization of an m-by-n
// matrix a, and partially tests Dorgql which forms the m-by-m
// orthogonal matrix Q.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqlt01( m, n, a, af, q, l, lda, tau, work, lwork,
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
// Dqlt01 tests Dgeqlf, which computes the QL factorization of an m-by-n
// matrix a, and partially tests Dorgql which forms the m-by-m
// orthogonal matrix Q.
//
// Dqlt01 compares L with Q'*a, and checks that Q is orthogonal.
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
//          Details of the QL factorization of a, as returned by Dgeqlf.
//          See Dgeqlf for further details.
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,M)
//          The m-by-m orthogonal matrix Q.
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
//          The leading dimension of the arrays a, af, Q and R.
//          lda >= max(m,N).
// \endverbatim
//
// \param[out] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (min(m,N))
//          The scalar factors of the elementary reflectors, as returned
//          by Dgeqlf.
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
//          rwork is DOUBLE PRECISION array, dimension (m)
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (2)
//          The test ratios:
//          result1 = norm( L - Q'*a) / ( m * norm(a) * eps)
//          result(2) = norm( I - Q'*Q) / ( m * eps)
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
func Dqlt01(m *int, n *int, a *[][]float64, af *[][]float64, q *[][]float64, l *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
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
	(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqlf"); return &y }()
	Dgeqlf((m), (n), (af), (lda), (tau), (work), (lwork), info)
	//
	//     Copy details of Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), rogue, rogue, (q), (lda))
	if (*(m)) >= (*(n)) {
		if (*(n)) < (*(m)) && (*(n)) > 0 {
			Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (*(m))-(*(n)), (n), (af), (lda), &((*(q))[0][(*(m))-(*(n))+0]), (lda))
		}
		if (*(n)) > 1 {
			Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (*(n))-1, (*(n))-1, &((*(af))[(*(m))-(*(n))+0][1]), (lda), &((*(q))[(*(m))-(*(n))+0][(*(m))-(*(n))+1]), (lda))
		}
	} else {
		if (*(m)) > 1 {
			Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (*(m))-1, (*(m))-1, &((*(af))[0][(*(n))-(*(m))+1]), (lda), &((*(q))[0][1]), (lda))
		}
	}
	//
	//     Generate the m-by-m matrix Q
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgql"); return &y }()
	Dorgql((m), (m), minmn, (q), (lda), (tau), (work), (lwork), info)
	//
	//     Copy L
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, (l), (lda))
	if (*(m)) >= (*(n)) {
		if (*(n)) > 0 {
			Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (n), (n), &((*(af))[(*(m))-(*(n))+0][0]), (lda), &((*(l))[(*(m))-(*(n))+0][0]), (lda))
		}
	} else {
		if (*(n)) > (*(m)) && (*(m)) > 0 {
			Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (*(n))-(*(m)), (af), (lda), (l), (lda))
		}
		if (*(m)) > 0 {
			Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (m), (m), &((*(af))[0][(*(n))-(*(m))+0]), (lda), &((*(l))[0][(*(n))-(*(m))+0]), (lda))
		}
	}
	//
	//     Compute L - Q'*a
	//
	Dgemm(func() *[]byte {y :=[]byte("Transpose"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), (m), (n), (m), -(*one), (q), (lda), (a), (lda), one, (l), (lda))
	//
	//     Compute norm( L - Q'*a) / ( m * norm(a) * eps) .
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (a), (lda), (rwork)))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (l), (lda), (rwork)))
	if (*anorm) > (*zero) {
		(*(result))[0] = (((*resid) / DBLE(MAX(1, (*(m))))) / (*anorm)) / (*eps)
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute I - Q'*Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (m), zero, one, (l), (lda))
	Dsyrk(func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (m), (m), -(*one), (q), (lda), one, (l), (lda))
	//
	//     Compute norm( I - Q'*Q) / ( m * eps) .
	//
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (m), (l), (lda), (rwork)))
	//
	(*(result))[1] = ((*resid) / DBLE(MAX(1, (*(m))))) / (*eps)
	//
	return
	//
	//     End of Dqlt01
	//
}
