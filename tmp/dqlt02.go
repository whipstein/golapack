package goblas

import 

// Dqlt02 tests Dorgql, which generates an m-by-n matrix Q with
// orthonornmal columns that is defined as the product of k elementary
// reflectors.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqlt02( m, n, k, a, af, q, l, lda, tau, work, lwork,
//                          rwork, result)
//
//       .. Scalar Arguments ..
//       intEGER            k, lda, lwork, m, N
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
// Dqlt02 tests Dorgql, which generates an m-by-n matrix Q with
// orthonornmal columns that is defined as the product of k elementary
// reflectors.
//
// Given the QL factorization of an m-by-n matrix a, Dqlt02 generates
// the orthogonal matrix Q defined by the factorization of the last k
// columns of A; it compares L(m-n+1:m,n-k+1:n) with
// q(1:m,m-n+1:m)'*a(1:m,n-k+1:n), and checks that the columns of Q are
// orthonormal.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix Q to be generated.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix Q to be generated.
//          M >= N >= 0.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is intEGER
//          The number of elementary reflectors whose product defines the
//          matrix Q. N >= K >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The m-by-n matrix A which was factorized by Dqlt01.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the QL factorization of a, as returned by Dgeqlf.
//          See Dgeqlf for further details.
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[out] L
// \verbatim
//          L is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the arrays a, af, Q and L. lda >= M.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (n)
//          The scalar factors of the elementary reflectors corresponding
//          to the QL factorization in af.
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
func Dqlt02(m *int, n *int, k *int, a *[][]float64, af *[][]float64, q *[][]float64, l *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
	zero := new(float64)
	one := new(float64)
	rogue := new(float64)
	info := new(int)
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
	//     Quick return if possible
	//
	if (*(m)) == 0 || (*(n)) == 0 || (*(k)) == 0 {
		(*(result))[0] = (*zero)
		(*(result))[1] = (*zero)
		return
	}
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     Copy the last k columns of the factorization to the array Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), rogue, rogue, (q), (lda))
	if (*(k)) < (*(m)) {
		Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (*(m))-(*(k)), (k), &((*(af))[0][(*(n))-(*(k))+0]), (lda), &((*(q))[0][(*(n))-(*(k))+0]), (lda))
	}
	if (*(k)) > 1 {
		Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (*(k))-1, (*(k))-1, &((*(af))[(*(m))-(*(k))+0][(*(n))-(*(k))+1]), (lda), &((*(q))[(*(m))-(*(k))+0][(*(n))-(*(k))+1]), (lda))
	}
	//
	//     Generate the last n columns of the matrix Q
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgql"); return &y }()
	Dorgql((m), (n), (k), (q), (lda), &((*(tau))[(*(n))-(*(k))+0]), (work), (lwork), info)
	//
	//     Copy L(m-n+1:m,n-k+1:n)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (k), zero, zero, &((*(l))[(*(m))-(*(n))+0][(*(n))-(*(k))+0]), (lda))
	Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (k), (k), &((*(af))[(*(m))-(*(k))+0][(*(n))-(*(k))+0]), (lda), &((*(l))[(*(m))-(*(k))+0][(*(n))-(*(k))+0]), (lda))
	//
	//     Compute L(m-n+1:m,n-k+1:n) - q(1:m,m-n+1:m)' * a(1:m,n-k+1:n)
	//
	Dgemm(func() *[]byte {y :=[]byte("Transpose"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), (n), (k), (m), -(*one), (q), (lda), &((*(a))[0][(*(n))-(*(k))+0]), (lda), one, &((*(l))[(*(m))-(*(n))+0][(*(n))-(*(k))+0]), (lda))
	//
	//     Compute norm( L - Q'*a) / ( m * norm(a) * eps) .
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (k), &((*(a))[0][(*(n))-(*(k))+0]), (lda), (rwork)))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (k), &((*(l))[(*(m))-(*(n))+0][(*(n))-(*(k))+0]), (lda), (rwork)))
	if (*anorm) > (*zero) {
		(*(result))[0] = (((*resid) / DBLE(MAX(1, (*(m))))) / (*anorm)) / (*eps)
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute I - Q'*Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, (l), (lda))
	Dsyrk(func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (n), (m), -(*one), (q), (lda), one, (l), (lda))
	//
	//     Compute norm( I - Q'*Q) / ( m * eps) .
	//
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (n), (l), (lda), (rwork)))
	//
	(*(result))[1] = ((*resid) / DBLE(MAX(1, (*(m))))) / (*eps)
	//
	return
	//
	//     End of Dqlt02
	//
}
