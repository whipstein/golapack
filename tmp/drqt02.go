package goblas

import 

// drqt02 tests Dorgrq, which generates an m-by-n matrix Q with
// orthonornmal rows that is defined as the product of k elementary
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
//       SUBROUTinE drqt02( m, n, k, a, af, q, r, lda, tau, work, lwork,
//                          rwork, result)
//
//       .. Scalar Arguments ..
//       inTEGER            k, lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), af( lda, *), q( lda, *),
//      $                   r( lda, *), result(*), rwork(*), tau(*),
//      $                   work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// drqt02 tests Dorgrq, which generates an m-by-n matrix Q with
// orthonornmal rows that is defined as the product of k elementary
// reflectors.
//
// Given the RQ factorization of an m-by-n matrix a, drqt02 generates
// the orthogonal matrix Q defined by the factorization of the last k
// rows of A; it compares r(m-k+1:m,n-m+1:n) with
// a(m-k+1:m,1:n)*q(n-m+1:n,1:n)', and checks that the rows of Q are
// orthonormal.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows of the matrix Q to be generated.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns of the matrix Q to be generated.
//          N >= M >= 0.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is inTEGER
//          The number of elementary reflectors whose product defines the
//          matrix Q. M >= K >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The m-by-n matrix A which was factorized by drqt01.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the RQ factorization of a, as returned by Dgerqf.
//          See Dgerqf for further details.
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[out] R
// \verbatim
//          R is DOUBLE PRECISION array, dimension (lda,M)
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the arrays a, af, Q and L. lda >= N.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (m)
//          The scalar factors of the elementary reflectors corresponding
//          to the RQ factorization in af.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (lwork)
// \endverbatim
//
// \param[in] lwork
// \verbatim
//          lwork is inTEGER
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
//          result1 = norm( R - A*Q') / ( N * norm(a) * eps)
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
func drqt02(m *int, n *int, k *int, a *[][]float64, af *[][]float64, q *[][]float64, r *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
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
	//     .. Scalars in Common ..
	//     ..
	//     .. Common blocks ..
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
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	//
	//     Copy the last k rows of the factorization to the array Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (n), rogue, rogue, (q), (lda))
	if (*(k)) < (*(n)) {
		Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), (k), (*(n))-(*(k)), &((*(af))[(*(m))-(*(k))+0][0]), (lda), &((*(q))[(*(m))-(*(k))+0][0]), (lda))
	}
	if (*(k)) > 1 {
		Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y}(), (*(k))-1, (*(k))-1, &((*(af))[(*(m))-(*(k))+1][(*(n))-(*(k))+0]), (lda), &((*(q))[(*(m))-(*(k))+1][(*(n))-(*(k))+0]), (lda))
	}
	//
	//     Generate the last n rows of the matrix Q
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgrq"); return &y}()
	Dorgrq((m), (n), (k), (q), (lda), &((*(tau))[(*(m))-(*(k))+0]), (work), (lwork), info)
	//
	//     Copy r(m-k+1:m,n-m+1:n)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (k), (m), zero, zero, &((*(r))[(*(m))-(*(k))+0][(*(n))-(*(m))+0]), (lda))
	Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y}(), (k), (k), &((*(af))[(*(m))-(*(k))+0][(*(n))-(*(k))+0]), (lda), &((*(r))[(*(m))-(*(k))+0][(*(n))-(*(k))+0]), (lda))
	//
	//     Compute r(m-k+1:m,n-m+1:n) - a(m-k+1:m,1:n) * q(n-m+1:n,1:n)'
	//
	Dgemm(func() *[]byte {y :=[]byte("No transpose"); return &y}(), func() *[]byte {y :=[]byte("Transpose"); return &y}(), (k), (m), (n), -(*one), &((*(a))[(*(m))-(*(k))+0][0]), (lda), (q), (lda), one, &((*(r))[(*(m))-(*(k))+0][(*(n))-(*(m))+0]), (lda))
	//
	//     Compute norm( R - A*Q') / ( N * norm(a) * eps) .
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (k), (n), &((*(a))[(*(m))-(*(k))+0][0]), (lda), (rwork)))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y}(), (k), (m), &((*(r))[(*(m))-(*(k))+0][(*(n))-(*(m))+0]), (lda), (rwork)))
	if (*anorm) > (*zero) {
		(*(result))[0] = (((*resid) / DBLE(MAX(1, (*(n))))) / (*anorm)) / (*eps)
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute I - Q*Q'
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (m), zero, one, (r), (lda))
	Dsyrk(func() *[]byte {y :=[]byte("Upper"); return &y}(), func() *[]byte {y :=[]byte("No transpose"); return &y}(), (m), (n), -(*one), (q), (lda), one, (r), (lda))
	//
	//     Compute norm( I - Q*Q') / ( N * eps) .
	//
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y}(), func() *[]byte {y :=[]byte("Upper"); return &y}(), (m), (r), (lda), (rwork)))
	//
	(*(result))[1] = ((*resid) / DBLE(MAX(1, (*(n))))) / (*eps)
	//
	return
	//
	//     End of drqt02
	//
}
