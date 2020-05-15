package goblas

import 

// Dqrt02 tests Dorgqr, which generates an m-by-n matrix Q with
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
//       SUBROUTinE Dqrt02( m, n, k, a, af, q, r, lda, tau, work, lwork,
//                          rwork, result)
//
//       .. Scalar Arguments ..
//       intEGER            k, lda, lwork, m, N
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
// Dqrt02 tests Dorgqr, which generates an m-by-n matrix Q with
// orthonornmal columns that is defined as the product of k elementary
// reflectors.
//
// Given the QR factorization of an m-by-n matrix a, Dqrt02 generates
// the orthogonal matrix Q defined by the factorization of the first k
// columns of A; it compares r(1:n,1:k) with q(1:m,1:n)'*a(1:m,1:k),
// and checks that the columns of Q are orthonormal.
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
//          The m-by-n matrix A which was factorized by Dqrt01.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the QR factorization of a, as returned by Dgeqrf.
//          See Dgeqrf for further details.
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[out] R
// \verbatim
//          R is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the arrays a, af, Q and R. lda >= M.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (n)
//          The scalar factors of the elementary reflectors corresponding
//          to the QR factorization in af.
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
//          result(1) = norm( R - Q'*a) / ( m * norm(a) * eps)
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
func Dqrt02(m *int, n *int, k *int, a *[][]float64, af *[][]float64, q *[][]float64, r *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
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
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     Copy the first k columns of the factorization to the array Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), rogue, rogue, (q), (lda))
	Dlacpy(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*(m))-1, (k), &((*(af))[1][0]), (lda), &((*(q))[1][0]), (lda))
	//
	//     Generate the first n columns of the matrix Q
	//
	(*srnamt) = *func() *[]byte {y :=[]byte("Dorgqr"); return &y }()
	Dorgqr((m), (n), (k), (q), (lda), (tau), (work), (lwork), info)
	//
	//     Copy r(1:n,1:k)
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (k), zero, zero, (r), (lda))
	Dlacpy(func() *[]byte {y :=[]byte("Upper"); return &y }(), (n), (k), (af), (lda), (r), (lda))
	//
	//     Compute r(1:n,1:k) - q(1:m,1:n)' * a(1:m,1:k)
	//
	Dgemm(func() *[]byte {y :=[]byte("Transpose"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), (n), (k), (m), -(*one), (q), (lda), (a), (lda), one, (r), (lda))
	//
	//     Compute norm( R - Q'*a) / ( m * norm(a) * eps) .
	//
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (k), (a), (lda), (rwork)))
	(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (n), (k), (r), (lda), (rwork)))
	if (*anorm) > (*zero) {
		(*(result))[0] = (((*resid) / DBLE(MAX(1, (*(m))))) / (*anorm)) / (*eps)
	} else {
		(*(result))[0] = (*zero)
	}
	//
	//     Compute I - Q'*Q
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, (r), (lda))
	Dsyrk(func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (n), (m), -(*one), (q), (lda), one, (r), (lda))
	//
	//     Compute norm( I - Q'*Q) / ( m * eps) .
	//
	(*resid) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), (n), (r), (lda), (rwork)))
	//
	(*(result))[1] = ((*resid) / DBLE(MAX(1, (*(m))))) / (*eps)
	//
	return
	//
	//     End of Dqrt02
	//
}
