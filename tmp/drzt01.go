package goblas

import 

// Drzt01 returns
//      || A - R*Q || / ( m * eps * ||A||)
// for an upper trapezoidal A that was factored with Dtzrzf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Drzt01( m, n, a, af, lda, tau, work,
//                        lwork)
//
//       .. Scalar Arguments ..
//       intEGER            lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), af( lda, *), tau(*),
//      $                   work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Drzt01 returns
//      || A - R*Q || / ( m * eps * ||A||)
// for an upper trapezoidal A that was factored with Dtzrzf.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrices A and af.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrices A and af.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original upper trapezoidal M by N matrix A.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          The output of Dtzrzf for input matrix A.
//          The lower triangle is not referenced.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the arrays A and af.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (m)
//          Details of the Householder transformations as returned by
//          Dtzrzf.
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
//          The length of the array work.  lwork >= m*n + m*nb.
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
func Drzt01(m *int, n *int, a *[][]float64, af *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int) (Drzt01Return *float64) {
	Drzt01Return = new(float64)
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	norma := new(float64)
	rwork := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
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
	//     .. Executable Statements ..
	//
	(*(Drzt01Return)) = (*zero)
	//
	if (*(lwork)) < (*(m))*(*(n))+(*(m)) {
		Xerbla(func() *[]byte {y :=[]byte("Drzt01"); return &y }(), func() *int {y := 8; return &y }())
		return
	}
	//
	//     Quick return if possible
	//
	if (*(m)) <= 0 || (*(n)) <= 0 {
		return
	}
	//
	(*norma) = (*Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y }(), (m), (n), (a), (lda), rwork))
	//
	//     Copy upper triangle R
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, (work), (m))
	for (*j) = 1; (*j) <= (*(m)); (*j)++ {
		for (*i) = 1; (*i) <= (*j); (*i)++ {
			(*(work))[((*j)-1)*(*(m))+(*i)-1] = (*(af))[(*i)-1][(*j)-1]
			//Label10:
		}
		//Label20:
	}
	//
	//     R = r * P1 * ... *P(m)
	//
	Dormrz(func() *[]byte {y :=[]byte("Right"); return &y }(), func() *[]byte {y :=[]byte("No tranpose"); return &y }(), (m), (n), (m), (*(n))-(*(m)), (af), (lda), (tau), (work), (m), &((*(work))[(*(m))*(*(n))+0]), (*(lwork))-(*(m))*(*(n)), info)
	//
	//     R = R - A
	//
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		Daxpy((m), -(*one), &((*(a))[0][(*i)-1]), func() *int {y := 1; return &y }(), &((*(work))[((*i)-1)*(*(m))+0]), func() *int {y := 1; return &y }())
		//Label30:
	}
	//
	(*(Drzt01Return)) = (*Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y }(), (m), (n), (work), (m), rwork))
	//
	(*(Drzt01Return)) = (*(Drzt01Return)) / (Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()) * DBLE(MAX((*(m)), (*(n)))))
	if (*norma) != (*zero) {
		(*(Drzt01Return)) = (*(Drzt01Return)) / (*norma)
	}
	//
	return
	//
	//     End of Drzt01
	//
}
