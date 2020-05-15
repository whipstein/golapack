package goblas

import 

// Drzt02 returns
//      || I - Q'*Q || / ( m * eps)
// where the matrix Q is defined by the Householder transformations
// generated by Dtzrzf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Drzt02( m, n, af, lda, tau, work,
//                        lwork)
//
//       .. Scalar Arguments ..
//       intEGER            lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   af( lda, *), tau(*), work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Drzt02 returns
//      || I - Q'*Q || / ( m * eps)
// where the matrix Q is defined by the Householder transformations
// generated by Dtzrzf.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix af.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix af.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          The output of Dtzrzf.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array af.
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
//          length of work array. lwork >= N*n+N*nb.
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
func Drzt02(m *int, n *int, af *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int) (Drzt02Return *float64) {
	Drzt02Return = new(float64)
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
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
	(*(Drzt02Return)) = (*zero)
	//
	if (*(lwork)) < (*(n))*(*(n))+(*(n)) {
		Xerbla(func() *[]byte {y :=[]byte("Drzt02"); return &y }(), func() *int {y := 7; return &y }())
		return
	}
	//
	//     Quick return if possible
	//
	if (*(m)) <= 0 || (*(n)) <= 0 {
		return
	}
	//
	//     q := I
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (n), (n), zero, one, (work), (n))
	//
	//     q := P1 * ... * P(m) * Q
	//
	Dormrz(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), (n), (n), (m), (*(n))-(*(m)), (af), (lda), (tau), (work), (n), &((*(work))[(*(n))*(*(n))+0]), (*(lwork))-(*(n))*(*(n)), info)
	//
	//     q := P(m) * ... * P1 * Q
	//
	Dormrz(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (n), (n), (m), (*(n))-(*(m)), (af), (lda), (tau), (work), (n), &((*(work))[(*(n))*(*(n))+0]), (*(lwork))-(*(n))*(*(n)), info)
	//
	//     q := Q - I
	//
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		(*(work))[((*i)-1)*(*(n))+(*i)-1] = (*(work))[((*i)-1)*(*(n))+(*i)-1] - (*one)
		//Label10:
	}
	//
	(*(Drzt02Return)) = Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y }(), (n), (n), (work), (n), rwork) / (Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()) * DBLE(MAX((*(m)), (*(n)))))
	return
	//
	//     End of Drzt02
	//
}
