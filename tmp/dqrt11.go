package goblas

import 

// Dqrt11 computes the test ratio
//
//       || Q'*Q - I || / (eps * m)
//
// where the orthogonal matrix Q is represented as a product of
// elementary transformations.  Each transformation has the form
//
//    H(k) = I - tau(k) v(k) v(k)'
//
// where tau(k) is stored in tau(k) and v(k) is an m-vector of the form
//[0 ... 0 1 x(k)]', where x(k) is a vector of length m-k stored
// in a(k+1:m,k).
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dqrt11( m, k, a, lda, tau, work, lwork)
//
//       .. Scalar Arguments ..
//       inTEGER            k, lda, lwork, M
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), tau(*), work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt11 computes the test ratio
//
//       || Q'*Q - I || / (eps * m)
//
// where the orthogonal matrix Q is represented as a product of
// elementary transformations.  Each transformation has the form
//
//    H(k) = I - tau(k) v(k) v(k)'
//
// where tau(k) is stored in tau(k) and v(k) is an m-vector of the form
//[0 ... 0 1 x(k)]', where x(k) is a vector of length m-k stored
// in a(k+1:m,k).
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows of the matrix A.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is inTEGER
//          The number of columns of A whose subdiagonal entries
//          contain information about orthogonal transformations.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,k)
//          The (possibly partial) output of a QR reduction routine.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (k)
//          The scaling factors tau for the elementary transformations as
//          computed by the QR factorization routine.
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
//          The length of the array work.  lwork >= M*M + M.
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
func Dqrt11(m *int, k *int, a *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int) (dqrt11Return *float64) {
	dqrt11Return = new(float64)
	zero := new(float64)
	one := new(float64)
	info := new(int)
	j := new(int)
	rdummy := func() *[]float64 {
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
	(*zero) = 0.0
	(*one) = 1.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Executable Statements ..
	//
	(*(dqrt11Return)) = (*zero)
	//
	//     Test for sufficient workspace
	//
	if (*(lwork)) < (*(m))*(*(m))+(*(m)) {
		Xerbla(func() *[]byte {y :=[]byte("Dqrt11"); return &y}(), func() *int {y := 7; return &y}())
		return
	}
	//
	//     Quick return if possible
	//
	if (*(m)) <= 0 {
		return
	}
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (m), zero, one, (work), (m))
	//
	//     Form Q
	//
	Dorm2r(func() *[]byte {y :=[]byte("Left"); return &y}(), func() *[]byte {y :=[]byte("No transpose"); return &y}(), (m), (m), (k), (a), (lda), (tau), (work), (m), &((*(work))[(*(m))*(*(m))+0]), info)
	//
	//     Form Q'*Q
	//
	Dorm2r(func() *[]byte {y :=[]byte("Left"); return &y}(), func() *[]byte {y :=[]byte("Transpose"); return &y}(), (m), (m), (k), (a), (lda), (tau), (work), (m), &((*(work))[(*(m))*(*(m))+0]), info)
	//
	for (*j) = 1; (*j) <= (*(m)); (*j)++ {
		(*(work))[((*j)-1)*(*(m))+(*j)-1] = (*(work))[((*j)-1)*(*(m))+(*j)-1] - (*one)
		//Label10:
	}
	//
	(*(dqrt11Return)) = Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y}(), (m), (m), (work), (m), rdummy) / (DBLE((*(m))) * Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	//
	return
	//
	//     End of Dqrt11
	//
}
