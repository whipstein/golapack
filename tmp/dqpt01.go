package goblas

import 

// Dqpt01 tests the QR-factorization with pivoting of a matrix A.  The
// array af contains the (possibly partial) QR-factorization of a, where
// the upper triangle of af(1:k,1:k) is a partial triangular factor,
// the entries below the diagonal in the first k columns are the
// Householder vectors, and the rest of af contains a partially updated
// matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dqpt01( m, n, k, a, af, lda, tau, jpvt,
//                        work, lwork)
//
//       .. Scalar Arguments ..
//       inTEGER            k, lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       inTEGER            JPVt(*)
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
// Dqpt01 tests the QR-factorization with pivoting of a matrix A.  The
// array af contains the (possibly partial) QR-factorization of a, where
// the upper triangle of af(1:k,1:k) is a partial triangular factor,
// the entries below the diagonal in the first k columns are the
// Householder vectors, and the rest of af contains a partially updated
// matrix.
//
// This function returns ||A*P - Q*R||/(||norm(a)||*eps*M)
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows of the matrices A and af.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns of the matrices A and af.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is inTEGER
//          The number of columns of af that have been reduced
//          to upper triangular form.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda, N)
//          The original matrix A.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          The (possibly partial) output of DGEQPF.  The upper triangle
//          of af(1:k,1:k) is a partial triangular factor, the entries
//          below the diagonal in the first k columns are the Householder
//          vectors, and the rest of af contains a partially updated
//          matrix.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the arrays A and af.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (k)
//          Details of the Householder transformations as returned by
//          DGEQPF.
// \endverbatim
//
// \param[in] jpvt
// \verbatim
//          jpvt is inTEGER array, dimension (n)
//          Pivot information as returned by DGEQPF.
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
//          The length of the array work.  lwork >= M*n+N.
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
func Dqpt01(m *int, n *int, k *int, a *[][]float64, af *[][]float64, lda *int, tau *[]float64, jpvt *[]int, work *[]float64, lwork *int) (dqpt01Return *float64) {
	dqpt01Return = new(float64)
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
	(*zero) = 0.0
	(*one) = 1.0
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
	(*(dqpt01Return)) = (*zero)
	//
	//     Test if there is enough workspace
	//
	if (*(lwork)) < (*(m))*(*(n))+(*(n)) {
		Xerbla(func() *[]byte {y :=[]byte("Dqpt01"); return &y}(), func() *int {y := 10; return &y}())
		return
	}
	//
	//     Quick return if possible
	//
	if (*(m)) <= 0 || (*(n)) <= 0 {
		return
	}
	//
	(*norma) = (*Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y}(), (m), (n), (a), (lda), rwork))
	//
	for (*j) = 1; (*j) <= (*(k)); (*j)++ {
		for (*i) = 1; (*i) <= (Min((*j), (*(m)))); (*i)++ {
			(*(work))[((*j)-1)*(*(m))+(*i)-1] = (*(af))[(*i)-1][(*j)-1]
			//Label10:
		}
		for (*i) = (*j) + 1; (*i) <= (*(m)); (*i)++ {
			(*(work))[((*j)-1)*(*(m))+(*i)-1] = (*zero)
			//Label20:
		}
		//Label30:
	}
	for (*j) = (*(k)) + 1; (*j) <= (*(n)); (*j)++ {
		Dcopy((m), &((*(af))[0][(*j)-1]), func() *int {y := 1; return &y}(), &((*(work))[((*j)-1)*(*(m))+0]), func() *int {y := 1; return &y}())
		//Label40:
	}
	//
	Dormqr(func() *[]byte {y :=[]byte("Left"); return &y}(), func() *[]byte {y :=[]byte("No transpose"); return &y}(), (m), (n), (k), (af), (lda), (tau), (work), (m), &((*(work))[(*(m))*(*(n))+0]), (*(lwork))-(*(m))*(*(n)), info)
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		//
		//        Compare i-th column of QR and jpvt(i)-th column of A
		//
		Daxpy((m), -(*one), &((*(a))[0][(*(jpvt))[(*j)-1]-1]), func() *int {y := 1; return &y}(), &((*(work))[((*j)-1)*(*(m))+0]), func() *int {y := 1; return &y}())
		//Label50:
	}
	//
	(*(dqpt01Return)) = Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y}(), (m), (n), (work), (m), rwork) / (DBLE(MAX((*(m)), (*(n)))) * Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	if (*norma) != (*zero) {
		(*(dqpt01Return)) = (*(dqpt01Return)) / (*norma)
	}
	//
	return
	//
	//     End of Dqpt01
	//
}
