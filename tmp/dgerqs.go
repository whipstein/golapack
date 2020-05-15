package goblas

import 

// Dgerqs computes a minimum-norm solution
//     min || A*X - B ||
// using the RQ factorization
//     A = R*Q
// computed by Dgerqf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dgerqs( m, n, nrhs, a, lda, tau, b, ldb, work, lwork,
//                          info)
//
//       .. Scalar Arguments ..
//       inTEGER            info, lda, ldb, lwork, m, n, nrhs
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), tau(*),
//      $                   work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Compute a minimum-norm solution
//     min || A*X - B ||
// using the RQ factorization
//     A = R*Q
// computed by Dgerqf.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows of the matrix A.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns of the matrix A.  N >= M >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The number of columns of B.  nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the RQ factorization of the original matrix A as
//          returned by Dgerqf.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.  lda >= M.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (m)
//          Details of the orthogonal matrix Q.
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          On entry, the right hand side vectors for the linear system.
//          On exit, the solution vectors X.  Each solution vector
//          is contained in rows 1:N of a column of B.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B. ldb >= max(1,N).
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
//          The length of the array work.  lwork must be at least nrhs,
//          and should be at least nrhs*nb, where nb is the block size
//          for this environment.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is inTEGER
//          = 0: successful exit
//          < 0: if info = -i, the i-th argument had an illegal value
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
func Dgerqs(m *int, n *int, nrhs *int, a *[][]float64, lda *int, tau *[]float64, b *[][]float64, ldb *int, work *[]float64, lwork *int, info *int) {
	zero := new(float64)
	one := new(float64)
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
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input parameters.
	//
	(*(info)) = 0
	if (*(m)) < 0 {
		(*(info)) = -1
	} else if (*(n)) < 0 || (*(m)) > (*(n)) {
		(*(info)) = -2
	} else if (*(nrhs)) < 0 {
		(*(info)) = -3
	} else if (*(lda)) < (MAX(1, (*(m)))) {
		(*(info)) = -5
	} else if (*(ldb)) < (MAX(1, (*(n)))) {
		(*(info)) = -8
	} else if (*(lwork)) < 1 || (*(lwork)) < (*(nrhs)) && (*(m)) > 0 && (*(n)) > 0 {
		(*(info)) = -10
	}
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), -(*(info)))
		return
	}
	//
	//     Quick return if possible
	//
	if (*(n)) == 0 || (*(nrhs)) == 0 || (*(m)) == 0 {
		return
	}
	//
	//     Solve R*X = B(n-m+1:n,:)
	//
	Dtrsm(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Non-unit"); return &y }(), (m), (nrhs), one, &((*(a))[0][(*(n))-(*(m))+0]), (lda), &((*(b))[(*(n))-(*(m))+0][0]), (ldb))
	//
	//     Set B(1:n-m,:) to zero
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (*(n))-(*(m)), (nrhs), zero, zero, (b), (ldb))
	//
	//     b := Q' * B
	//
	Dormrq(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (n), (nrhs), (m), (a), (lda), (tau), (b), (ldb), (work), (lwork), (info))
	//
	return
	//
	//     End of Dgerqs
	//
}
