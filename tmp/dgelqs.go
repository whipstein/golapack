package goblas

import 

// Dgelqs computes a minimum-norm solution
//     min || A*X - B ||
// using the LQ factorization
//     A = L*Q
// computed by Dgelqf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dgelqs( m, n, nrhs, a, lda, tau, b, ldb, work, lwork,
//                          info)
//
//       .. Scalar Arguments ..
//       intEGER            info, lda, ldb, lwork, m, n, nrhs
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
// using the LQ factorization
//     A = L*Q
// computed by Dgelqf.
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
//          The number of columns of the matrix A.  N >= M >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of B.  nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the LQ factorization of the original matrix A as
//          returned by Dgelqf.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
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
//          On entry, the m-by-nrhs right hand side matrix B.
//          On exit, the n-by-nrhs solution matrix X.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B. ldb >= N.
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
//          The length of the array work.  lwork must be at least nrhs,
//          and should be at least nrhs*nb, where nb is the block size
//          for this environment.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
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
func Dgelqs(m *int, n *int, nrhs *int, a *[][]float64, lda *int, tau *[]float64, b *[][]float64, ldb *int, work *[]float64, lwork *int, info *int) {
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
		Xerbla(func() *[]byte {y :=[]byte("Dgelqs"); return &y }(), -(*(info)))
		return
	}
	//
	//     Quick return if possible
	//
	if (*(n)) == 0 || (*(nrhs)) == 0 || (*(m)) == 0 {
		return
	}
	//
	//     Solve L*X = B(1:m,:)
	//
	Dtrsm(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Lower"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Non-unit"); return &y }(), (m), (nrhs), one, (a), (lda), (b), (ldb))
	//
	//     Set B(m+1:n,:) to zero
	//
	if (*(m)) < (*(n)) {
		Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (*(n))-(*(m)), (nrhs), zero, zero, &((*(b))[(*(m))+0][0]), (ldb))
	}
	//
	//     b := Q' * B
	//
	Dormlq(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (n), (nrhs), (m), (a), (lda), (tau), (b), (ldb), (work), (lwork), (info))
	//
	return
	//
	//     End of Dgelqs
	//
}
