package goblas

import 

// Dgeqrs solves the least squares problem
//     min || A*X - B ||
// using the QR factorization
//     A = Q*R
// computed by Dgeqrf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dgeqrs( m, n, nrhs, a, lda, tau, b, ldb, work, lwork,
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
// Solve the least squares problem
//     min || A*X - B ||
// using the QR factorization
//     A = Q*R
// computed by Dgeqrf.
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
//          The number of columns of the matrix A.  M >= N >= 0.
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
//          Details of the QR factorization of the original matrix A as
//          returned by Dgeqrf.
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
//          tau is DOUBLE PRECISION array, dimension (n)
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
//          The leading dimension of the array B. ldb >= M.
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
func Dgeqrs(m *int, n *int, nrhs *int, a *[][]float64, lda *int, tau *[]float64, b *[][]float64, ldb *int, work *[]float64, lwork *int, info *int) {
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
	(*one) = 1.0e+0
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input arguments.
	//
	(*(info)) = 0
	if (*(m)) < 0 {
		(*(info)) = -1
	} else if (*(n)) < 0 || (*(n)) > (*(m)) {
		(*(info)) = -2
	} else if (*(nrhs)) < 0 {
		(*(info)) = -3
	} else if (*(lda)) < (MAX(1, (*(m)))) {
		(*(info)) = -5
	} else if (*(ldb)) < (MAX(1, (*(m)))) {
		(*(info)) = -8
	} else if (*(lwork)) < 1 || (*(lwork)) < (*(nrhs)) && (*(m)) > 0 && (*(n)) > 0 {
		(*(info)) = -10
	}
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dgeqrs"); return &y }(), -(*(info)))
		return
	}
	//
	//     Quick return if possible
	//
	if (*(n)) == 0 || (*(nrhs)) == 0 || (*(m)) == 0 {
		return
	}
	//
	//     b := Q' * B
	//
	Dormqr(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Transpose"); return &y }(), (m), (nrhs), (n), (a), (lda), (tau), (b), (ldb), (work), (lwork), (info))
	//
	//     Solve R*X = B(1:n,:)
	//
	Dtrsm(func() *[]byte {y :=[]byte("Left"); return &y }(), func() *[]byte {y :=[]byte("Upper"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Non-unit"); return &y }(), (n), (nrhs), one, (a), (lda), (b), (ldb))
	//
	return
	//
	//     End of Dgeqrs
	//
}
