package goblas

import 

// Dqrt15 generates a matrix with full or deficient rank and of various
// norms.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqrt15( scale, rksel, m, n, nrhs, a, lda, b, ldb, S,
//                          rank, norma, normb, iseed, work, lwork)
//
//       .. Scalar Arguments ..
//       inTEGER            lda, ldb, lwork, m, n, nrhs, rank, rksel, scale
//       DOUBLE PRECISION   norma, normb
//       ..
//       .. Array Arguments ..
//       inTEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), S(*), work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt15 generates a matrix with full or deficient rank and of various
// norms.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] scale
// \verbatim
//          scale is inTEGER
//          scale = 1: normally scaled matrix
//          scale = 2: matrix scaled up
//          scale = 3: matrix scaled down
// \endverbatim
//
// \param[in] rksel
// \verbatim
//          rksel is inTEGER
//          rksel = 1: full rank matrix
//          rksel = 2: rank-deficient matrix
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number of rows of the matrix A.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns of A.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The number of columns of B.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The M-by-N matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb, nrhs)
//          A matrix that is in the range space of matrix A.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B.
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension Min(m,N)
//          Singular values of A.
// \endverbatim
//
// \param[out] rank
// \verbatim
//          rank is inTEGER
//          number of nonzero singular values of A.
// \endverbatim
//
// \param[out] norma
// \verbatim
//          norma is DOUBLE PRECISION
//          one-norm of A.
// \endverbatim
//
// \param[out] normb
// \verbatim
//          normb is DOUBLE PRECISION
//          one-norm of B.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is integer array, dimension (4)
//          seed for random number generator.
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
//          length of work space required.
//          lwork >= MAX(M+Min(m,N),nrhs*Min(m,N),2*n+M)
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
func Dqrt15(scale *int, rksel *int, m *int, n *int, nrhs *int, a *[][]float64, lda *int, b *[][]float64, ldb *int, s *[]float64, rank *int, norma *float64, normb *float64, iseed *[]int, work *[]float64, lwork *int) {
	zero := new(float64)
	one := new(float64)
	two := new(float64)
	svmin := new(float64)
	info := new(int)
	j := new(int)
	mn := new(int)
	bignum := new(float64)
	eps := new(float64)
	smlnum := new(float64)
	temp := new(float64)
	dummy := func() *[]float64 {
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
	(*two) = 2.0
	(*svmin) = 0.1
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
	(*mn) = (Min((*(m)), (*(n))))
	if (*(lwork)) < (*MAX((*(m))+(*mn), (*mn)*(*(nrhs)), 2*(*(n))+(*(m)))) {
		Xerbla(func() *[]byte {y :=[]byte("Dqrt15"); return &y}(), func() *int {y := 16; return &y}())
		return
	}
	//
	(*smlnum) = (*Dlamch(func() *[]byte {y :=[]byte("Safe minimum"); return &y}()))
	(*bignum) = (*one) / (*smlnum)
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	(*smlnum) = ((*smlnum) / (*eps)) / (*eps)
	(*bignum) = (*one) / (*smlnum)
	//
	//     Determine rank and (unscaled) singular values
	//
	if (*(rksel)) == 1 {
		(*(rank)) = (*mn)
	} else if (*(rksel)) == 2 {
		(*(rank)) = (3 * (*mn)) / 4
		for (*j) = (*(rank)) + 1; (*j) <= (*mn); (*j)++ {
			(*(s))[(*j)-1] = (*zero)
			//Label10:
		}
	} else {
		Xerbla(func() *[]byte {y :=[]byte("Dqrt15"); return &y}(), func() *int {y := 2; return &y}())
	}
	//
	if (*(rank)) > 0 {
		//
		//        Nontrivial case
		//
		(*(s))[0] = (*one)
		for (*j) = 2; (*j) <= (*(rank)); (*j)++ {
		Label20:
			;
			(*temp) = (*Dlarnd(func() *int {y := 1; return &y}(), (iseed)))
			if (*temp) > (*svmin) {
				(*(s))[(*j)-1] = (ABS((*temp)))
			} else {
				goto Label20
			}
			//Label30:
		}
		Dlaord(func() *[]byte {y :=[]byte("Decreasing"); return &y}(), (rank), (s), func() *int {y := 1; return &y}())
		//
		//        Generate 'rank' columns of a random orthogonal matrix in A
		//
		Dlarnv(func() *int {y := 2; return &y}(), (iseed), (m), (work))
		Dscal((m), (*one)/Dnrm2((m), (work), func() *int {y := 1; return &y}()), (work), func() *int {y := 1; return &y}())
		Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (rank), zero, one, (a), (lda))
		dlarf(func() *[]byte {y :=[]byte("Left"); return &y}(), (m), (rank), (work), func() *int {y := 1; return &y}(), two, (a), (lda), &((*(work))[(*(m))+0]))
		//
		//        workspace used: m+mn
		//
		//        Generate consistent rhs in the range space of A
		//
		Dlarnv(func() *int {y := 2; return &y}(), (iseed), (*(rank))*(*(nrhs)), (work))
		Dgemm(func() *[]byte {y :=[]byte("No transpose"); return &y}(), func() *[]byte {y :=[]byte("No transpose"); return &y}(), (m), (nrhs), (rank), one, (a), (lda), (work), (rank), zero, (b), (ldb))
		//
		//        work space used: <= mn *nrhs
		//
		//        generate (unscaled) matrix A
		//
		for (*j) = 1; (*j) <= (*(rank)); (*j)++ {
			Dscal((m), &((*(s))[(*j)-1]), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y}())
			//Label40:
		}
		if (*(rank)) < (*(n)) {
			Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (*(n))-(*(rank)), zero, zero, &((*(a))[0][(*(rank))+0]), (lda))
		}
		Dlaror(func() *[]byte {y :=[]byte("Right"); return &y}(), func() *[]byte {y :=[]byte("No initialization"); return &y}(), (m), (n), (a), (lda), (iseed), (work), info)
		//
	} else {
		//
		//        work space used 2*n+m
		//
		//        Generate null matrix and rhs
		//
		for (*j) = 1; (*j) <= (*mn); (*j)++ {
			(*(s))[(*j)-1] = (*zero)
			//Label50:
		}
		Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (n), zero, zero, (a), (lda))
		Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (m), (nrhs), zero, zero, (b), (ldb))
		//
	}
	//
	//     Scale the matrix
	//
	if (*(scale)) != 1 {
		(*(norma)) = (*Dlange(func() *[]byte {y :=[]byte("Max"); return &y}(), (m), (n), (a), (lda), dummy))
		if (*(norma)) != (*zero) {
			if (*(scale)) == 2 {
				//
				//              matrix scaled up
				//
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), bignum, (m), (n), (a), (lda), info)
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), bignum, mn, func() *int {y := 1; return &y}(), (s), mn, info)
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), bignum, (m), (nrhs), (b), (ldb), info)
			} else if (*(scale)) == 3 {
				//
				//              matrix scaled down
				//
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), smlnum, (m), (n), (a), (lda), info)
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), smlnum, mn, func() *int {y := 1; return &y}(), (s), mn, info)
				dlascl(func() *[]byte {y :=[]byte("General"); return &y}(), func() *int {y := 0; return &y}(), func() *int {y := 0; return &y}(), (norma), smlnum, (m), (nrhs), (b), (ldb), info)
			} else {
				Xerbla(func() *[]byte {y :=[]byte("Dqrt15"); return &y}(), func() *int {y := 1; return &y}())
				return
			}
		}
	}
	//
	(*(norma)) = (*Dasum(mn, (s), func() *int {y := 1; return &y}()))
	(*(normb)) = (*Dlange(func() *[]byte {y :=[]byte("one-norm"); return &y}(), (m), (nrhs), (b), (ldb), dummy))
	//
	return
	//
	//     End of Dqrt15
	//
}
