package goblas

import 

// dlatm6 generates test matrices for the generalized eigenvalue
// problem, their corresponding right and left eigenvector matrices,
// and also reciprocal condition numbers for all eigenvalues and
// the reciprocal condition numbers of eigenvectors corresponding to
// the 1th and 5th eigenvalues.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE dlatm6( _type, n, a, lda, b, x, ldx, Y, ldy, alpha,
//                          beta, wx, wy, S, dif)
//
//       .. Scalar Arguments ..
//       inTEGER            lda, ldx, ldy, n, _type
//       DOUBLE PRECISION   alpha, beta, wx, wy
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( lda, *), dif(*), S(*),
//      $                   X( ldx, *), Y( ldy, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// dlatm6 generates test matrices for the generalized eigenvalue
// problem, their corresponding right and left eigenvector matrices,
// and also reciprocal condition numbers for all eigenvalues and
// the reciprocal condition numbers of eigenvectors corresponding to
// the 1th and 5th eigenvalues.
//
// Test Matrices
// =============
//
// Two kinds of test matrix pairs
//
//       (A, B) = inverse(YH) * (Da, Db) * inverse(x)
//
// are used in the tests:
//
// Type 1:
//    Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
//          0   2+a   0    0    0         0   1   0   0   0
//          0    0   3+a   0    0         0   0   1   0   0
//          0    0    0   4+a   0         0   0   0   1   0
//          0    0    0    0   5+a,      0   0   0   0   1, and
//
// Type 2:
//    Da =  1   -1    0    0    0    Db = 1   0   0   0   0
//          1    1    0    0    0         0   1   0   0   0
//          0    0    1    0    0         0   0   1   0   0
//          0    0    0   1+a  1+b        0   0   0   1   0
//          0    0    0  -1-b  1+a,      0   0   0   0   1 .
//
// In both cases the same inverse(YH) and inverse(x) are used to compute
// (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
//
// YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
//         0    1   -y    y   -y         0   1   x  -x  -x
//         0    0    1    0    0         0   0   1   0   0
//         0    0    0    1    0         0   0   0   1   0
//         0    0    0    0    1,        0   0   0   0   1 ,
//
// where a, b, x and y will have all values independently of each other.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] _type
// \verbatim
//          _type is inTEGER
//          Specifies the problem type (see further details).
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          Size of the matrices A and B.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda, N).
//          On exit A N-by-N is initialized according to _type.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of A and of B.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (lda, N).
//          On exit B N-by-N is initialized according to _type.
// \endverbatim
//
// \param[out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx, N).
//          On exit X is the N-by-N matrix of right eigenvectors.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is inTEGER
//          The leading dimension of X.
// \endverbatim
//
// \param[out] Y
// \verbatim
//          Y is DOUBLE PRECISION array, dimension (ldy, N).
//          On exit Y is the N-by-N matrix of left eigenvectors.
// \endverbatim
//
// \param[in] ldy
// \verbatim
//          ldy is inTEGER
//          The leading dimension of Y.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is DOUBLE PRECISION
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is DOUBLE PRECISION
//
//          Weighting _constants for matrix A.
// \endverbatim
//
// \param[in] wx
// \verbatim
//          wx is DOUBLE PRECISION
//          _constant for right eigenvector matrix.
// \endverbatim
//
// \param[in] wy
// \verbatim
//          wy is DOUBLE PRECISION
//          _constant for left eigenvector matrix.
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension (n)
//          S(i) is the reciprocal condition number for eigenvalue i.
// \endverbatim
//
// \param[out] dif
// \verbatim
//          dif is DOUBLE PRECISION array, dimension (n)
//          dif(i) is the reciprocal condition number for eigenvector i.
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
// \ingroup double_matgen
//
//  =====================================================================
func dlatm6(_type *int, n *int, a *[][]float64, lda *int, b *[][]float64, x *[][]float64, ldx *int, y *[][]float64, ldy *int, alpha *float64, beta *float64, wx *float64, wy *float64, s *[]float64, dif *[]float64) {
	zero := new(float64)
	one := new(float64)
	two := new(float64)
	three := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	work := func() *[]float64 {
		arr := make([]float64, 100)
		return &arr
	}()
	z := func() *[][]float64 {
		arr := make([][]float64, 12)
		for u := 0; u < 12; u++ {
			arr[u] = make([]float64, 12)
		}
		return &arr
	}()
	//
	//  -- lapACK computational routine (version 3.7.0) --
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
	(*two) = 2.0e+0
	(*three) = 3.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Generate test problem ...
	//     (Da, Db) ...
	//
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			//
			if (*i) == (*j) {
				(*(a))[(*i)-1][(*i)-1] = DBLE((*i)) + (*(alpha))
				(*(b))[(*i)-1][(*i)-1] = (*one)
			} else {
				(*(a))[(*i)-1][(*j)-1] = (*zero)
				(*(b))[(*i)-1][(*j)-1] = (*zero)
			}
			//
			//Label10:
		}
		//Label20:
	}
	//
	//     Form X and Y
	//
	Dlacpy(func() *byte {y := byte('F'); return &y}(), (n), (n), (b), (lda), (y), (ldy))
	(*(y))[2][0] = -(*(wy))
	(*(y))[3][0] = (*(wy))
	(*(y))[4][0] = -(*(wy))
	(*(y))[2][1] = -(*(wy))
	(*(y))[3][1] = (*(wy))
	(*(y))[4][1] = -(*(wy))
	//
	Dlacpy(func() *byte {y := byte('F'); return &y}(), (n), (n), (b), (lda), (x), (ldx))
	(*(x))[0][2] = -(*(wx))
	(*(x))[0][3] = -(*(wx))
	(*(x))[0][4] = (*(wx))
	(*(x))[1][2] = (*(wx))
	(*(x))[1][3] = -(*(wx))
	(*(x))[1][4] = -(*(wx))
	//
	//     Form (A, B)
	//
	(*(b))[0][2] = (*(wx)) + (*(wy))
	(*(b))[1][2] = -(*(wx)) + (*(wy))
	(*(b))[0][3] = (*(wx)) - (*(wy))
	(*(b))[1][3] = (*(wx)) - (*(wy))
	(*(b))[0][4] = -(*(wx)) + (*(wy))
	(*(b))[1][4] = (*(wx)) + (*(wy))
	if (*(_type)) == 1 {
		(*(a))[0][2] = (*(wx))*(*(a))[0][0] + (*(wy))*(*(a))[2][2]
		(*(a))[1][2] = -(*(wx))*(*(a))[1][1] + (*(wy))*(*(a))[2][2]
		(*(a))[0][3] = (*(wx))*(*(a))[0][0] - (*(wy))*(*(a))[3][3]
		(*(a))[1][3] = (*(wx))*(*(a))[1][1] - (*(wy))*(*(a))[3][3]
		(*(a))[0][4] = -(*(wx))*(*(a))[0][0] + (*(wy))*(*(a))[4][4]
		(*(a))[1][4] = (*(wx))*(*(a))[1][1] + (*(wy))*(*(a))[4][4]
	} else if (*(_type)) == 2 {
		(*(a))[0][2] = (*two)*(*(wx)) + (*(wy))
		(*(a))[1][2] = (*(wy))
		(*(a))[0][3] = -(*(wy)) * ((*two) + (*(alpha)) + (*(beta)))
		(*(a))[1][3] = (*two)*(*(wx)) - (*(wy))*((*two)+(*(alpha))+(*(beta)))
		(*(a))[0][4] = -(*two)*(*(wx)) + (*(wy))*((*(alpha))-(*(beta)))
		(*(a))[1][4] = (*(wy)) * ((*(alpha)) - (*(beta)))
		(*(a))[0][0] = (*one)
		(*(a))[0][1] = -(*one)
		(*(a))[1][0] = (*one)
		(*(a))[1][1] = (*(a))[0][0]
		(*(a))[2][2] = (*one)
		(*(a))[3][3] = (*one) + (*(alpha))
		(*(a))[3][4] = (*one) + (*(beta))
		(*(a))[4][3] = -(*(a))[3][4]
		(*(a))[4][4] = (*(a))[3][3]
	}
	//
	//     Compute condition numbers
	//
	if (*(_type)) == 1 {
		//
		(*(s))[0] = (*one) / SQRt(((*one)+(*three)*(*(wy))*(*(wy)))/((*one)+(*(a))[0][0]*(*(a))[0][0]))
		(*(s))[1] = (*one) / SQRt(((*one)+(*three)*(*(wy))*(*(wy)))/((*one)+(*(a))[1][1]*(*(a))[1][1]))
		(*(s))[2] = (*one) / SQRt(((*one)+(*two)*(*(wx))*(*(wx)))/((*one)+(*(a))[2][2]*(*(a))[2][2]))
		(*(s))[3] = (*one) / SQRt(((*one)+(*two)*(*(wx))*(*(wx)))/((*one)+(*(a))[3][3]*(*(a))[3][3]))
		(*(s))[4] = (*one) / SQRt(((*one)+(*two)*(*(wx))*(*(wx)))/((*one)+(*(a))[4][4]*(*(a))[4][4]))
		//
		Dlakf2(func() *int {y := 1; return &y}(), func() *int {y := 4; return &y}(), (a), (lda), &((*(a))[1][1]), (b), &((*(b))[1][1]), Z, func() *int {y := 12; return &y}())
		Dgesvd(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 8; return &y}(), func() *int {y := 8; return &y}(), Z, func() *int {y := 12; return &y}(), work, &((*work)[8]), func() *int {y := 1; return &y}(), &((*work)[9]), func() *int {y := 1; return &y}(), &((*work)[10]), func() *int {y := 40; return &y}(), info)
		(*(dif))[0] = (*work)[7]
		//
		Dlakf2(func() *int {y := 4; return &y}(), func() *int {y := 1; return &y}(), (a), (lda), &((*(a))[4][4]), (b), &((*(b))[4][4]), Z, func() *int {y := 12; return &y}())
		Dgesvd(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 8; return &y}(), func() *int {y := 8; return &y}(), Z, func() *int {y := 12; return &y}(), work, &((*work)[8]), func() *int {y := 1; return &y}(), &((*work)[9]), func() *int {y := 1; return &y}(), &((*work)[10]), func() *int {y := 40; return &y}(), info)
		(*(dif))[4] = (*work)[7]
		//
	} else if (*(_type)) == 2 {
		//
		(*(s))[0] = (*one) / SQRt((*one)/(*three)+(*(wy))*(*(wy)))
		(*(s))[1] = (*(s))[0]
		(*(s))[2] = (*one) / SQRt((*one)/(*two)+(*(wx))*(*(wx)))
		(*(s))[3] = (*one) / SQRt(((*one)+(*two)*(*(wx))*(*(wx)))/((*one)+((*one)+(*(alpha)))*((*one)+(*(alpha)))+((*one)+(*(beta)))*((*one)+(*(beta)))))
		(*(s))[4] = (*(s))[3]
		//
		Dlakf2(func() *int {y := 2; return &y}(), func() *int {y := 3; return &y}(), (a), (lda), &((*(a))[2][2]), (b), &((*(b))[2][2]), Z, func() *int {y := 12; return &y}())
		Dgesvd(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 12; return &y}(), func() *int {y := 12; return &y}(), Z, func() *int {y := 12; return &y}(), work, &((*work)[12]), func() *int {y := 1; return &y}(), &((*work)[13]), func() *int {y := 1; return &y}(), &((*work)[14]), func() *int {y := 60; return &y}(), info)
		(*(dif))[0] = (*work)[11]
		//
		Dlakf2(func() *int {y := 3; return &y}(), func() *int {y := 2; return &y}(), (a), (lda), &((*(a))[3][3]), (b), &((*(b))[3][3]), Z, func() *int {y := 12; return &y}())
		Dgesvd(func() *byte {y := byte('N'); return &y}(), func() *byte {y := byte('N'); return &y}(), func() *int {y := 12; return &y}(), func() *int {y := 12; return &y}(), Z, func() *int {y := 12; return &y}(), work, &((*work)[12]), func() *int {y := 1; return &y}(), &((*work)[13]), func() *int {y := 1; return &y}(), &((*work)[14]), func() *int {y := 60; return &y}(), info)
		(*(dif))[4] = (*work)[11]
		//
	}
	//
	return
	//
	//     End of dlatm6
	//
}
