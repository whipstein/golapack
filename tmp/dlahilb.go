package goblas

import 

// Dlahilb generates an N by N scaled Hilbert matrix in A along with
// nrhs right-hand sides in B and solutions in X such that A*X=B.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlahilb( n, nrhs, a, lda, x, ldx, b, ldb, work, info)
//
//       .. Scalar Arguments ..
//       inTEGER n, nrhs, lda, ldx, ldb, info
//       .. Array Arguments ..
//       DOUBLE PRECISION a(lda, N), X(ldx, nrhs), B(ldb, nrhs), work(n)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlahilb generates an N by N scaled Hilbert matrix in A along with
// nrhs right-hand sides in B and solutions in X such that A*X=B.
//
// The Hilbert matrix is scaled by M = LCM(1, 2, ..., 2*n-1) so that all
// entries are integers.  The right-hand sides are the first nrhs
// columns of m * the identity matrix, and the solutions are the
// first nrhs columns of the inverse Hilbert matrix.
//
// The condition number of the Hilbert matrix grows exponentially with
// its size, roughly as O(e ** (3.5*n)).  Additionally, the inverse
// Hilbert matrices beyond a relatively small dimension cannot be
// generated exactly without extra precision.  Precision is exhausted
// when the largest entry in the inverse Hilbert matrix is greater than
// 2 to the power of the number of bits in the fraction of the data type
// used plus one, which is 24 for single precision.
//
// In single, the generated solution is exact for N <= 6 and has
// small componentwise error for 7 <= N <= 11.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The dimension of the matrix A.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The requested number of right-hand sides.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda, N)
//          The generated scaled Hilbert matrix.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.  lda >= N.
// \endverbatim
//
// \param[out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx, nrhs)
//          The generated exact solutions.  Currently, the first nrhs
//          columns of the inverse Hilbert matrix.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is inTEGER
//          The leading dimension of the array X.  ldx >= N.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb, nrhs)
//          The generated right-hand sides.  Currently, the first nrhs
//          columns of LCM(1, 2, ..., 2*n-1) * the identity matrix.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B.  ldb >= N.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is inTEGER
//          = 0: successful exit
//          = 1: N is too large; the data is still generated but may not
//               be not exact.
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
// \date June 2017
//
// \ingroup double_lin
//
//  =====================================================================
func Dlahilb(n *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, b *[][]float64, ldb *int, work *[]float64, info *int) {
	tm := new(int)
	ti := new(int)
	r := new(int)
	m := new(int)
	i := new(int)
	j := new(int)
	nmaxExact := new(int)
	nmaxApprox := new(int)
	//
	//  -- lapACK test routine (version 3.7.1) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2017
	//
	//     .. Scalar Arguments ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//     .. Local Scalars ..
	//     ..
	//     .. Parameters ..
	//     nmaxExact   the largest dimension where the generated data is
	//                  exact.
	//     nmaxApprox  the largest dimension where the generated data has
	//                  a small componentwise relative error.
	(*nmaxExact) = 6
	(*nmaxApprox) = 11
	//     ..
	//     .. External Functions
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input arguments
	//
	(*(info)) = 0
	if (*(n)) < 0 || (*(n)) > (*nmaxApprox) {
		(*(info)) = -1
	} else if (*(nrhs)) < 0 {
		(*(info)) = -2
	} else if (*(lda)) < (*(n)) {
		(*(info)) = -4
	} else if (*(ldx)) < (*(n)) {
		(*(info)) = -6
	} else if (*(ldb)) < (*(n)) {
		(*(info)) = -8
	}
	if (*(info)) < 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dlahilb"); return &y}(), -(*(info)))
		return
	}
	if (*(n)) > (*nmaxExact) {
		(*(info)) = 1
	}
	//
	//     Compute M = the LCM of the integers[1, 2*n-1].  The largest
	//     reasonable N is small enough that integers suffice (up to N = 11).
	(*m) = 1
	for (*i) = 2; (*i) <= (2*(*(n)) - 1); (*i)++ {
		(*tm) = (*m)
		(*TI) = (*i)
		(*R) = (MOD((*tm), (*TI)))
		for (*R) != 0 {
			(*tm) = (*TI)
			(*TI) = (*R)
			(*R) = (MOD((*tm), (*TI)))
		}
		(*m) = ((*m) / (*TI)) * (*i)
	}
	//
	//     Generate the scaled Hilbert matrix in A
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(a))[(*i)-1][(*j)-1] = DBLE((*m)) / ((*i) + (*j) - 1)
		}
	}
	//
	//     Generate matrix B as simply the first nrhs columns of m * the
	//     identity.
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), (n), (nrhs), func() *float64 {y := 0.0e+0; return &y}(), DBLE((*m)), (b), (ldb))
	//     Generate the true solutions in X.  Because B = the first nrhs
	//     columns of M*i, the true solutions are just the first nrhs columns
	//     of the inverse Hilbert matrix.
	(*(work))[0] = (*(n))
	for (*j) = 2; (*j) <= (*(n)); (*j)++ {
		(*(work))[(*j)-1] = ((((*(work))[(*j)-0] / ((*j) - 1)) * ((*j) - 1 - (*(n)))) / ((*j) - 1)) * ((*(n)) + (*j) - 1)
	}
	//
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(x))[(*i)-1][(*j)-1] = ((*(work))[(*i)-1] * (*(work))[(*j)-1]) / ((*i) + (*j) - 1)
		}
	}
	//
}
