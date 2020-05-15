package goblas

import 

// Dlagsy generates a real symmetric matrix a, by pre- and post-
// multiplying a real diagonal matrix D with a random orthogonal matrix:
// A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
// orthogonal transformations.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlagsy( n, k, d, a, lda, iseed, work, info)
//
//       .. Scalar Arguments ..
//       intEGER            info, k, lda, N
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), d(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlagsy generates a real symmetric matrix a, by pre- and post-
// multiplying a real diagonal matrix D with a random orthogonal matrix:
// A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
// orthogonal transformations.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is intEGER
//          The number of nonzero subdiagonals within the band of A.
//          0 <= K <= N-1.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (n)
//          The diagonal elements of the diagonal matrix D.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The generated n by n symmetric matrix A (the full matrix is
//          stored).
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= N.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is intEGER array, dimension (4)
//          On entry, the seed of the random number generator; the array
//          elements must be between 0 and 4095, and iseed(4) must be
//          odd.
//          On exit, the seed is updated.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (2*n)
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
// \ingroup double_matgen
//
//  =====================================================================
func Dlagsy(n *int, k *int, d *[]float64, a *[][]float64, lda *int, iseed *[]int, work *[]float64, info *int) {
	zero := new(float64)
	one := new(float64)
	half := new(float64)
	i := new(int)
	j := new(int)
	alpha := new(float64)
	tau := new(float64)
	wa := new(float64)
	wb := new(float64)
	wn := new(float64)
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
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
	(*half) = 0.5e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input arguments
	//
	(*(info)) = 0
	if (*(n)) < 0 {
		(*(info)) = -1
	} else if (*(k)) < 0 || (*(k)) > (*(n))-1 {
		(*(info)) = -2
	} else if (*(lda)) < (MAX(1, (*(n)))) {
		(*(info)) = -5
	}
	if (*(info)) < 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dlagsy"); return &y }(), -(*(info)))
		return
	}
	//
	//     initialize lower triangle of A to diagonal matrix
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
			(*(a))[(*i)-1][(*j)-1] = (*zero)
			//Label10:
		}
		//Label20:
	}
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		(*(a))[(*i)-1][(*i)-1] = (*(d))[(*i)-1]
		//Label30:
	}
	//
	//     Generate lower triangle of symmetric matrix
	//
	for (*i) = (*(n)) - 1; (*i) <= 1; (*i) += -1 {
		//
		//        generate random reflection
		//
		Dlarnv(func() *int {y := 3; return &y }(), (iseed), (*(n))-(*i)+1, (work))
		(*wn) = (*Dnrm2((*(n))-(*i)+1, (work), func() *int {y := 1; return &y }()))
		(*wa) = (*SIGN(wn, &((*(work))[0])))
		if (*wn) == (*zero) {
			(*tau) = (*zero)
		} else {
			(*wb) = (*(work))[0] + (*wa)
			Dscal((*(n))-(*i), (*one)/(*wb), &((*(work))[1]), func() *int {y := 1; return &y }())
			(*(work))[0] = (*one)
			(*tau) = (*wb) / (*wa)
		}
		//
		//        apply random reflection to a(i:n,i:n) from the left
		//        and the right
		//
		//        compute  y := tau * A * u
		//
		Dsymv(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*(n))-(*i)+1, tau, &((*(a))[(*i)-1][(*i)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }())
		//
		//        compute  v := y - 1/2 * tau * ( y, u) * u
		//
		(*alpha) = -(*half) * (*tau) * Ddot((*(n))-(*i)+1, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		Daxpy((*(n))-(*i)+1, alpha, (work), func() *int {y := 1; return &y }(), &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }())
		//
		//        apply the transformation as a rank-2 update to a(i:n,i:n)
		//
		Dsyr2(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*(n))-(*i)+1, -(*one), (work), func() *int {y := 1; return &y }(), &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][(*i)-1]), (lda))
		//Label40:
	}
	//
	//     Reduce number of subdiagonals to K
	//
	for (*i) = 1; (*i) <= (*(n))-1-(*(k)); (*i)++ {
		//
		//        generate reflection to annihilate a(k+i+1:n,i)
		//
		(*wn) = (*Dnrm2((*(n))-(*(k))-(*i)+1, &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }()))
		(*wa) = (*SIGN(wn, &((*(a))[(*(k))+(*i)-1][(*i)-1])))
		if (*wn) == (*zero) {
			(*tau) = (*zero)
		} else {
			(*wb) = (*(a))[(*(k))+(*i)-1][(*i)-1] + (*wa)
			Dscal((*(n))-(*(k))-(*i), (*one)/(*wb), &((*(a))[(*(k))+(*i)+0][(*i)-1]), func() *int {y := 1; return &y }())
			(*(a))[(*(k))+(*i)-1][(*i)-1] = (*one)
			(*tau) = (*wb) / (*wa)
		}
		//
		//        apply reflection to a(k+i:n,i+1:k+i-1) from the left
		//
		Dgemv(func() *[]byte {y :=[]byte("Transpose"); return &y }(), (*(n))-(*(k))-(*i)+1, (*(k))-1, one, &((*(a))[(*(k))+(*i)-1][(*i)+0]), (lda), &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), zero, (work), func() *int {y := 1; return &y }())
		Dger((*(n))-(*(k))-(*i)+1, (*(k))-1, -(*tau), &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*(k))+(*i)-1][(*i)+0]), (lda))
		//
		//        apply reflection to a(k+i:n,k+i:n) from the left and the right
		//
		//        compute  y := tau * A * u
		//
		Dsymv(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*(n))-(*(k))-(*i)+1, tau, &((*(a))[(*(k))+(*i)-1][(*(k))+(*i)-1]), (lda), &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), zero, (work), func() *int {y := 1; return &y }())
		//
		//        compute  v := y - 1/2 * tau * ( y, u) * u
		//
		(*alpha) = -(*half) * (*tau) * Ddot((*(n))-(*(k))-(*i)+1, (work), func() *int {y := 1; return &y }(), &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }())
		Daxpy((*(n))-(*(k))-(*i)+1, alpha, &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
		//
		//        apply symmetric rank-2 update to a(k+i:n,k+i:n)
		//
		Dsyr2(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*(n))-(*(k))-(*i)+1, -(*one), &((*(a))[(*(k))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*(k))+(*i)-1][(*(k))+(*i)-1]), (lda))
		//
		(*(a))[(*(k))+(*i)-1][(*i)-1] = -(*wa)
		for (*j) = (*(k)) + (*i) + 1; (*j) <= (*(n)); (*j)++ {
			(*(a))[(*j)-1][(*i)-1] = (*zero)
			//Label50:
		}
		//Label60:
	}
	//
	//     Store full symmetric matrix
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
			(*(a))[(*j)-1][(*i)-1] = (*(a))[(*i)-1][(*j)-1]
			//Label70:
		}
		//Label80:
	}
	return
	//
	//     End of Dlagsy
	//
}
