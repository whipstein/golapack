package goblas

import 

// Dlarge pre- and post-multiplies a real general n by n matrix A
// with a random orthogonal matrix: A = U*D*U'.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlarge( n, a, lda, iseed, work, info)
//
//       .. Scalar Arguments ..
//       intEGER            info, lda, N
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlarge pre- and post-multiplies a real general n by n matrix A
// with a random orthogonal matrix: A = U*D*U'.
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
// \param[in,out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          On entry, the original n by n matrix A.
//          On exit, A is overwritten by U*a*U' for some random
//          orthogonal matrix U.
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
func Dlarge(n *int, a *[][]float64, lda *int, iseed *[]int, work *[]float64, info *int) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
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
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input arguments
	//
	(*(info)) = 0
	if (*(n)) < 0 {
		(*(info)) = -1
	} else if (*(lda)) < (MAX(1, (*(n)))) {
		(*(info)) = -3
	}
	if (*(info)) < 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dlarge"); return &y }(), -(*(info)))
		return
	}
	//
	//     pre- and post-multiply A by random orthogonal matrix
	//
	for (*i) = (*(n)); (*i) <= 1; (*i) += -1 {
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
		//        multiply a(i:n,1:n) by random reflection from the left
		//
		Dgemv(func() *[]byte {y :=[]byte("Transpose"); return &y }(), (*(n))-(*i)+1, (n), one, &((*(a))[(*i)-1][0]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }())
		Dger((*(n))-(*i)+1, (n), -(*tau), (work), func() *int {y := 1; return &y }(), &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][0]), (lda))
		//
		//        multiply a(1:n,i:n) by random reflection from the right
		//
		Dgemv(func() *[]byte {y :=[]byte("No transpose"); return &y }(), (n), (*(n))-(*i)+1, one, &((*(a))[0][(*i)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }())
		Dger((n), (*(n))-(*i)+1, -(*tau), &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[0][(*i)-1]), (lda))
		//Label10:
	}
	return
	//
	//     End of Dlarge
	//
}
