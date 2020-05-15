package goblas

import 

// Dlagge generates a real general m by n matrix a, by pre- and post-
// multiplying a real diagonal matrix D with random orthogonal matrices:
// A = U*D*V. The lower and upper bandwidths may then be reduced to
// kl and ku by additional orthogonal transformations.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlagge( m, n, kl, ku, d, a, lda, iseed, work, info)
//
//       .. Scalar Arguments ..
//       intEGER            info, kl, ku, lda, m, N
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
// Dlagge generates a real general m by n matrix a, by pre- and post-
// multiplying a real diagonal matrix D with random orthogonal matrices:
// A = U*D*V. The lower and upper bandwidths may then be reduced to
// kl and ku by additional orthogonal transformations.
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
//          The number of columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is intEGER
//          The number of nonzero subdiagonals within the band of A.
//          0 <= kl <= M-1.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//          The number of nonzero superdiagonals within the band of A.
//          0 <= ku <= N-1.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (min(m,N))
//          The diagonal elements of the diagonal matrix D.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The generated m by n matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= M.
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
//          work is DOUBLE PRECISION array, dimension (M+N)
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
func Dlagge(m *int, n *int, kl *int, ku *int, d *[]float64, a *[][]float64, lda *int, iseed *[]int, work *[]float64, info *int) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
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
	if (*(m)) < 0 {
		(*(info)) = -1
	} else if (*(n)) < 0 {
		(*(info)) = -2
	} else if (*(kl)) < 0 || (*(kl)) > (*(m))-1 {
		(*(info)) = -3
	} else if (*(ku)) < 0 || (*(ku)) > (*(n))-1 {
		(*(info)) = -4
	} else if (*(lda)) < (MAX(1, (*(m)))) {
		(*(info)) = -7
	}
	if (*(info)) < 0 {
		Xerbla(func() *[]byte {y :=[]byte("Dlagge"); return &y }(), -(*(info)))
		return
	}
	//
	//     initialize A to diagonal matrix
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			(*(a))[(*i)-1][(*j)-1] = (*zero)
			//Label10:
		}
		//Label20:
	}
	for (*i) = 1; (*i) <= (Min((*(m)), (*(n)))); (*i)++ {
		(*(a))[(*i)-1][(*i)-1] = (*(d))[(*i)-1]
		//Label30:
	}
	//
	//     Quick exit if the user wants a diagonal matrix
	//
	if ((*(kl)) == 0) && ((*(ku)) == 0) {
		return
	}
	//
	//     pre- and post-multiply A by random orthogonal matrices
	//
	for (*i) = Min((*(m)), (*(n))); (*i) <= 1; (*i) += -1 {
		if (*i) < (*(m)) {
			//
			//           generate random reflection
			//
			Dlarnv(func() *int {y := 3; return &y }(), (iseed), (*(m))-(*i)+1, (work))
			(*wn) = (*Dnrm2((*(m))-(*i)+1, (work), func() *int {y := 1; return &y }()))
			(*wa) = (*SIGN(wn, &((*(work))[0])))
			if (*wn) == (*zero) {
				(*tau) = (*zero)
			} else {
				(*wb) = (*(work))[0] + (*wa)
				Dscal((*(m))-(*i), (*one)/(*wb), &((*(work))[1]), func() *int {y := 1; return &y }())
				(*(work))[0] = (*one)
				(*tau) = (*wb) / (*wa)
			}
			//
			//           multiply a(i:m,i:n) by random reflection from the left
			//
			Dgemv(func() *[]byte {y :=[]byte("Transpose"); return &y }(), (*(m))-(*i)+1, (*(n))-(*i)+1, one, &((*(a))[(*i)-1][(*i)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*(m))+0]), func() *int {y := 1; return &y }())
			Dger((*(m))-(*i)+1, (*(n))-(*i)+1, -(*tau), (work), func() *int {y := 1; return &y }(), &((*(work))[(*(m))+0]), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][(*i)-1]), (lda))
		}
		if (*i) < (*(n)) {
			//
			//           generate random reflection
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
			//           multiply a(i:m,i:n) by random reflection from the right
			//
			Dgemv(func() *[]byte {y :=[]byte("No transpose"); return &y }(), (*(m))-(*i)+1, (*(n))-(*i)+1, one, &((*(a))[(*i)-1][(*i)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }())
			Dger((*(m))-(*i)+1, (*(n))-(*i)+1, -(*tau), &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][(*i)-1]), (lda))
		}
		//Label40:
	}
	//
	//     Reduce number of subdiagonals to kl and number of superdiagonals
	//     to ku
	//
	for (*i) = 1; (*i) <= (MAX((*(m))-1-(*(kl)), (*(n))-1-(*(ku)))); (*i)++ {
		if (*(kl)) <= (*(ku)) {
			//
			//           annihilate subdiagonal elements first (necessary if kl = 0)
			//
			if (*i) <= (Min((*(m))-1-(*(kl)), (*(n)))) {
				//
				//              generate reflection to annihilate a(kl+i+1:m,i)
				//
				(*wn) = (*Dnrm2((*(m))-(*(kl))-(*i)+1, &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }()))
				(*wa) = (*SIGN(wn, &((*(a))[(*(kl))+(*i)-1][(*i)-1])))
				if (*wn) == (*zero) {
					(*tau) = (*zero)
				} else {
					(*wb) = (*(a))[(*(kl))+(*i)-1][(*i)-1] + (*wa)
					Dscal((*(m))-(*(kl))-(*i), (*one)/(*wb), &((*(a))[(*(kl))+(*i)+0][(*i)-1]), func() *int {y := 1; return &y }())
					(*(a))[(*(kl))+(*i)-1][(*i)-1] = (*one)
					(*tau) = (*wb) / (*wa)
				}
				//
				//              apply reflection to a(kl+i:m,i+1:n) from the left
				//
				Dgemv(func() *[]byte {y :=[]byte("Transpose"); return &y }(), (*(m))-(*(kl))-(*i)+1, (*(n))-(*i), one, &((*(a))[(*(kl))+(*i)-1][(*i)+0]), (lda), &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), zero, (work), func() *int {y := 1; return &y }())
				Dger((*(m))-(*(kl))-(*i)+1, (*(n))-(*i), -(*tau), &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*(kl))+(*i)-1][(*i)+0]), (lda))
				(*(a))[(*(kl))+(*i)-1][(*i)-1] = -(*wa)
			}
			//
			if (*i) <= (Min((*(n))-1-(*(ku)), (*(m)))) {
				//
				//              generate reflection to annihilate a(i,ku+i+1:n)
				//
				(*wn) = (*Dnrm2((*(n))-(*(ku))-(*i)+1, &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda)))
				(*wa) = (*SIGN(wn, &((*(a))[(*i)-1][(*(ku))+(*i)-1])))
				if (*wn) == (*zero) {
					(*tau) = (*zero)
				} else {
					(*wb) = (*(a))[(*i)-1][(*(ku))+(*i)-1] + (*wa)
					Dscal((*(n))-(*(ku))-(*i), (*one)/(*wb), &((*(a))[(*i)-1][(*(ku))+(*i)+0]), (lda))
					(*(a))[(*i)-1][(*(ku))+(*i)-1] = (*one)
					(*tau) = (*wb) / (*wa)
				}
				//
				//              apply reflection to a(i+1:m,ku+i:n) from the right
				//
				Dgemv(func() *[]byte {y :=[]byte("No transpose"); return &y }(), (*(m))-(*i), (*(n))-(*(ku))-(*i)+1, one, &((*(a))[(*i)+0][(*(ku))+(*i)-1]), (lda), &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda), zero, (work), func() *int {y := 1; return &y }())
				Dger((*(m))-(*i), (*(n))-(*(ku))-(*i)+1, -(*tau), (work), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda), &((*(a))[(*i)+0][(*(ku))+(*i)-1]), (lda))
				(*(a))[(*i)-1][(*(ku))+(*i)-1] = -(*wa)
			}
		} else {
			//
			//           annihilate superdiagonal elements first (necessary if
			//           ku = 0)
			//
			if (*i) <= (Min((*(n))-1-(*(ku)), (*(m)))) {
				//
				//              generate reflection to annihilate a(i,ku+i+1:n)
				//
				(*wn) = (*Dnrm2((*(n))-(*(ku))-(*i)+1, &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda)))
				(*wa) = (*SIGN(wn, &((*(a))[(*i)-1][(*(ku))+(*i)-1])))
				if (*wn) == (*zero) {
					(*tau) = (*zero)
				} else {
					(*wb) = (*(a))[(*i)-1][(*(ku))+(*i)-1] + (*wa)
					Dscal((*(n))-(*(ku))-(*i), (*one)/(*wb), &((*(a))[(*i)-1][(*(ku))+(*i)+0]), (lda))
					(*(a))[(*i)-1][(*(ku))+(*i)-1] = (*one)
					(*tau) = (*wb) / (*wa)
				}
				//
				//              apply reflection to a(i+1:m,ku+i:n) from the right
				//
				Dgemv(func() *[]byte {y :=[]byte("No transpose"); return &y }(), (*(m))-(*i), (*(n))-(*(ku))-(*i)+1, one, &((*(a))[(*i)+0][(*(ku))+(*i)-1]), (lda), &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda), zero, (work), func() *int {y := 1; return &y }())
				Dger((*(m))-(*i), (*(n))-(*(ku))-(*i)+1, -(*tau), (work), func() *int {y := 1; return &y }(), &((*(a))[(*i)-1][(*(ku))+(*i)-1]), (lda), &((*(a))[(*i)+0][(*(ku))+(*i)-1]), (lda))
				(*(a))[(*i)-1][(*(ku))+(*i)-1] = -(*wa)
			}
			//
			if (*i) <= (Min((*(m))-1-(*(kl)), (*(n)))) {
				//
				//              generate reflection to annihilate a(kl+i+1:m,i)
				//
				(*wn) = (*Dnrm2((*(m))-(*(kl))-(*i)+1, &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }()))
				(*wa) = (*SIGN(wn, &((*(a))[(*(kl))+(*i)-1][(*i)-1])))
				if (*wn) == (*zero) {
					(*tau) = (*zero)
				} else {
					(*wb) = (*(a))[(*(kl))+(*i)-1][(*i)-1] + (*wa)
					Dscal((*(m))-(*(kl))-(*i), (*one)/(*wb), &((*(a))[(*(kl))+(*i)+0][(*i)-1]), func() *int {y := 1; return &y }())
					(*(a))[(*(kl))+(*i)-1][(*i)-1] = (*one)
					(*tau) = (*wb) / (*wa)
				}
				//
				//              apply reflection to a(kl+i:m,i+1:n) from the left
				//
				Dgemv(func() *[]byte {y :=[]byte("Transpose"); return &y }(), (*(m))-(*(kl))-(*i)+1, (*(n))-(*i), one, &((*(a))[(*(kl))+(*i)-1][(*i)+0]), (lda), &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), zero, (work), func() *int {y := 1; return &y }())
				Dger((*(m))-(*(kl))-(*i)+1, (*(n))-(*i), -(*tau), &((*(a))[(*(kl))+(*i)-1][(*i)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*(kl))+(*i)-1][(*i)+0]), (lda))
				(*(a))[(*(kl))+(*i)-1][(*i)-1] = -(*wa)
			}
		}
		//
		if (*i) <= (*(n)) {
			for (*j) = (*(kl)) + (*i) + 1; (*j) <= (*(m)); (*j)++ {
				(*(a))[(*j)-1][(*i)-1] = (*zero)
				//Label50:
			}
		}
		//
		if (*i) <= (*(m)) {
			for (*j) = (*(ku)) + (*i) + 1; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*i)-1][(*j)-1] = (*zero)
				//Label60:
			}
		}
		//Label70:
	}
	return
	//
	//     End of Dlagge
	//
}
