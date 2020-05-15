package goblas

import 

// Dget01 re_constructs a matrix A from its L*U factorization and
// computes the residual
//    norm(L*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dget01( m, n, a, lda, afac, ldafac, ipiv, rwork,
//                          resid)
//
//       .. Scalar Arguments ..
//       intEGER            lda, ldafac, m, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       intEGER            ipiv(*)
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *), rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dget01 re_constructs a matrix A from its L*U factorization and
// computes the residual
//    norm(L*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon.
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
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original M x N matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,M).
// \endverbatim
//
// \param[in,out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factored form of the matrix A.  afac contains the factors
//          L and U from the L*U factorization as computed by DGETRF.
//          Overwritten with the re_constructed matrix, and then with the
//          difference L*U - A.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is intEGER
//          The leading dimension of the array afac.  ldafac >= max(1,M).
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is intEGER array, dimension (n)
//          The pivot indices from DGETRF.
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (m)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(L*U - A) / ( N * norm(a) * eps)
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
func Dget01(m *int, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, ipiv *[]int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
	k := new(int)
	anorm := new(float64)
	eps := new(float64)
	t := new(float64)
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
	//
	//     .. Parameters ..
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Quick exit if M = 0 or N = 0.
	//
	if (*(m)) <= 0 || (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Determine eps and the norm of A.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*anorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (a), (lda), (rwork)))
	//
	//     Compute the product L*U and overwrite afac with the result.
	//     A column at a time of the product is obtained, starting with
	//     column N.
	//
	for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
		if (*k) > (*(m)) {
			Dtrmv(func() *[]byte {y :=[]byte("Lower"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Unit"); return &y }(), (m), (afac), (ldafac), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }())
		} else {
			//
			//           Compute elements (K+1:M,k)
			//
			(*t) = (*(afac))[(*k)-(1)][(*k)-(1)]
			if (*k)+1 <= (*(m)) {
				Dscal((*(m))-(*k), t, &((*(afac))[(*k)+0][(*k)-(1)]), func() *int {y := 1; return &y }())
				Dgemv(func() *[]byte {y :=[]byte("No transpose"); return &y }(), (*(m))-(*k), (*k)-1, one, &((*(afac))[(*k)+0][0]), (ldafac), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }(), one, &((*(afac))[(*k)+0][(*k)-(1)]), func() *int {y := 1; return &y }())
			}
			//
			//           Compute the (K,k) element
			//
			(*(afac))[(*k)-(1)][(*k)-(1)] = (*t) + Ddot((*k)-1, &((*(afac))[(*k)-(1)][0]), (ldafac), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }())
			//
			//           Compute elements (1:K-1,k)
			//
			Dtrmv(func() *[]byte {y :=[]byte("Lower"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *[]byte {y :=[]byte("Unit"); return &y }(), (*k)-1, (afac), (ldafac), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }())
		}
		//Label10:
	}
	Dlaswp((n), (afac), (ldafac), func() *int {y := 1; return &y }(), Min((*(m)), (*(n))), (ipiv), -1)
	//
	//     Compute the difference  L*U - A  and store in afac.
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			(*(afac))[(*i)-(1)][(*j)-(1)] = (*(afac))[(*i)-(1)][(*j)-(1)] - (*(a))[(*i)-(1)][(*j)-(1)]
			//Label20:
		}
		//Label30:
	}
	//
	//     Compute norm( L*U - A) / ( N * norm(a) * eps)
	//
	(*(resid)) = (*Dlange(func() *byte {y := byte('1'); return &y }(), (m), (n), (afac), (ldafac), (rwork)))
	//
	if (*anorm) <= (*zero) {
		if (*(resid)) != (*zero) {
			(*(resid)) = (*one) / (*eps)
		}
	} else {
		(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	}
	//
	return
	//
	//     End of Dget01
	//
}
