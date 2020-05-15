package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dpot01 re_constructs a symmetric positive definite matrix  A  from
// its L*l' or U'*U factorization and computes the residual
//    norm( L*l' - A) / ( N * norm(a) * eps) or
//    norm( U'*U - A) / ( N * norm(a) * eps),
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
//       SUBROUTinE Dpot01( uplo, n, a, lda, afac, ldafac, rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            lda, ldafac, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *), rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dpot01 re_constructs a symmetric positive definite matrix  A  from
// its L*l' or U'*U factorization and computes the residual
//    norm( L*l' - A) / ( N * norm(a) * eps) or
//    norm( U'*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the upper or lower triangular part of the
//          symmetric matrix A is stored:
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original symmetric matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,N)
// \endverbatim
//
// \param[in,out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          On entry, the factor L or U from the L*l' or U'*U
//          factorization of A.
//          Overwritten with the re_constructed matrix, and then with the
//          difference L*l' - A (or U'*U - A).
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is intEGER
//          The leading dimension of the array afac.  ldafac >= max(1,N).
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          If uplo = 'L', norm(L*l' - A) / ( N * norm(a) * eps)
//          If uplo = 'U', norm(U'*U - A) / ( N * norm(a) * eps)
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
func Dpot01(uplo *byte, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, rwork *[]float64, resid *float64) {
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
	//     Quick exit if N = 0.
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (lda), (rwork)))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute the product U'*U, overwriting U.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//
			//           Compute the (K,k) element of the result.
			//
			(*t) = (*Ddot(K, &((*(afac))[0][(*k)-1]), func() *int {y := 1; return &y }(), &((*(afac))[0][(*k)-1]), func() *int {y := 1; return &y }()))
			(*(afac))[(*k)-1][(*k)-1] = (*t)
			//
			//           Compute the rest of column K.
			//
			Dtrmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (*k)-1, (afac), (ldafac), &((*(afac))[0][(*k)-1]), func() *int {y := 1; return &y }())
			//
			//Label10:
		}
		//
		//     Compute the product L*l', overwriting L.
		//
	} else {
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//
			//           Add a multiple of column K of the factor L to each of
			//           columns K+1 through N.
			//
			if (*k)+1 <= (*(n)) {
				Dsyr(func() *[]byte {y := []byte("Lower"); return &y }(), (*(n))-(*k), one, &((*(afac))[(*k)+0][(*k)-1]), func() *int {y := 1; return &y }(), &((*(afac))[(*k)+0][(*k)+0]), (ldafac))
			}
			//
			//           Scale column K by the diagonal element.
			//
			(*t) = (*(afac))[(*k)-1][(*k)-1]
			Dscal((*(n))-(*k)+1, t, &((*(afac))[(*k)-1][(*k)-1]), func() *int {y := 1; return &y }())
			//
			//Label20:
		}
	}
	//
	//     Compute the difference  L*l' - A (or U'*U - A).
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j); (*i)++ {
				(*(afac))[(*i)-1][(*j)-1] = (*(afac))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
				//Label30:
			}
			//Label40:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
				(*(afac))[(*i)-1][(*j)-1] = (*(afac))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
				//Label50:
			}
			//Label60:
		}
	}
	//
	//     Compute norm( L*U - A) / ( N * norm(a) * eps)
	//
	(*(resid)) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (afac), (ldafac), (rwork)))
	//
	(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	//
	return
	//
	//     End of Dpot01
	//
}
