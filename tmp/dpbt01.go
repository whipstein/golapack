package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dpbt01 re_constructs a symmetric positive definite band matrix A from
// its L*L' or U'*U factorization and computes the residual
//    norm( L*L' - A) / ( N * norm(a) * eps) or
//    norm( U'*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon, L' is the conjugate transpose of
// L, and U' is the conjugate transpose of U.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dpbt01( uplo, n, kd, a, lda, afac, ldafac, rwork,
//                          resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       inTEGER            kd, lda, ldafac, N
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
// Dpbt01 re_constructs a symmetric positive definite band matrix A from
// its L*L' or U'*U factorization and computes the residual
//    norm( L*L' - A) / ( N * norm(a) * eps) or
//    norm( U'*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon, L' is the conjugate transpose of
// L, and U' is the conjugate transpose of U.
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
//          N is inTEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kd
// \verbatim
//          kd is inTEGER
//          The number of super-diagonals of the matrix A if uplo = 'U',
//          or the number of sub-diagonals if uplo = 'L'.  kd >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original symmetric band matrix A.  If uplo = 'U', the
//          upper triangular part of A is stored as a band matrix; if
//          uplo = 'L', the lower triangular part of A is stored.  The
//          columns of the appropriate triangle are stored in the columns
//          of A and the diagonals of the triangle are stored in the rows
//          of A.  See DPBTRF for further details.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER.
//          The leading dimension of the array A.  lda >= max(1,kd+1).
// \endverbatim
//
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factored form of the matrix A.  afac contains the factor
//          L or U from the L*L' or U'*U factorization in band storage
//          format, as computed by DPBTRF.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is inTEGER
//          The leading dimension of the array afac.
//          ldafac >= max(1,kd+1).
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
//          If uplo = 'L', norm(L*L' - A) / ( N * norm(a) * eps)
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
func Dpbt01(uplo *byte, n *int, kd *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	j := new(int)
	k := new(int)
	kc := new(int)
	klen := new(int)
	ml := new(int)
	mu := new(int)
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
	(*anorm) = (*Dlansb(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (kd), (a), (lda), (rwork)))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute the product U'*U, overwriting U.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			(*kc) = (MAX(1, (*(kd))+2-(*k)))
			(*klen) = (*(kd)) + 1 - (*kc)
			//
			//           Compute the (K,k) element of the result.
			//
			(*t) = (*Ddot((*klen)+1, &((*(afac))[(*kc)-1][(*k)-1]), func() *int {y := 1; return &y }(), &((*(afac))[(*kc)-1][(*k)-1]), func() *int {y := 1; return &y }()))
			(*(afac))[(*(kd))+0][(*k)-1] = (*t)
			//
			//           Compute the rest of column K.
			//
			if (*klen) > 0 {
				Dtrmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), klen, &((*(afac))[(*(kd))+0][(*k)-(*klen)-1]), (*(ldafac))-1, &((*(afac))[(*kc)-1][(*k)-1]), func() *int {y := 1; return &y }())
			}
			//
			//Label10:
		}
		//
		//     uplo = 'L':  Compute the product L*L', overwriting L.
		//
	} else {
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			(*klen) = (Min((*(kd)), (*(n))-(*k)))
			//
			//           Add a multiple of column K of the factor L to each of
			//           columns K+1 through N.
			//
			if (*klen) > 0 {
				Dsyr(func() *[]byte {y := []byte("Lower"); return &y }(), klen, one, &((*(afac))[1][(*k)-1]), func() *int {y := 1; return &y }(), &((*(afac))[0][(*k)+0]), (*(ldafac))-1)
			}
			//
			//           Scale column K by the diagonal element.
			//
			(*t) = (*(afac))[0][(*k)-1]
			Dscal((*klen)+1, t, &((*(afac))[0][(*k)-1]), func() *int {y := 1; return &y }())
			//
			//Label20:
		}
	}
	//
	//     Compute the difference  L*L' - A  or  U'*U - A.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			(*mu) = (MAX(1, (*(kd))+2-(*j)))
			for (*i) = (*mu); (*i) <= (*(kd))+1; (*i)++ {
				(*(afac))[(*i)-1][(*j)-1] = (*(afac))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
				//Label30:
			}
			//Label40:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			(*ML) = (Min((*(kd))+1, (*(n))-(*j)+1))
			for (*i) = 1; (*i) <= (*ML); (*i)++ {
				(*(afac))[(*i)-1][(*j)-1] = (*(afac))[(*i)-1][(*j)-1] - (*(a))[(*i)-1][(*j)-1]
				//Label50:
			}
			//Label60:
		}
	}
	//
	//     Compute norm( L*L' - A) / ( N * norm(a) * eps)
	//
	(*(resid)) = (*Dlansb(func() *byte {y := byte('I'); return &y }(), (uplo), (n), (kd), (afac), (ldafac), (rwork)))
	//
	(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	//
	return
	//
	//     End of Dpbt01
	//
}
