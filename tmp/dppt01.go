package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dppt01 re_constructs a symmetric positive definite packed matrix A
// from its L*l' or U'*U factorization and computes the residual
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
//       SUBROUTinE Dppt01( uplo, n, a, afac, rwork, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          uplo
//       intEGER            N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a(*), afac(*), rwork(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dppt01 re_constructs a symmetric positive definite packed matrix A
// from its L*l' or U'*U factorization and computes the residual
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
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The original symmetric matrix a, stored as a packed
//          triangular matrix.
// \endverbatim
//
// \param[in,out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          On entry, the factor L or U from the L*l' or U'*U
//          factorization of a, stored as a packed triangular matrix.
//          Overwritten with the re_constructed matrix, and then with the
//          difference L*l' - A (or U'*U - A).
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
func Dppt01(uplo *byte, n *int, a *[]float64, afac *[]float64, rwork *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	k := new(int)
	kc := new(int)
	npp := new(int)
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
	//     Quick exit if N = 0
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (a), (rwork)))
	if (*anorm) <= (*zero) {
		(*(resid)) = (*one) / (*eps)
		return
	}
	//
	//     Compute the product U'*U, overwriting U.
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		(*kc) = ((*(n))*((*(n))-1))/2 + 1
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//
			//           Compute the (K,k) element of the result.
			//
			(*t) = (*Ddot(K, &((*(afac))[(*kc)-1]), func() *int {y := 1; return &y }(), &((*(afac))[(*kc)-1]), func() *int {y := 1; return &y }()))
			(*(afac))[(*kc)+(*k)-0] = (*t)
			//
			//           Compute the rest of column K.
			//
			if (*k) > 1 {
				Dtpmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (*k)-1, (afac), &((*(afac))[(*kc)-1]), func() *int {y := 1; return &y }())
				(*kc) = (*kc) - ((*k) - 1)
			}
			//Label10:
		}
		//
		//     Compute the product L*l', overwriting L.
		//
	} else {
		(*kc) = ((*(n)) * ((*(n)) + 1)) / 2
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//
			//           Add a multiple of column K of the factor L to each of
			//           columns K+1 through N.
			//
			if (*k) < (*(n)) {
				Dspr(func() *[]byte {y := []byte("Lower"); return &y }(), (*(n))-(*k), one, &((*(afac))[(*kc)+0]), func() *int {y := 1; return &y }(), &((*(afac))[(*kc)+(*(n))-(*k)+0]))
			}
			//
			//           Scale column K by the diagonal element.
			//
			(*t) = (*(afac))[(*kc)-1]
			Dscal((*(n))-(*k)+1, t, &((*(afac))[(*kc)-1]), func() *int {y := 1; return &y }())
			//
			(*kc) = (*kc) - ((*(n)) - (*k) + 2)
			//Label20:
		}
	}
	//
	//     Compute the difference  L*l' - A (or U'*U - A).
	//
	(*npp) = (*(n)) * ((*(n)) + 1) / 2
	for (*i) = 1; (*i) <= (*npp); (*i)++ {
		(*(afac))[(*i)-1] = (*(afac))[(*i)-1] - (*(a))[(*i)-1]
		//Label30:
	}
	//
	//     Compute norm( L*U - A) / ( N * norm(a) * eps)
	//
	(*(resid)) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (afac), (rwork)))
	//
	(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	//
	return
	//
	//     End of Dppt01
	//
}
