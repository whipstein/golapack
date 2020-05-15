package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dpst01 re_constructs a symmetric positive semidefinite matrix A
// from its L or U factors and the permutation matrix P and computes
// the residual
//    norm( P*l*l'*P' - A) / ( N * norm(a) * eps) or
//    norm( P*U'*U*P' - A) / ( N * norm(a) * eps),
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
//       SUBROUTinE Dpst01( uplo, n, a, lda, afac, ldafac, perm, ldperm,
//                          piv, rwork, resid, rank)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   resid
//       intEGER            lda, ldafac, ldperm, n, rank
//       CHARACTER          uplo
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *),
//      $                   perm( ldperm, *), rwork(*)
//       intEGER            piv(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dpst01 re_constructs a symmetric positive semidefinite matrix A
// from its L or U factors and the permutation matrix P and computes
// the residual
//    norm( P*l*l'*P' - A) / ( N * norm(a) * eps) or
//    norm( P*U'*U*P' - A) / ( N * norm(a) * eps),
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
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factor L or U from the L*l' or U'*U
//          factorization of A.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is intEGER
//          The leading dimension of the array afac.  ldafac >= max(1,N).
// \endverbatim
//
// \param[out] perm
// \verbatim
//          perm is DOUBLE PRECISION array, dimension (ldperm,N)
//          Overwritten with the re_constructed matrix, and then with the
//          difference P*l*l'*P' - A (or P*U'*U*P' - A)
// \endverbatim
//
// \param[in] ldperm
// \verbatim
//          ldperm is intEGER
//          The leading dimension of the array perm.
//          ldaperm >= max(1,N).
// \endverbatim
//
// \param[in] piv
// \verbatim
//          piv is intEGER array, dimension (n)
//          piv is such that the nonzero entries are
//          P( piv( K), K) = 1.
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
// \param[in] rank
// \verbatim
//          rank is intEGER
//          number of nonzero singular values of A.
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
func Dpst01(uplo *byte, n *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, perm *[][]float64, ldperm *int, piv *[]int, rwork *[]float64, resid *float64, rank *int) {
	zero := new(float64)
	one := new(float64)
	anorm := new(float64)
	eps := new(float64)
	t := new(float64)
	i := new(int)
	j := new(int)
	k := new(int)
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
		//
		if (*(rank)) < (*(n)) {
			for (*j) = (*(rank)) + 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*(rank)) + 1; (*i) <= (*j); (*i)++ {
					(*(afac))[(*i)-(1)][(*j)-(1)] = (*zero)
					//Label100:
				}
				//Label110:
			}
		}
		//
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//
			//           Compute the (K,k) element of the result.
			//
			(*t) = (*Ddot(K, &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }(), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }()))
			(*(afac))[(*k)-(1)][(*k)-(1)] = (*t)
			//
			//           Compute the rest of column K.
			//
			Dtrmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("Transpose"); return &y }(), func() *[]byte {y := []byte("Non-unit"); return &y }(), (*k)-1, (afac), (ldafac), &((*(afac))[0][(*k)-(1)]), func() *int {y := 1; return &y }())
			//
			//Label120:
		}
		//
		//     Compute the product L*l', overwriting L.
		//
	} else {
		//
		if (*(rank)) < (*(n)) {
			for (*j) = (*(rank)) + 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
					(*(afac))[(*i)-(1)][(*j)-(1)] = (*zero)
					//Label130:
				}
				//Label140:
			}
		}
		//
		for (*k) = (*(n)); (*k) <= 1; (*k) += -1 {
			//           Add a multiple of column K of the factor L to each of
			//           columns K+1 through N.
			//
			if (*k)+1 <= (*(n)) {
				Dsyr(func() *[]byte {y := []byte("Lower"); return &y }(), (*(n))-(*k), one, &((*(afac))[(*k)+0][(*k)-(1)]), func() *int {y := 1; return &y }(), &((*(afac))[(*k)+0][(*k)+0]), (ldafac))
			}
			//
			//           Scale column K by the diagonal element.
			//
			(*t) = (*(afac))[(*k)-(1)][(*k)-(1)]
			Dscal((*(n))-(*k)+1, t, &((*(afac))[(*k)-(1)][(*k)-(1)]), func() *int {y := 1; return &y }())
			//Label150:
		}
		//
	}
	//
	//        Form P*l*l'*P' or P*U'*U*P'
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				if (*(piv))[(*i)-(1)] <= (*(piv))[(*j)-(1)] {
					if (*i) <= (*j) {
						(*(perm))[(piv)(i)-(1)][(*(piv))[(*j)-(1)]-(1)] = (*(afac))[(*i)-(1)][(*j)-(1)]
					} else {
						(*(perm))[(piv)(i)-(1)][(*(piv))[(*j)-(1)]-(1)] = (*(afac))[(*j)-(1)][(*i)-(1)]
					}
				}
				//Label160:
			}
			//Label170:
		}
		//
		//
	} else {
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				if (*(piv))[(*i)-(1)] >= (*(piv))[(*j)-(1)] {
					if (*i) >= (*j) {
						(*(perm))[(piv)(i)-(1)][(*(piv))[(*j)-(1)]-(1)] = (*(afac))[(*i)-(1)][(*j)-(1)]
					} else {
						(*(perm))[(piv)(i)-(1)][(*(piv))[(*j)-(1)]-(1)] = (*(afac))[(*j)-(1)][(*i)-(1)]
					}
				}
				//Label180:
			}
			//Label190:
		}
		//
	}
	//
	//     Compute the difference  P*l*l'*P' - A (or P*U'*U*P' - A).
	//
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*j); (*i)++ {
				(*(perm))[(*i)-(1)][(*j)-(1)] = (*(perm))[(*i)-(1)][(*j)-(1)] - (*(a))[(*i)-(1)][(*j)-(1)]
				//Label200:
			}
			//Label210:
		}
	} else {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
				(*(perm))[(*i)-(1)][(*j)-(1)] = (*(perm))[(*i)-(1)][(*j)-(1)] - (*(a))[(*i)-(1)][(*j)-(1)]
				//Label220:
			}
			//Label230:
		}
	}
	//
	//     Compute norm( P*l*l'P - A) / ( N * norm(a) * eps), or
	//     ( P*U'*U*P' - A)/ ( N * norm(a) * eps).
	//
	(*(resid)) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), (uplo), (n), (perm), (ldafac), (rwork)))
	//
	(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	//
	return
	//
	//     End of Dpst01
	//
}
