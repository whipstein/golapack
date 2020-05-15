package goblas

import 

//  Dtrt06 computes a test ratio comparing rcond (the reciprocal
//  condition number of a triangular matrix A) and rcondc, the estimate
//  computed by Dtrcon.  Information about the triangular matrix A is
//  used if one estimate is zero and the other is non-zero to decide if
//  underflow in the estimate is justified.
//
//   =========== DOCUMENTATION ===========
//
//  Online html documentation available at
//             http://www.netlib.org/lapack/explore-html/
//
//   Definition:
//   ===========
//
//        SUBROUTinE Dtrt06( rcond, rcondc, uplo, diag, n, a, lda, work,
//                           rat)
//
//        .. Scalar Arguments ..
//        CHARACTER          diag, uplo
//        intEGER            lda, N
//        DOUBLE PRECISION   rat, rcond, rcondc
//        ..
//        .. Array Arguments ..
//        DOUBLE PRECISION   a( lda, *), work(*)
//        ..
//
//
//  \par Purpose:
//   =============
//
//  \verbatim
//
//  DTRT06 computes a test ratio comparing rcond (the reciprocal
//  condition number of a triangular matrix A) and rcondc, the estimate
//  computed by Dtrcon.  Information about the triangular matrix A is
//  used if one estimate is zero and the other is non-zero to decide if
//  underflow in the estimate is justified.
//  \endverbatim
//
//   Arguments:
//   ==========
//
//  \param[in] rcond
//  \verbatim
//           rcond is DOUBLE PRECISION
//           The estimate of the reciprocal condition number obtained by
//           forming the explicit inverse of the matrix A and computing
//           rcond = 1/( norm(a) * norm(inv(a))).
//  \endverbatim
//
//  \param[in] rcondc
//  \verbatim
//           rcondc is DOUBLE PRECISION
//           The estimate of the reciprocal condition number computed by
//           Dtrcon.
//  \endverbatim
//
//  \param[in] uplo
//  \verbatim
//           uplo is CHARACTER
//           Specifies whether the matrix A is upper or lower triangular.
//           = 'U':  Upper triangular
//           = 'L':  Lower triangular
//  \endverbatim
//
//  \param[in] diag
//  \verbatim
//           diag is CHARACTER
//           Specifies whether or not the matrix A is unit triangular.
//           = 'N':  Non-unit triangular
//           = 'U':  Unit triangular
//  \endverbatim
//
//  \param[in] N
//  \verbatim
//           N is intEGER
//           The order of the matrix A.  N >= 0.
//  \endverbatim
//
//  \param[in] A
//  \verbatim
//           A is DOUBLE PRECISION array, dimension (lda,N)
//           The triangular matrix A.  If uplo = 'U', the leading n by n
//           upper triangular part of the array A contains the upper
//           triangular matrix, and the strictly lower triangular part of
//           A is not referenced.  If uplo = 'L', the leading n by n lower
//           triangular part of the array A contains the lower triangular
//           matrix, and the strictly upper triangular part of A is not
//           referenced.  If diag = 'U', the diagonal elements of A are
//           also not referenced and are assumed to be 1.
//  \endverbatim
//
//  \param[in] lda
//  \verbatim
//           lda is intEGER
//           The leading dimension of the array A.  lda >= max(1,N).
//  \endverbatim
//
//  \param[out] work
//  \verbatim
//           work is DOUBLE PRECISION array, dimension (n)
//  \endverbatim
//
//  \param[out] rat
//  \verbatim
//           rat is DOUBLE PRECISION
//           The test ratio.  If both rcond and rcondc are nonzero,
//              rat = MAX( rcond, rcondc)/Min( rcond, rcondc) - 1.
//           If rat = 0, the two estimates are exactly the same.
//  \endverbatim
//
//   Authors:
//   ========
//
//  \author Univ. of Tennessee
//  \author Univ. of California Berkeley
//  \author Univ. of Colorado Denver
//  \author NAG Ltd.
//
//  \date December 2016
//
//  \ingroup double_lin
//
//   =====================================================================
func Dtrt06(rcond *float64, rcondc *float64, uplo *byte, diag *byte, n *int, a *[][]float64, lda *int, work *[]float64, rat *float64) {
	zero := new(float64)
	one := new(float64)
	anorm := new(float64)
	bignum := new(float64)
	eps := new(float64)
	rmax := new(float64)
	rmin := new(float64)
	smlnum := new(float64)
	//
	// -- lapACK test routine (version 3.7.0) --
	// -- lapACK is a software package provided by Univ. of Tennessee,    --
	// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//    December 2016
	//
	//    .. Scalar Arguments ..
	//    ..
	//    .. Array Arguments ..
	//    ..
	//
	// =====================================================================
	//
	//    .. Parameters ..
	(*zero) = 0.0e+0
	(*one) = 1.0e+0
	//  ..
	//  .. Local Scalars ..
	//  ..
	//  .. External Functions ..
	//  ..
	//  .. Intrinsic Functions ..
	//  ..
	//  .. External Subroutines ..
	//  ..
	//  .. Executable Statements ..
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	(*rmax) = (MAX((*(rcond)), (*(rcondc))))
	(*rmin) = (Min((*(rcond)), (*(rcondc))))
	//
	//  Do the easy cases first.
	//
	if (*rmin) < (*zero) {
		//
		// Invalid value for rcond or rcondc, return 1/eps.
		//
		(*(rat)) = (*one) / (*eps)

	} else if (*rmin) > (*zero) {
		//
		// Both estimates are positive, return rmax/rmin - 1.
		//
		(*(rat)) = (*rmax)/(*rmin) - (*one)

	} else if (*rmax) == (*zero) {
		//
		// Both estimates zero.
		//
		(*(rat)) = (*zero)

	} else {
		//
		// one estimate is zero, the other is non-zero.  If the matrix is
		// ill-conditioned, return the nonzero estimate multiplied by
		// 1/eps; if the matrix is badly scaled, return the nonzero
		// estimate multiplied by bignum/tmAX, where tmAX is the maximum
		// element in absolute value in A.
		//
		(*smlnum) = (*Dlamch(func() *[]byte {y :=[]byte("Safe minimum"); return &y}()))
		(*bignum) = (*one) / (*smlnum)
		Dlabad(smlnum, bignum)
		(*anorm) = (*Dlantr(func() *byte {y := byte('M'); return &y}(), (uplo), (diag), (n), (n), (a), (lda), (work)))

		(*(rat)) = (*rmax) * (Min((*bignum)/MAX((*one), (*anorm)), (*one)/(*eps)))
	}

	return
	//
	//  End of DTRT06
	//
}
