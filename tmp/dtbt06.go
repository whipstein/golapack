package goblas

import 

// Dtbt06 computes a test ratio comparing rcond (the reciprocal
// condition number of a triangular matrix A) and rcondc, the estimate
// computed by Dtbcon.  Information about the triangular matrix A is
// used if one estimate is zero and the other is non-zero to decide if
// underflow in the estimate is justified.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dtbt06( rcond, rcondc, uplo, diag, n, kd, AB, ldab,
//                          work, rat)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, uplo
//       intEGER            kd, ldab, N
//       DOUBLE PRECISION   rat, rcond, rcondc
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   AB( ldab, *), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtbt06 computes a test ratio comparing rcond (the reciprocal
// condition number of a triangular matrix A) and rcondc, the estimate
// computed by Dtbcon.  Information about the triangular matrix A is
// used if one estimate is zero and the other is non-zero to decide if
// underflow in the estimate is justified.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The estimate of the reciprocal condition number obtained by
//          forming the explicit inverse of the matrix A and computing
//          rcond = 1/( norm(a) * norm(inv(a))).
// \endverbatim
//
// \param[in] rcondc
// \verbatim
//          rcondc is DOUBLE PRECISION
//          The estimate of the reciprocal condition number computed by
//          Dtbcon.
// \endverbatim
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER
//          Specifies whether the matrix A is upper or lower triangular.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER
//          Specifies whether or not the matrix A is unit triangular.
//          = 'N':  Non-unit triangular
//          = 'U':  Unit triangular
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kd
// \verbatim
//          kd is intEGER
//          The number of superdiagonals or subdiagonals of the
//          triangular band matrix A.  kd >= 0.
// \endverbatim
//
// \param[in] AB
// \verbatim
//          AB is DOUBLE PRECISION array, dimension (ldab,N)
//          The upper or lower triangular band matrix a, stored in the
//          first kd+1 rows of the array. The j-th column of A is stored
//          in the j-th column of the array AB as follows:
//          if uplo = 'U', AB(kd+1+i-j,j) = a(i,j) for max(1,j-kd)<=i<=j;
//          if uplo = 'L', AB(1+i-j,j)    = a(i,j) for j<=i<=min(n,j+kd).
// \endverbatim
//
// \param[in] ldab
// \verbatim
//          ldab is intEGER
//          The leading dimension of the array AB.  ldab >= kd+1.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] rat
// \verbatim
//          rat is DOUBLE PRECISION
//          The test ratio.  If both rcond and rcondc are nonzero,
//             rat = MAX( rcond, rcondc)/Min( rcond, rcondc) - 1.
//          If rat = 0, the two estimates are exactly the same.
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
func Dtbt06(rcond *float64, rcondc *float64, uplo *byte, diag *byte, n *int, kd *int, ab *[][]float64, ldab *int, work *[]float64, rat *float64) {
	zero := new(float64)
	one := new(float64)
	anorm := new(float64)
	bignum := new(float64)
	eps := new(float64)
	rmax := new(float64)
	rmin := new(float64)
	smlnum := new(float64)
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
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*rmax) = (MAX((*(rcond)), (*(rcondc))))
	(*rmin) = (Min((*(rcond)), (*(rcondc))))
	//
	//     Do the easy cases first.
	//
	if (*rmin) < (*zero) {
		//
		//        Invalid value for rcond or rcondc, return 1/eps.
		//
		(*(rat)) = (*one) / (*eps)
		//
	} else if (*rmin) > (*zero) {
		//
		//        Both estimates are positive, return rmax/rmin - 1.
		//
		(*(rat)) = (*rmax)/(*rmin) - (*one)
		//
	} else if (*rmax) == (*zero) {
		//
		//        Both estimates zero.
		//
		(*(rat)) = (*zero)
		//
	} else {
		//
		//        one estimate is zero, the other is non-zero.  If the matrix is
		//        ill-conditioned, return the nonzero estimate multiplied by
		//        1/eps; if the matrix is badly scaled, return the nonzero
		//        estimate multiplied by bignum/tmAX, where tmAX is the maximum
		//        element in absolute value in A.
		//
		(*smlnum) = (*Dlamch(func() *[]byte {y :=[]byte("Safe minimum"); return &y }()))
		(*bignum) = (*one) / (*smlnum)
		Dlabad(smlnum, bignum)
		(*anorm) = (*Dlantb(func() *byte {y := byte('M'); return &y }(), (uplo), (diag), (n), (kd), (ab), (ldab), (work)))
		//
		(*(rat)) = (*rmax) * (Min((*bignum)/MAX((*one), (*anorm)), (*one)/(*eps)))
	}
	//
	return
	//
	//     End of Dtbt06
	//
}
