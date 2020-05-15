package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtpt01 computes the residual for a triangular matrix A times its
// inverse when A is stored in packed format:
//    resid = norm(A*ainv - I) / ( N * norm(a) * norm(ainv) * eps),
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
//       SUBROUTinE Dtpt01( uplo, diag, n, AP, ainvp, rcond, work, resid)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, uplo
//       intEGER            N
//       DOUBLE PRECISION   rcond, resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   ainvp(*), AP(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtpt01 computes the residual for a triangular matrix A times its
// inverse when A is stored in packed format:
//    resid = norm(A*ainv - I) / ( N * norm(a) * norm(ainv) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the matrix A is upper or lower triangular.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
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
// \param[in] AP
// \verbatim
//          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The original upper or lower triangular matrix a, packed
//          columnwise in a linear array.  The j-th column of A is stored
//          in the array AP as follows:
//          if uplo = 'U', AP((j-1)*j/2 + i) = a(i,j) for 1<=i<=j;
//          if uplo = 'L',
//             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = a(i,j) for j<=i<=n.
// \endverbatim
//
// \param[in,out] ainvp
// \verbatim
//          ainvp is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          On entry, the (triangular) inverse of the matrix a, packed
//          columnwise in a linear array as in AP.
//          On exit, the contents of ainvp are destroyed.
// \endverbatim
//
// \param[out] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The reciprocal condition number of a, computed as
//          1/(norm(a) * norm(ainv)).
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(A*ainv - I) / ( N * norm(a) * norm(ainv) * eps)
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
func Dtpt01(uplo *byte, diag *byte, n *int, ap *[]float64, ainvp *[]float64, rcond *float64, work *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	unitd := new(bool)
	j := new(int)
	jc := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	eps := new(float64)
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
		(*(rcond)) = (*one)
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if anorm = 0 or ainvnm = 0.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*anorm) = (*Dlantp(func() *byte {y := byte('1'); return &y }(), (uplo), (diag), (n), (ap), (work)))
	(*ainvnm) = (*Dlantp(func() *byte {y := byte('1'); return &y }(), (uplo), (diag), (n), (ainvp), (work)))
	if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
		(*(rcond)) = (*zero)
		(*(resid)) = (*one) / (*eps)
		return
	}
	(*(rcond)) = ((*one) / (*anorm)) / (*ainvnm)
	//
	//     Compute A * ainv, overwriting ainv.
	//
	(*unitd) = (*blas.Lsame((diag), func() *byte {y := byte('U'); return &y }()))
	if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
		(*jc) = 1
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			if *unitd {
				(*(ainvp))[(*jc)+(*j)-0] = (*one)
			}
			//
			//           Form the j-th column of A*ainv
			//
			Dtpmv(func() *[]byte {y := []byte("Upper"); return &y }(), func() *[]byte {y := []byte("No transpose"); return &y }(), (diag), j, (ap), &((*(ainvp))[(*jc)-(1)]), func() *int {y := 1; return &y }())
			//
			//           Subtract 1 from the diagonal
			//
			(*(ainvp))[(*jc)+(*j)-0] = (*(ainvp))[(*jc)+(*j)-0] - (*one)
			(*jc) = (*jc) + (*j)
			//Label10:
		}
	} else {
		(*jc) = 1
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			if *unitd {
				(*(ainvp))[(*jc)-(1)] = (*one)
			}
			//
			//           Form the j-th column of A*ainv
			//
			Dtpmv(func() *[]byte {y := []byte("Lower"); return &y }(), func() *[]byte {y := []byte("No transpose"); return &y }(), (diag), (*(n))-(*j)+1, &((*(ap))[(*jc)-(1)]), &((*(ainvp))[(*jc)-(1)]), func() *int {y := 1; return &y }())
			//
			//           Subtract 1 from the diagonal
			//
			(*(ainvp))[(*jc)-(1)] = (*(ainvp))[(*jc)-(1)] - (*one)
			(*jc) = (*jc) + (*(n)) - (*j) + 1
			//Label20:
		}
	}
	//
	//     Compute norm(A*ainv - I) / (n * norm(a) * norm(ainv) * eps)
	//
	(*(resid)) = (*Dlantp(func() *byte {y := byte('1'); return &y }(), (uplo), func() *[]byte {y := []byte("Non-unit"); return &y }(), (n), (ainvp), (work)))
	//
	(*(resid)) = (((*(resid)) * (*(rcond))) / DBLE((*(n)))) / (*eps)
	//
	return
	//
	//     End of Dtpt01
	//
}
