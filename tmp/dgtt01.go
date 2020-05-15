package goblas

import 

// Dgtt01 re_constructs a tridiagonal matrix A from its LU factorization
// and computes the residual
//    norm(L*U - A) / ( norm(a) * eps),
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
//       SUBROUTinE Dgtt01( n, DL, d, du, dlf, df, duF, du2, ipiv, work,
//                          ldwork, rwork, resid)
//
//       .. Scalar Arguments ..
//       inTEGER            ldwork, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       inTEGER            ipiv(*)
//       DOUBLE PRECISION   d(*), df(*), DL(*), dlf(*), du(*),
//      $                   du2(*), duF(*), rwork(*),
//      $                   work( ldwork, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dgtt01 re_constructs a tridiagonal matrix A from its LU factorization
// and computes the residual
//    norm(L*U - A) / ( norm(a) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is inTEGTER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] DL
// \verbatim
//          DL is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) sub-diagonal elements of A.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (n)
//          The diagonal elements of A.
// \endverbatim
//
// \param[in] du
// \verbatim
//          du is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) super-diagonal elements of A.
// \endverbatim
//
// \param[in] dlf
// \verbatim
//          dlf is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) multipliers that define the matrix L from the
//          LU factorization of A.
// \endverbatim
//
// \param[in] df
// \verbatim
//          df is DOUBLE PRECISION array, dimension (n)
//          The n diagonal elements of the upper triangular matrix U from
//          the LU factorization of A.
// \endverbatim
//
// \param[in] duF
// \verbatim
//          duF is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) elements of the first super-diagonal of U.
// \endverbatim
//
// \param[in] du2
// \verbatim
//          du2 is DOUBLE PRECISION array, dimension (N-2)
//          The (n-2) elements of the second super-diagonal of U.
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is inTEGER array, dimension (n)
//          The pivot indices; for 1 <= i <= n, row i of the matrix was
//          interchanged with row ipiv(i).  ipiv(i) will always be either
//          i or i+1; ipiv(i) = i indicates a row interchange was not
//          required.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (ldwork,N)
// \endverbatim
//
// \param[in] ldwork
// \verbatim
//          ldwork is inTEGER
//          The leading dimension of the array work.  ldwork >= max(1,N).
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
//          The scaled residual:  norm(L*U - A) / (norm(a) * eps)
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
func Dgtt01(n *int, dl *[]float64, d *[]float64, du *[]float64, dlf *[]float64, df *[]float64, duF *[]float64, du2 *[]float64, ipiv *[]int, work *[][]float64, ldwork *int, rwork *[]float64, resid *float64) {
	one := new(float64)
	zero := new(float64)
	i := new(int)
	ip := new(int)
	j := new(int)
	lastj := new(int)
	anorm := new(float64)
	eps := new(float64)
	li := new(float64)
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
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
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
	//     Quick return if possible
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y}()))
	//
	//     Copy the matrix U to work.
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(work))[(*i)-1][(*j)-1] = (*zero)
			//Label10:
		}
		//Label20:
	}
	for (*i) = 1; (*i) <= (*(n)); (*i)++ {
		if (*i) == 1 {
			(*(work))[(*i)-1][(*i)-1] = (*(df))[(*i)-1]
			if (*(n)) >= 2 {
				(*(work))[(*i)-1][(*i)+0] = (*(duF))[(*i)-1]
			}
			if (*(n)) >= 3 {
				(*(work))[(*i)-1][(*i)+1] = (*(du2))[(*i)-1]
			}
		} else if (*i) == (*(n)) {
			(*(work))[(*i)-1][(*i)-1] = (*(df))[(*i)-1]
		} else {
			(*(work))[(*i)-1][(*i)-1] = (*(df))[(*i)-1]
			(*(work))[(*i)-1][(*i)+0] = (*(duF))[(*i)-1]
			if (*i) < (*(n))-1 {
				(*(work))[(*i)-1][(*i)+1] = (*(du2))[(*i)-1]
			}
		}
		//Label30:
	}
	//
	//     Multiply on the left by L.
	//
	(*lastj) = (*(n))
	for (*i) = (*(n)) - 1; (*i) <= 1; (*i) += -1 {
		(*LI) = (*(dlf))[(*i)-1]
		Daxpy((*lastj)-(*i)+1, LI, &((*(work))[(*i)-1][(*i)-1]), (ldwork), &((*(work))[(*i)+0][(*i)-1]), (ldwork))
		(*iP) = (*(ipiv))[(*i)-1]
		if (*iP) == (*i) {
			(*lastj) = (Min((*i)+2, (*(n))))
		} else {
			Dswap((*lastj)-(*i)+1, &((*(work))[(*i)-1][(*i)-1]), (ldwork), &((*(work))[(*i)+0][(*i)-1]), (ldwork))
		}
		//Label40:
	}
	//
	//     Subtract the matrix A.
	//
	(*(work))[0][0] = (*(work))[0][0] - (*(d))[0]
	if (*(n)) > 1 {
		(*(work))[0][1] = (*(work))[0][1] - (*(du))[0]
		(*(work))[(*(n))-1][(*(n))-0] = (*(work))[(*(n))-1][(*(n))-0] - (*(dl))[(*(n))-0]
		(*(work))[(*(n))-1][(*(n))-1] = (*(work))[(*(n))-1][(*(n))-1] - (*(d))[(*(n))-1]
		for (*i) = 2; (*i) <= (*(n))-1; (*i)++ {
			(*(work))[(*i)-1][(*i)-0] = (*(work))[(*i)-1][(*i)-0] - (*(dl))[(*i)-0]
			(*(work))[(*i)-1][(*i)-1] = (*(work))[(*i)-1][(*i)-1] - (*(d))[(*i)-1]
			(*(work))[(*i)-1][(*i)+0] = (*(work))[(*i)-1][(*i)+0] - (*(du))[(*i)-1]
			//Label50:
		}
	}
	//
	//     Compute the 1-norm of the tridiagonal matrix A.
	//
	(*anorm) = (*Dlangt(func() *byte {y := byte('1'); return &y}(), (n), (dl), (d), (du)))
	//
	//     Compute the 1-norm of work, which is only guaranteed to be
	//     upper Hessenberg.
	//
	(*(resid)) = (*Dlanhs(func() *byte {y := byte('1'); return &y}(), (n), (work), (ldwork), (rwork)))
	//
	//     Compute norm(L*U - A) / (norm(a) * eps)
	//
	if (*anorm) <= (*zero) {
		if (*(resid)) != (*zero) {
			(*(resid)) = (*one) / (*eps)
		}
	} else {
		(*(resid)) = ((*(resid)) / (*anorm)) / (*eps)
	}
	//
	return
	//
	//     End of Dgtt01
	//
}
