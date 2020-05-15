package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dget07 tests the error bounds from iterative refinement for the
// computed solution to a system of equations op(a)*X = b, where A is a
// general n by n matrix and op(a) = A or A**T, depending on trans.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dget07( trans, n, nrhs, a, lda, b, ldb, x, ldx, xact,
//                          ldxact, ferr, chkferr, berr, reslts)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       LOGICAL            chkferr
//       inTEGER            lda, ldb, ldx, ldxact, n, nrhs
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), berr(*), ferr(*),
//      $                   reslts(*), X( ldx, *), XACt( ldxact, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dget07 tests the error bounds from iterative refinement for the
// computed solution to a system of equations op(a)*X = b, where A is a
// general n by n matrix and op(a) = A or A**T, depending on trans.
//
// reslts1 = test of the error bound
//           = norm(X - xact) / ( norm(x) * ferr)
//
// A large value is returned if this ratio is not less than one.
//
// reslts(2) = residual from the iterative refinement routine
//           = the maximum of berr / ( (n+1)*eps + (*)), where
//             (*) = (n+1)*unfl / (min_i (abs(op(a))*abs(x) +abs(b))_i)
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the form of the system of equations.
//          = 'N':  A * X = B     (No transpose)
//          = 'T':  A**T * X = B  (Transpose)
//          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of rows of the matrices X and xact.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The number of columns of the matrices X and xact.  nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original n by n matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.  lda >= max(1,N).
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          The right hand side vectors for the system of linear
//          equations.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B.  ldb >= max(1,N).
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The computed solution vectors.  Each vector is stored as a
//          column of the matrix X.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is inTEGER
//          The leading dimension of the array X.  ldx >= max(1,N).
// \endverbatim
//
// \param[in] xact
// \verbatim
//          xact is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The exact solution vectors.  Each vector is stored as a
//          column of the matrix xact.
// \endverbatim
//
// \param[in] ldxact
// \verbatim
//          ldxact is inTEGER
//          The leading dimension of the array xact.  ldxact >= max(1,N).
// \endverbatim
//
// \param[in] ferr
// \verbatim
//          ferr is DOUBLE PRECISION array, dimension (nrhs)
//          The estimated forward error bounds for each solution vector
//          X.  If xtRUE is the true solution, ferr bounds the magnitude
//          of the largest entry in (X - xtRUE) divided by the magnitude
//          of the largest entry in X.
// \endverbatim
//
// \param[in] chkferr
// \verbatim
//          chkferr is LOGICAL
//          Set to .TRUE. to check ferr, .FALSE. not to check ferr.
//          When the test system is ill-conditioned, the "true"
//          solution in xact may be incorrect.
// \endverbatim
//
// \param[in] berr
// \verbatim
//          berr is DOUBLE PRECISION array, dimension (nrhs)
//          The componentwise relative backward error of each solution
//          vector (i.e., the smallest relative change in any entry of A
//          or B that makes X an exact solution).
// \endverbatim
//
// \param[out] reslts
// \verbatim
//          reslts is DOUBLE PRECISION array, dimension (2)
//          The maximum over the nrhs solution vectors of the ratios:
//          reslts1 = norm(X - xact) / ( norm(x) * ferr)
//          reslts(2) = berr / ( (n+1)*eps + (*))
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
func Dget07(trans *byte, n *int, nrhs *int, a *[][]float64, lda *int, b *[][]float64, ldb *int, x *[][]float64, ldx *int, xact *[][]float64, ldxact *int, ferr *[]float64, chkferr *bool, berr *[]float64, reslts *[]float64) {
	zero := new(float64)
	one := new(float64)
	notran := new(bool)
	i := new(int)
	imax := new(int)
	j := new(int)
	k := new(int)
	axbi := new(float64)
	diff := new(float64)
	eps := new(float64)
	errbnd := new(float64)
	ovfl := new(float64)
	tmp := new(float64)
	unfl := new(float64)
	xnorm := new(float64)
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
	//     .. Executable Statements ..
	//
	//     Quick exit if N = 0 or nrhs = 0.
	//
	if (*(n)) <= 0 || (*(nrhs)) <= 0 {
		(*(reslts))[0] = (*zero)
		(*(reslts))[1] = (*zero)
		return
	}
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*unfl) = (*Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()))
	(*ovfl) = (*one) / (*unfl)
	(*notran) = (*blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()))
	//
	//     Test 1:  Compute the maximum of
	//        norm(X - xact) / ( norm(x) * ferr)
	//     over all the vectors X and xact using the infinity-norm.
	//
	(*errbnd) = (*zero)
	if *(chkferr) {
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			(*imax) = (*Idamax((n), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }()))
			(*xnorm) = (MAX(ABS(((*(x))[(*imax)-1][(*j)-1])), (*unfl)))
			(*diff) = (*zero)
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*diff) = (MAX((*diff), ABS((*(x))[(*i)-1][(*j)-1]-(*(xact))[(*i)-1][(*j)-1])))
				//Label10:
			}
			//
			if (*xnorm) > (*one) {
				goto Label20
			} else if (*diff) <= (*ovfl)*(*xnorm) {
				goto Label20
			} else {
				(*errbnd) = (*one) / (*eps)
				goto Label30
			}
			//
		Label20:
			;
			if (*diff)/(*xnorm) <= (*(ferr))[(*j)-1] {
				(*errbnd) = (MAX((*errbnd), ((*diff)/(*xnorm))/(*(ferr))[(*j)-1]))
			} else {
				(*errbnd) = (*one) / (*eps)
			}
		Label30:
		}
	}
	(*(reslts))[0] = (*errbnd)
	//
	//     Test 2:  Compute the maximum of berr / ( (n+1)*eps + (*)), where
	//     (*) = (n+1)*unfl / (min_i (abs(op(a))*abs(x) +abs(b))_i)
	//
	for (*k) = 1; (*k) <= (*(nrhs)); (*k)++ {
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*tmp) = (ABS(((*(b))[(*i)-1][(*k)-1])))
			if *notran {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					(*tmp) = (*tmp) + ABS(((*(a))[(*i)-1][(*j)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
					//Label40:
				}
			} else {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					(*tmp) = (*tmp) + ABS(((*(a))[(*j)-1][(*i)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
					//Label50:
				}
			}
			if (*i) == 1 {
				(*axbi) = (*tmp)
			} else {
				(*axbi) = (Min((*axbi), (*tmp)))
			}
			//Label60:
		}
		(*tmp) = (*(berr))[(*k)-1] / (((*(n))+1)*(*eps) + ((*(n))+1)*(*unfl)/MAX((*axbi), ((*(n))+1)*(*unfl)))
		if (*k) == 1 {
			(*(reslts))[1] = (*tmp)
		} else {
			(*(reslts))[1] = (MAX(((*(reslts))[1]), (*tmp)))
		}
		//Label70:
	}
	//
	return
	//
	//     End of Dget07
	//
}
