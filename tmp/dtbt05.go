package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dtbt05 tests the error bounds from iterative refinement for the
// computed solution to a system of equations A*X = b, where A is a
// triangular band matrix.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dtbt05( uplo, trans, diag, n, kd, nrhs, AB, ldab, b,
//                          ldb, x, ldx, xact, ldxact, ferr, berr, reslts)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       intEGER            kd, ldab, ldb, ldx, ldxact, n, nrhs
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   AB( ldab, *), B( ldb, *), berr(*),
//      $                   ferr(*), reslts(*), X( ldx, *),
//      $                   xact( ldxact, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dtbt05 tests the error bounds from iterative refinement for the
// computed solution to a system of equations A*X = b, where A is a
// triangular band matrix.
//
// reslts1 = test of the error bound
//           = norm(X - xact) / ( norm(x) * ferr)
//
// A large value is returned if this ratio is not less than one.
//
// reslts(2) = residual from the iterative refinement routine
//           = the maximum of berr / ( nz*eps + (*)), where
//             (*) = nz*unfl / (min_i (abs(a)*abs(x) +abs(b))_i)
//             and nz = max. number of nonzeros in any row of a, plus 1
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
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the form of the system of equations.
//          = 'N':  A * X = B  (No transpose)
//          = 'T':  A'* X = B  (Transpose)
//          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
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
//          The number of rows of the matrices x, b, and xact, and the
//          order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kd
// \verbatim
//          kd is intEGER
//          The number of super-diagonals of the matrix A if uplo = 'U',
//          or the number of sub-diagonals if uplo = 'L'.  kd >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of the matrices x, b, and xact.
//          nrhs >= 0.
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
//          If diag = 'U', the diagonal elements of A are not referenced
//          and are assumed to be 1.
// \endverbatim
//
// \param[in] ldab
// \verbatim
//          ldab is intEGER
//          The leading dimension of the array AB.  ldab >= kd+1.
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
//          ldb is intEGER
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
//          ldx is intEGER
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
//          ldxact is intEGER
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
//          reslts(2) = berr / ( nz*eps + (*))
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
func Dtbt05(uplo *byte, trans *byte, diag *byte, n *int, kd *int, nrhs *int, ab *[][]float64, ldab *int, b *[][]float64, ldb *int, x *[][]float64, ldx *int, xact *[][]float64, ldxact *int, ferr *[]float64, berr *[]float64, reslts *[]float64) {
	zero := new(float64)
	one := new(float64)
	notran := new(bool)
	unit := new(bool)
	upper := new(bool)
	i := new(int)
	ifu := new(int)
	imax := new(int)
	j := new(int)
	k := new(int)
	nz := new(int)
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
	(*upper) = (*blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()))
	(*notran) = (*blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()))
	(*unit) = (*blas.Lsame((diag), func() *byte {y := byte('U'); return &y }()))
	(*nz) = Min((*(kd)), (*(n))-1) + 1
	//
	//     Test 1:  Compute the maximum of
	//        norm(X - xact) / ( norm(x) * ferr)
	//     over all the vectors X and xact using the infinity-norm.
	//
	(*errbnd) = (*zero)
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
	(*(reslts))[0] = (*errbnd)
	//
	//     Test 2:  Compute the maximum of berr / ( nz*eps + (*)), where
	//     (*) = nz*unfl / (min_i (abs(a)*abs(x) +abs(b))_i)
	//
	(*ifu) = 0
	if *unit {
		(*ifu) = 1
	}
	for (*k) = 1; (*k) <= (*(nrhs)); (*k)++ {
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*tmp) = (ABS(((*(b))[(*i)-1][(*k)-1])))
			if *upper {
				if !(*notran) {
					for (*j) = MAX((*i)-(*(kd)), 1); (*j) <= (*i)-(*ifu); (*j)++ {
						(*tmp) = (*tmp) + ABS(((*(ab))[(*(kd))+1-(*i)+(*j)-1][(*i)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
						//Label40:
					}
					if *unit {
						(*tmp) = (*tmp) + ABS(((*(x))[(*i)-1][(*k)-1]))
					}
				} else {
					if *unit {
						(*tmp) = (*tmp) + ABS(((*(x))[(*i)-1][(*k)-1]))
					}
					for (*j) = (*i) + (*ifu); (*j) <= (Min((*i)+(*(kd)), (*(n)))); (*j)++ {
						(*tmp) = (*tmp) + ABS(((*(ab))[(*(kd))+1+(*i)-(*j)-1][(*j)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
						//Label50:
					}
				}
			} else {
				if *notran {
					for (*j) = MAX((*i)-(*(kd)), 1); (*j) <= (*i)-(*ifu); (*j)++ {
						(*tmp) = (*tmp) + ABS(((*(ab))[1+(*i)-(*j)-1][(*j)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
						//Label60:
					}
					if *unit {
						(*tmp) = (*tmp) + ABS(((*(x))[(*i)-1][(*k)-1]))
					}
				} else {
					if *unit {
						(*tmp) = (*tmp) + ABS(((*(x))[(*i)-1][(*k)-1]))
					}
					for (*j) = (*i) + (*ifu); (*j) <= (Min((*i)+(*(kd)), (*(n)))); (*j)++ {
						(*tmp) = (*tmp) + ABS(((*(ab))[1+(*j)-(*i)-1][(*i)-1]))*ABS(((*(x))[(*j)-1][(*k)-1]))
						//Label70:
					}
				}
			}
			if (*i) == 1 {
				(*axbi) = (*tmp)
			} else {
				(*axbi) = (Min((*axbi), (*tmp)))
			}
			//Label80:
		}
		(*tmp) = (*(berr))[(*k)-1] / ((*nz)*(*eps) + (*nz)*(*unfl)/MAX((*axbi), (*nz)*(*unfl)))
		if (*k) == 1 {
			(*(reslts))[1] = (*tmp)
		} else {
			(*(reslts))[1] = (MAX(((*(reslts))[1]), (*tmp)))
		}
		//Label90:
	}
	//
	return
	//
	//     End of Dtbt05
	//
}
