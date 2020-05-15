package goblas

import 

// Dget04 computes the difference between a computed solution and the
// true solution to a system of linear equations.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dget04( n, nrhs, x, ldx, xact, ldxact, rcond, resid)
//
//       .. Scalar Arguments ..
//       intEGER            ldx, ldxact, n, nrhs
//       DOUBLE PRECISION   rcond, resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   X( ldx, *), xact( ldxact, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dget04 computes the difference between a computed solution and the
// true solution to a system of linear equations.
//
// resid =  ( norm(X-xact) * rcond) / ( norm(xact) * eps),
// where rcond is the reciprocal of the condition number and eps is the
// machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of rows of the matrices X and xact.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of the matrices X and xact.  nrhs >= 0.
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
//          xact is DOUBLE PRECISION array, dimension( ldx, nrhs)
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
// \param[in] rcond
// \verbatim
//          rcond is DOUBLE PRECISION
//          The reciprocal of the condition number of the coefficient
//          matrix in the system of equations.
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          The maximum over the nrhs solution vectors of
//          ( norm(X-xact) * rcond) / ( norm(xact) * eps)
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
func Dget04(n *int, nrhs *int, x *[][]float64, ldx *int, xact *[][]float64, ldxact *int, rcond *float64, resid *float64) {
	zero := new(float64)
	i := new(int)
	ix := new(int)
	j := new(int)
	diffnm := new(float64)
	eps := new(float64)
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
		(*(resid)) = (*zero)
		return
	}
	//
	//     Exit with resid = 1/eps if rcond is invalid.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	if (*(rcond)) < (*zero) {
		(*(resid)) = 1.0 / (*eps)
		return
	}
	//
	//     Compute the maximum of
	//        norm(X - xact) / ( norm(xact) * eps)
	//     over all the vectors X and xact .
	//
	(*(resid)) = (*zero)
	for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
		(*ix) = (*Idamax((n), &((*(xact))[0][(*j)-(1)]), func() *int {y := 1; return &y }()))
		(*xnorm) = (ABS(((*(xact))[(*ix)-(1)][(*j)-(1)])))
		(*diffnm) = (*zero)
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*diffnm) = (MAX((*diffnm), ABS((*(x))[(*i)-(1)][(*j)-(1)]-(*(xact))[(*i)-(1)][(*j)-(1)])))
			//Label10:
		}
		if (*xnorm) <= (*zero) {
			if (*diffnm) > (*zero) {
				(*(resid)) = 1.0 / (*eps)
			}
		} else {
			(*(resid)) = (MAX((*(resid)), ((*diffnm)/(*xnorm))*(*(rcond))))
		}
		//Label20:
	}
	if (*(resid))*(*eps) < 1.0 {
		(*(resid)) = (*(resid)) / (*eps)
	}
	//
	return
	//
	//     End of Dget04
	//
}
