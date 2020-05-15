package goblas

// Dlaptm multiplies an N by nrhs matrix X by a symmetric tridiagonal
// matrix A and stores the result in a matrix B.  The operation has the
// form
//
//    b := alpha * A * X + beta * B
//
// where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlaptm( n, nrhs, alpha, d, e, x, ldx, beta, b, ldb)
//
//       .. Scalar Arguments ..
//       intEGER            ldb, ldx, n, nrhs
//       DOUBLE PRECISION   alpha, beta
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   B( ldb, *), d(*), E(*), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlaptm multiplies an N by nrhs matrix X by a symmetric tridiagonal
// matrix A and stores the result in a matrix B.  The operation has the
// form
//
//    b := alpha * A * X + beta * B
//
// where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The order of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of columns
//          of the matrices X and B.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is DOUBLE PRECISION
//          The scalar alpha.  alpha must be 1. or -1.; otherwise,
//          it is assumed to be 0.
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (n)
//          The n diagonal elements of the tridiagonal matrix A.
// \endverbatim
//
// \param[in] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) subdiagonal or superdiagonal elements of A.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          The N by nrhs matrix X.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is intEGER
//          The leading dimension of the array X.  ldx >= max(n,1).
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is DOUBLE PRECISION
//          The scalar beta.  beta must be 0., 1., or -1.; otherwise,
//          it is assumed to be 1.
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          On entry, the N by nrhs matrix B.
//          On exit, B is overwritten by the matrix expression
//          b := alpha * A * X + beta * B.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B.  ldb >= max(n,1).
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
func Dlaptm(n *int, nrhs *int, alpha *float64, d *[]float64, e *[]float64, x *[][]float64, ldx *int, beta *float64, b *[][]float64, ldb *int) {
	one := new(float64)
	zero := new(float64)
	i := new(int)
	j := new(int)
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
	//     .. Executable Statements ..
	//
	if (*(n)) == 0 {
		return
	}
	//
	//     Multiply B by beta if beta.NE.1.
	//
	if (*(beta)) == (*zero) {
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*(b))[(*i)-(1)][(*j)-(1)] = (*zero)
				//Label10:
			}
			//Label20:
		}
	} else if (*(beta)) == -(*one) {
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*(b))[(*i)-(1)][(*j)-(1)] = -(*(b))[(*i)-(1)][(*j)-(1)]
				//Label30:
			}
			//Label40:
		}
	}
	//
	if (*(alpha)) == (*one) {
		//
		//        Compute b := B + A*X
		//
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			if (*(n)) == 1 {
				(*(b))[0][(*j)-(1)] = (*(b))[0][(*j)-(1)] + (*(d))[0]*(*(x))[0][(*j)-(1)]
			} else {
				(*(b))[0][(*j)-(1)] = (*(b))[0][(*j)-(1)] + (*(d))[0]*(*(x))[0][(*j)-(1)] + (*(e))[0]*(*(x))[1][(*j)-(1)]
				(*(b))[(*(n))-(1)][(*j)-(1)] = (*(b))[(*(n))-(1)][(*j)-(1)] + (*(e))[(*(n))-0]*(*(x))[(*(n))-0][(*j)-(1)] + (*(d))[(*(n))-(1)]*(*(x))[(*(n))-(1)][(*j)-(1)]
				for (*i) = 2; (*i) <= (*(n))-1; (*i)++ {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*(b))[(*i)-(1)][(*j)-(1)] + (*(e))[(*i)-0]*(*(x))[(*i)-0][(*j)-(1)] + (*(d))[(*i)-(1)]*(*(x))[(*i)-(1)][(*j)-(1)] + (*(e))[(*i)-(1)]*(*(x))[(*i)+0][(*j)-(1)]
					//Label50:
				}
			}
			//Label60:
		}
	} else if (*(alpha)) == -(*one) {
		//
		//        Compute b := B - A*X
		//
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			if (*(n)) == 1 {
				(*(b))[0][(*j)-(1)] = (*(b))[0][(*j)-(1)] - (*(d))[0]*(*(x))[0][(*j)-(1)]
			} else {
				(*(b))[0][(*j)-(1)] = (*(b))[0][(*j)-(1)] - (*(d))[0]*(*(x))[0][(*j)-(1)] - (*(e))[0]*(*(x))[1][(*j)-(1)]
				(*(b))[(*(n))-(1)][(*j)-(1)] = (*(b))[(*(n))-(1)][(*j)-(1)] - (*(e))[(*(n))-0]*(*(x))[(*(n))-0][(*j)-(1)] - (*(d))[(*(n))-(1)]*(*(x))[(*(n))-(1)][(*j)-(1)]
				for (*i) = 2; (*i) <= (*(n))-1; (*i)++ {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*(b))[(*i)-(1)][(*j)-(1)] - (*(e))[(*i)-0]*(*(x))[(*i)-0][(*j)-(1)] - (*(d))[(*i)-(1)]*(*(x))[(*i)-(1)][(*j)-(1)] - (*(e))[(*i)-(1)]*(*(x))[(*i)+0][(*j)-(1)]
					//Label70:
				}
			}
			//Label80:
		}
	}
	return
	//
	//     End of Dlaptm
	//
}
