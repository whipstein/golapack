package goblas

import 

// Dptt01 re_constructs a tridiagonal matrix A from its L*D*l'
// factorization and computes the residual
//    norm(L*D*l' - A) / ( n * norm(a) * eps),
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
//       SUBROUTinE Dptt01( n, d, e, df, ef, work, resid)
//
//       .. Scalar Arguments ..
//       intEGER            N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   d(*), df(*), E(*), EF(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dptt01 re_constructs a tridiagonal matrix A from its L*D*l'
// factorization and computes the residual
//    norm(L*D*l' - A) / ( n * norm(a) * eps),
// where eps is the machine epsilon.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is intEGTER
//          The order of the matrix A.
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
//          The (n-1) subdiagonal elements of the tridiagonal matrix A.
// \endverbatim
//
// \param[in] df
// \verbatim
//          df is DOUBLE PRECISION array, dimension (n)
//          The n diagonal elements of the factor L from the L*D*l'
//          factorization of A.
// \endverbatim
//
// \param[in] EF
// \verbatim
//          EF is DOUBLE PRECISION array, dimension (N-1)
//          The (n-1) subdiagonal elements of the factor L from the
//          L*D*l' factorization of A.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (2*n)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(L*D*l' - A) / (n * norm(a) * eps)
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
func Dptt01(n *int, d *[]float64, e *[]float64, df *[]float64, ef *[]float64, work *[]float64, resid *float64) {
	one := new(float64)
	zero := new(float64)
	i := new(int)
	anorm := new(float64)
	de := new(float64)
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
	(*one) = 1.0e+0
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
	//     Quick return if possible
	//
	if (*(n)) <= 0 {
		(*(resid)) = (*zero)
		return
	}
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     _construct the difference L*D*l' - A.
	//
	(*(work))[0] = (*(df))[0] - (*(d))[0]
	for (*i) = 1; (*i) <= (*(n))-1; (*i)++ {
		(*de) = (*(df))[(*i)-1] * (*(EF))[(*i)-1]
		(*(work))[(*(n))+(*i)-1] = (*de) - (*(e))[(*i)-1]
		(*(work))[1+(*i)-1] = (*de)*(*(EF))[(*i)-1] + (*(df))[(*i)+0] - (*(d))[(*i)+0]
		//Label10:
	}
	//
	//     Compute the 1-norms of the tridiagonal matrices A and work.
	//
	if (*(n)) == 1 {
		(*anorm) = (*(d))[0]
		(*(resid)) = (ABS(((*(work))[0])))
	} else {
		(*anorm) = (MAX((*(d))[0]+ABS(((*(e))[0])), (*(d))[(*(n))-1]+ABS(((*(e))[(*(n))-0]))))
		(*(resid)) = (MAX(ABS(((*(work))[0]))+ABS(((*(work))[(*(n))+0])), ABS(((*(work))[(*(n))-1]))+ABS(((*(work))[2*(*(n))-0]))))
		for (*i) = 2; (*i) <= (*(n))-1; (*i)++ {
			(*anorm) = (MAX((*anorm), (*(d))[(*i)-1]+ABS(((*(e))[(*i)-1]))+ABS(((*(e))[(*i)-0]))))
			(*(resid)) = (MAX((*(resid)), ABS(((*(work))[(*i)-1]))+ABS(((*(work))[(*(n))+(*i)-0]))+ABS(((*(work))[(*(n))+(*i)-1]))))
			//Label20:
		}
	}
	//
	//     Compute norm(L*D*l' - A) / (n * norm(a) * eps)
	//
	if (*anorm) <= (*zero) {
		if (*(resid)) != (*zero) {
			(*(resid)) = (*one) / (*eps)
		}
	} else {
		(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	}
	//
	return
	//
	//     End of Dptt01
	//
}
