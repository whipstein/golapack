package goblas

import 

// \brief \b Idamax
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION Idamax(N,DX,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Idamax finds the index of the first element having maximum absolute value.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] DX
// \verbatim
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of DX
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
// \date November 2017
//
// \ingroup aux_blas
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Idamax(n *int, dx *[]float64, incx *int) (idamaxReturn *int) {
	idamaxreturn := new(int)
	dmax := new(float64)
	i := new(int)
	ix := new(int)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	(*idamax) = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	(*idamax) = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		(*dmax) = (ABS(((*dx)[0])))
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if ABS(((*dx)[(*i)-1])) > (*dmax) {
				(*idamax) = (*i)
				(*dmax) = (ABS(((*dx)[(*i)-1])))
			}
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*ix) = 1
		(*dmax) = (ABS(((*dx)[0])))
		(*ix) = (*ix) + (*incx)
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if ABS(((*dx)[(*ix)-1])) > (*dmax) {
				(*idamax) = (*i)
				(*dmax) = (ABS(((*dx)[(*ix)-1])))
			}
			(*ix) = (*ix) + (*incx)
		}
	}
	return
}
