package goblas

import 

// \brief \b Isamax
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION Isamax(N,SX,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       REAL SX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Isamax finds the index of the first element having maximum absolute value.
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
// \param[in] SX
// \verbatim
//          SX is REAL array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of SX
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
func Isamax(n *int, sx *[]float64, incx *int) (isamaxReturn *int) {
	isamaxreturn := new(int)
	smax := new(float64)
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
	(*isamax) = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	(*isamax) = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		(*smax) = (ABS(((*sx)[0])))
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if ABS(((*sx)[(*i)-1])) > (*smax) {
				(*isamax) = (*i)
				(*smax) = (ABS(((*sx)[(*i)-1])))
			}
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*ix) = 1
		(*smax) = (ABS(((*sx)[0])))
		(*ix) = (*ix) + (*incx)
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if ABS(((*sx)[(*ix)-1])) > (*smax) {
				(*isamax) = (*i)
				(*smax) = (ABS(((*sx)[(*ix)-1])))
			}
			(*ix) = (*ix) + (*incx)
		}
	}
	return
}
