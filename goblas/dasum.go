package goblas

import 

// \brief \b Dasum
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dasum(N,DX,incx)
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
//    Dasum takes the sum of the absolute values.
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
// \ingroup double_blas_level1
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
func Dasum(n *int, dx *[]float64, incx *int) (dasumReturn *float64) {
	dasumreturn := new(float64)
	dtemp := new(float64)
	i := new(int)
	m := new(int)
	mp1 := new(int)
	nincx := new(int)
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
	(*dasum) = 0.0
	(*dtemp) = 0.0
	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//*        code for increment equal to 1
		//*
		//*
		//*        clean-up loop
		//*
		(*m) = (MOD((*n), int(6)))
		if (*m) != 0 {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*dtemp) = (*dtemp) + ABS(((*dx)[(*i)-1]))
			}
			if (*n) < 6 {
				(*dasum) = (*dtemp)
				return
			}
		}
		(*mp1) = (*m) + 1
		for (*i) = (*mp1); (*i) <= (*n); (*i) += 6 {
			(*dtemp) = (*dtemp) + ABS(((*dx)[(*i)-1])) + ABS(((*dx)[(*i)+0])) + ABS(((*dx)[(*i)+1])) + ABS(((*dx)[(*i)+2])) + ABS(((*dx)[(*i)+3])) + ABS(((*dx)[(*i)+4]))
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*nincx) = (*n) * (*incx)
		for (*i) = 1; (*i) <= (*nincx); (*i) += (*incx) {
			(*dtemp) = (*dtemp) + ABS(((*dx)[(*i)-1]))
		}
	}
	(*dasum) = (*dtemp)
	return
}
