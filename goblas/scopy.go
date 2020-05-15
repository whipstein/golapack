package goblas

import 
// \brief \b Scopy
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Scopy(N,SX,incx,SY,incy)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       REAL SX(//),SY(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Scopy copies a vector, x, to a vector, y.
//    uses unrolled loops for increments equal to 1.
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
// \param[out] SY
// \verbatim
//          SY is REAL array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of SY
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
// \ingroup single_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Scopy(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	i := new(int)
	ix := new(int)
	iy := new(int)
	m := new(int)
	mp1 := new(int)
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
	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*        code for both increments equal to 1
		//*
		//*
		//*        clean-up loop
		//*
		(*m) = (MOD((*n), int(7)))
		if (*m) != 0 {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*sy)[(*i)-1] = (*sx)[(*i)-1]
			}
			if (*n) < 7 {
				return
			}
		}
		(*mp1) = (*m) + 1
		for (*i) = (*mp1); (*i) <= (*n); (*i) += 7 {
			(*sy)[(*i)-1] = (*sx)[(*i)-1]
			(*sy)[(*i)+0] = (*sx)[(*i)+0]
			(*sy)[(*i)+1] = (*sx)[(*i)+1]
			(*sy)[(*i)+2] = (*sx)[(*i)+2]
			(*sy)[(*i)+3] = (*sx)[(*i)+3]
			(*sy)[(*i)+4] = (*sx)[(*i)+4]
			(*sy)[(*i)+5] = (*sx)[(*i)+5]
		}
	} else {
		//*
		//*        code for unequal increments or equal increments
		//*          not equal to 1
		//*
		(*ix) = 1
		(*iy) = 1
		if (*incx) < 0 {
			(*ix) = (-(*n)+1)*(*incx) + 1
		}
		if (*incy) < 0 {
			(*iy) = (-(*n)+1)*(*incy) + 1
		}
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*sy)[(*iy)-1] = (*sx)[(*ix)-1]
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	return
}
