package goblas

import 

// \brief \b Saxpy
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Saxpy(N,SA,SX,incx,SY,incy)
//
//       .. Scalar Arguments ..
//       REAL SA
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
//    Saxpy constant times a vector plus a vector.
//    uses unrolled loops for increments equal to one.
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
// \param[in] SA
// \verbatim
//          SA is REAL
//           On entry, SA specifies the scalar alpha.
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
// \param[in,out] SY
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
func Saxpy(n *int, sa *float64, sx *[]float64, incx *int, sy *[]float64, incy *int) {
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
	if (*sa) == 0.0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*        code for both increments equal to 1
		//*
		//*
		//*        clean-up loop
		//*
		(*m) = (MOD((*n), int(4)))
		if (*m) != 0 {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*sy)[(*i)-1] = (*sy)[(*i)-1] + (*sa)*(*sx)[(*i)-1]
			}
		}
		if (*n) < 4 {
			return
		}
		(*mp1) = (*m) + 1
		for (*i) = (*mp1); (*i) <= (*n); (*i) += 4 {
			(*sy)[(*i)-1] = (*sy)[(*i)-1] + (*sa)*(*sx)[(*i)-1]
			(*sy)[(*i)+0] = (*sy)[(*i)+0] + (*sa)*(*sx)[(*i)+0]
			(*sy)[(*i)+1] = (*sy)[(*i)+1] + (*sa)*(*sx)[(*i)+1]
			(*sy)[(*i)+2] = (*sy)[(*i)+2] + (*sa)*(*sx)[(*i)+2]
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
			(*sy)[(*iy)-1] = (*sy)[(*iy)-1] + (*sa)*(*sx)[(*ix)-1]
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	return
}
