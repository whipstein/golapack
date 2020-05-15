package goblas

// \brief \b Srot
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Srot(N,SX,incx,SY,incy,C,S)
//
//       .. Scalar Arguments ..
//       REAL C,S
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
//    applies a plane rotation.
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
// \param[in,out] SX
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
// \param[in] C
// \verbatim
//          C is REAL
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is REAL
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
func Srot(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, c *float64, s *float64) {
	stemp := new(float64)
	i := new(int)
	ix := new(int)
	iy := new(int)
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
	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*       code for both increments equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*stemp) = (*c)*(*sx)[(*i)-1] + (*s)*(*sy)[(*i)-1]
			(*sy)[(*i)-1] = (*c)*(*sy)[(*i)-1] - (*s)*(*sx)[(*i)-1]
			(*sx)[(*i)-1] = (*stemp)
		}
	} else {
		//*
		//*       code for unequal increments or equal increments not equal
		//*         to 1
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
			(*stemp) = (*c)*(*sx)[(*ix)-1] + (*s)*(*sy)[(*iy)-1]
			(*sy)[(*iy)-1] = (*c)*(*sy)[(*iy)-1] - (*s)*(*sx)[(*ix)-1]
			(*sx)[(*ix)-1] = (*stemp)
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	return
}
