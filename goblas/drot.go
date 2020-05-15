package goblas

// \brief \b Drot
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drot(N,DX,incx,DY,incy,C,S)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION C,S
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DX(//),DY(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Drot applies a plane rotation.
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
// \param[in,out] DX
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
// \param[in,out] DY
// \verbatim
//          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of DY
// \endverbatim
//
// \param[in] C
// \verbatim
//          C is DOUBLE PRECISION
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is DOUBLE PRECISION
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
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Drot(n *int, dx *[]float64, incx *int, dy *[]float64, incy *int, c *float64, s *float64) {
	dtemp := new(float64)
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
	if (*incx) == 1 .(*and).(*incy) == 1 {
		//*
		//*       code for both increments equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*dtemp) = (*c)*(*dx)[(*i)-1] + (*s)*(*dy)[(*i)-1]
			(*dy)[(*i)-1] = (*c)*(*dy)[(*i)-1] - (*s)*(*dx)[(*i)-1]
			(*dx)[(*i)-1] = (*dtemp)
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
			(*dtemp) = (*c)*(*dx)[(*ix)-1] + (*s)*(*dy)[(*iy)-1]
			(*dy)[(*iy)-1] = (*c)*(*dy)[(*iy)-1] - (*s)*(*dx)[(*ix)-1]
			(*dx)[(*ix)-1] = (*dtemp)
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	return
}
