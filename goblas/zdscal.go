package goblas

// \brief \b Zdscal
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zdscal(N,DA,ZX,incx)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION DA
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 ZX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Zdscal scales a vector by a constant.
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
// \param[in] DA
// \verbatim
//          DA is DOUBLE PRECISION
//           On entry, DA specifies the scalar alpha.
// \endverbatim
//
// \param[in,out] ZX
// \verbatim
//          ZX is COMPLEX//16 array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of ZX
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
// \ingroup complex16_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Zdscal(n *int, da *float64, zx *[]complex128, incx *int) {
	i := new(int)
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
	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*zx)[(*i)-1] = DCMPLX(da, func() *float64 { y := 0.0; return &y }()) * (*zx)[(*i)-1]
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*nincx) = (*n) * (*incx)
		for (*i) = 1; (*i) <= (*nincx); (*i) += (*incx) {
			(*zx)[(*i)-1] = DCMPLX(da, func() *float64 { y := 0.0; return &y }()) * (*zx)[(*i)-1]
		}
	}
	return
}
