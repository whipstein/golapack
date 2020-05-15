package goblas

// \brief \b Zdrot
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zdrot( N, CX, incx, CY, incy, C, S)
//
//       .. Scalar Arguments ..
//       INTEGER            incx, incy, N
//       DOUBLE PRECISION   C, S
//       ..
//       .. Array Arguments ..
//       COMPLEX//16         CX( //), CY( //)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Applies a plane rotation, where the cos and sin (c and s) are real
// and the vectors cx and cy are complex.
// jack dongarra, linpack, 3/11/78.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the vectors cx and cy.
//           N must be at least zero.
// \endverbatim
//
// \param[in,out] CX
// \verbatim
//          CX is COMPLEX//16 array, dimension at least
//           ( 1 + ( N - 1)//abs( incx)).
//           Before entry, the incremented array CX must contain the n
//           element vector cx. On exit, CX is overwritten by the updated
//           vector cx.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           CX. incx must not be zero.
// \endverbatim
//
// \param[in,out] CY
// \verbatim
//          CY is COMPLEX//16 array, dimension at least
//           ( 1 + ( N - 1)//abs( incy)).
//           Before entry, the incremented array CY must contain the n
//           element vector cy. On exit, CY is overwritten by the updated
//           vector cy.
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//           On entry, incy specifies the increment for the elements of
//           CY. incy must not be zero.
// \endverbatim
//
// \param[in] C
// \verbatim
//          C is DOUBLE PRECISION
//           On entry, C specifies the cosine, cos.
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is DOUBLE PRECISION
//           On entry, S specifies the sine, sin.
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
// \ingroup complex16_blas_level1
//
//  =====================================================================
func Zdrot(n *int, cx *[]complex128, incx *int, cy *[]complex128, incy *int, c *float64, s *float64) {
	i := new(int)
	ix := new(int)
	iy := new(int)
	ctemp := new(complex128)
	//*
	//*  -- Reference BLAS level1 routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     December 2016
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//* =====================================================================
	//*
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Executable Statements ..
	//*
	if (*n) <= 0 {
		return
	}
	if (*incx) == 1 && (*incy) == 1 {
		//*
		//*        code for both increments equal to 1
		//*
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ctemp) = (*c)*(*cx)[(*i)-1] + (*s)*(*cy)[(*i)-1]
			(*cy)[(*i)-1] = (*c)*(*cy)[(*i)-1] - (*s)*(*cx)[(*i)-1]
			(*cx)[(*i)-1] = (*ctemp)
		}
	} else {
		//*
		//*        code for unequal increments or equal increments not equal
		//*          to 1
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
			(*ctemp) = (*c)*(*cx)[(*ix)-1] + (*s)*(*cy)[(*iy)-1]
			(*cy)[(*iy)-1] = (*c)*(*cy)[(*iy)-1] - (*s)*(*cx)[(*ix)-1]
			(*cx)[(*ix)-1] = (*ctemp)
			(*ix) = (*ix) + (*incx)
			(*iy) = (*iy) + (*incy)
		}
	}
	return
}
