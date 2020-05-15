package goblas

import 

// \brief \b Sdsdot
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION Sdsdot(N,SB,SX,incx,SY,incy)
//
//       .. Scalar Arguments ..
//       REAL SB
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       REAL SX(//),SY(//)
//       ..
//
// \par Purpose:
//  =============
//
// \verbatim
//
//   Compute the inner product of two vectors with extended
//   precision accumulation.
//
//   Returns S.P. result with dot product accumulated in D.P.
//   Sdsdot = SB + sum for I = 0 to N-1 of SX(LX+I//incx)//SY(LY+I//incy),
//   where LX = 1 if incx .GE. 0, else LX = 1+(1-N)//incx, and LY is
//   defined in a similar way using incy.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is INTEGER
//          number of elements in input vector(s)
// \endverbatim
//
// \param[in] SB
// \verbatim
//          SB is REAL
//          single precision scalar to be added to inner product
// \endverbatim
//
// \param[in] SX
// \verbatim
//          SX is REAL array, dimension ( 1 + ( N - 1)//abs( incx))
//          single precision vector with N elements
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//          storage spacing between elements of SX
// \endverbatim
//
// \param[in] SY
// \verbatim
//          SY is REAL array, dimension ( 1 + ( N - 1)//abs( incx))
//          single precision vector with N elements
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//          storage spacing between elements of SY
// \endverbatim
//
//  Authors:
//  ========
//
// \author Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
// \author Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
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
//    REFERENCES
//
//    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
//    Krogh, Basic linear algebra subprograms for Fortran
//    usage, Algorithm No. 539, Transactions on Mathematical
//    Software 5, 3 (September 1979), pp. 308-323.
//
//    REVISION HISTORY  (YYMMDD)
//
//    791001  DATE WRITTEN
//    890531  Changed all specific intrinsics to generic.  (WRB)
//    890831  Modified array declarations.  (WRB)
//    890831  REVISION DATE from Version 3.2
//    891214  Prologue converted to Version 4.0 format.  (BAB)
//    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
//    920501  Reformatted the REFERENCES section.  (WRB)
//    070118  Reformat to LAPACK coding style
// \endverbatim
//
//  =====================================================================
func Sdsdot(n *int, sb *float64, sx *[]float64, incx *int, sy *[]float64, incy *int) (sdsdotReturn *float64) {
	sdsdotreturn := new(float64)
	dsdot := new(float64)
	i := new(int)
	kx := new(int)
	ky := new(int)
	ns := new(int)
	//*
	//*  -- Reference BLAS level1 routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     November 2017
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	(*dsdot) = (*sb)
	if (*n) <= 0 {
		(*sdsdot) = (*dsdot)
		return
	}
	if (*incx) == (*incy) && (*incx) > 0 {
		//*
		//*     Code for equal and positive increments.
		//*
		(*ns) = (*n) * (*incx)
		for (*i) = 1; (*i) <= (*ns); (*i) += (*incx) {
			(*dsdot) = (*dsdot) + DBLE(((*sx)[(*i)-1]))*DBLE(((*sy)[(*i)-1]))
		}
	} else {
		//*
		//*     Code for unequal or nonpositive increments.
		//*
		(*kx) = 1
		(*ky) = 1
		if (*incx) < 0 {
			(*kx) = 1 + (1-(*n))*(*incx)
		}
		if (*incy) < 0 {
			(*ky) = 1 + (1-(*n))*(*incy)
		}
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*dsdot) = (*dsdot) + DBLE(((*sx)[(*kx)-1]))*DBLE(((*sy)[(*ky)-1]))
			(*kx) = (*kx) + (*incx)
			(*ky) = (*ky) + (*incy)
		}
	}
	(*sdsdot) = (*dsdot)
	return
}
