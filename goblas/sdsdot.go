package goblas

// Sdsdot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       REAL FUNCTION SDSDOT(n,sb,sx,incx,sy,incy)
//
//       .. Scalar Arguments ..
//       REAL sb
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       REAL sx(*),sy(*)
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
//   SDSDOT = sb + sum for i = 0 to n-1 of sx(LX+i*incx)*sy(LY+i*incy),
//   where LX = 1 if incx .GE. 0, else LX = 1+(1-n)*incx, and LY is
//   defined in a similar way using incy.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is INTEGER
//          number of elements in input vector(s)
// \endverbatim
//
// \param[in] sb
// \verbatim
//          sb is REAL
//          single precision scalar to be added to inner product
// \endverbatim
//
// \param[in] sx
// \verbatim
//          sx is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
//          single precision vector with n elements
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//          storage spacing between elements of sx
// \endverbatim
//
// \param[in] sy
// \verbatim
//          sy is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
//          single precision vector with n elements
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//          storage spacing between elements of sy
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
func Sdsdot(major *byte, n *int, sb *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) float32 {
	var dsdot float32
	var i, kx, ky, ns int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	dsdot = *sb
	if *n <= 0 {
		return dsdot
	}
	if *incx == *incy && *incx > 0 {
		//
		//     Code for equal and positive increments.
		//
		ns = (*n) * (*incx)
		for i = 1; i <= ns; i += *incx {
			dsdot += (*sx)[i-1] * (*sy)[i-1]
		}
	} else {
		//
		//     Code for unequal or nonpositive increments.
		//
		kx = 1
		ky = 1
		if *incx < 0 {
			kx = 1 + (1-(*n))*(*incx)
		}
		if *incy < 0 {
			ky = 1 + (1-(*n))*(*incy)
		}
		for i = 1; i <= *n; i++ {
			dsdot += (*sx)[kx-1] * (*sy)[ky-1]
			kx += *incx
			ky += *incy
		}
	}
	return dsdot
}
