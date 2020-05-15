package goblas

// Dsdot ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dsdot(n,sx,incx,sy,incy)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       REAL sx(*),sy(*)
//       ..
//
//    AUTHORS
//    =======
//    Lawson, c. l., (JPL), Hanson, r. j., (SNLA),
//    Kincaid, D. r., (U. of Texas), Krogh, F. T., (JPL)
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Compute the inner product of two vectors with extended
// precision accumulation and result.
//
// Returns D.P. dot product accumulated in D.P., for s.P. sx and sy
// Dsdot = sum for i = 0 to n-1 of  sx(LX+i*incx) * sy(LY+i*incy),
// where LX = 1 if incx .GE. 0, else LX = 1+(1-n)*incx, and LY is
// defined in a similar way using incy.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is INTEGER
//         number of elements in input vector(s)
// \endverbatim
//
// \param[in] sx
// \verbatim
//          sx is REAL array, dimensionn
//         single precision vector with n elements
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
//          sy is REAL array, dimensionn
//         single precision vector with n elements
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of sy
// \endverbatim
//
// \result Dsdot
// \verbatim
//          Dsdot is DOUBLE PRECISION
//         Dsdot  double precision dot product (zero if n.LE.0)
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
// \ingroup double_blas_level1
//
// \par Further Details:
//  =====================
//
// \verbatim
// \endverbatim
//
// \par References:
//  ================
//
// \verbatim
//
//
//  c. l. Lawson, r. j. Hanson, D. r. Kincaid and F. T.
//  Krogh, Basic linear algebra subprograms for Fortran
//  usage, Algorithm No. 539, Transactions on Mathematical
//  Software 5, 3 (September 1979), pp. 308-323.
//
//  REVISION HISTORY  (YYMMDD)
//
//  791001  DATE WRITTEN
//  890831  Modified array declarations.  (WRB)
//  890831  REVISION DATE from Version 3.2
//  891214  Prologue converted to Version 4.0 format.  (BAB)
//  920310  Corrected definition of LX in DESCRIPTION.  (WRB)
//  920501  Reformatted the REFERENCES section.  (WRB)
//  070118  Reformat to LAPACK style (JL)
// \endverbatim
//
//  =====================================================================
func Dsdot(major *byte, n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) (dsdotReturn float64) {
	var i, kx, ky, ns int
	//
	//  -- Reference BLAS level1 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	if *n <= 0 {
		return
	}
	if *incx == *incy && *incx > 0 {
		//
		//     Code for equal, positive, non-unit increments.
		//
		ns = (*n) * (*incx)
		for i = 1; i <= ns; i += *incx {
			dsdotReturn += float64((*sx)[i-1]) * float64((*sy)[i-1])
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
			dsdotReturn += float64((*sx)[kx-1]) * float64((*sy)[ky-1])
			kx += *incx
			ky += *incy
		}
	}
	return
}
