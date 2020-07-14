package golapack

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
//       DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       REAL SX(*),SY(*)
//       ..
//
//    AUTHORS
//    =======
//    Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
//    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
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
// Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
// DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
// where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
// defined in a similar way using INCY.
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
//          SX is REAL array, dimension(N)
//         single precision vector with N elements
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//          storage spacing between elements of SX
// \endverbatim
//
// \param[in] SY
// \verbatim
//          SY is REAL array, dimension(N)
//         single precision vector with N elements
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of SY
// \endverbatim
//
// \result DSDOT
// \verbatim
//          DSDOT is DOUBLE PRECISION
//         DSDOT  double precision dot product (zero if N.LE.0)
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
//  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
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
func Dsdot(n *int, sx *[]float32, sxoff, incx *int, sy *[]float32, syoff, incy *int) (dsdotReturn float64) {
	var i, kx, ky, ns int

	dsdotReturn = 0.0
	if (*n) <= 0 {
		return
	}
	if (*incx) == (*incy) && (*incx) > 0 {
		//
		//     Code for equal, positive, non-unit increments.
		//
		ns = (*n) * (*incx)
		for i = 1; i <= ns; i += (*incx) {
			dsdotReturn = dsdotReturn + float64((*sx)[i-1+(*sxoff)])*float64((*sy)[i-1+(*syoff)])
		}
	} else {
		//
		//     Code for unequal or nonpositive increments.
		//
		kx = 1
		ky = 1
		if (*incx) < 0 {
			kx = 1 + (1-(*n))*(*incx)
		}
		if (*incy) < 0 {
			ky = 1 + (1-(*n))*(*incy)
		}
		for i = 1; i <= (*n); i++ {
			dsdotReturn = dsdotReturn + float64((*sx)[kx-1+(*sxoff)])*float64((*sy)[ky-1+(*syoff)])
			kx = kx + (*incx)
			ky = ky + (*incy)
		}
	}
	return
}
