package goblas

// Srotm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SROTM(n,sx,incx,sy,incy,sparam)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       REAL sparam(5),sx(*),sy(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY n MATRIX
//
//    (sx**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF sx ARE IN
//    (sx**T)
//
//    sx(LX+i*incx), i = 0 TO n-1, WHERE LX = 1 IF incx .GE. 0, ELSE
//    LX = (-incx)*n, AND SIMILARLY FOR sy USING USING LY AND incy.
//    WITH sparam1=sflag, H HAS ONE OF THE FOLLOWING FORMS..
//
//    sflag=-1.E0     sflag=0.E0        sflag=1.E0     sflag=-2.E0
//
//      (sh11  sh12)    (1.E0  sh12)    (sh11  1.E0)    (1.E0  0.E0)
//    H=(          )    (          )    (          )    (          )
//      (sh21  sh22),   (sh21  1.E0),   (-1.E0 sh22),   (0.E0  1.E0).
//    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN sparam.
//
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
// \param[in,out] sx
// \verbatim
//          sx is REAL array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of sx
// \endverbatim
//
// \param[in,out] sy
// \verbatim
//          sy is REAL array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of sy
// \endverbatim
//
// \param[in] sparam
// \verbatim
//          sparam is REAL array, dimension (5)
//     sparam1=sflag
//     sparam(2)=sh11
//     sparam(3)=sh21
//     sparam(4)=sh12
//     sparam(5)=sh22
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
//  =====================================================================
func Srotm(major *byte, n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, sparam *[]float32) {
	var sflag, sh11, sh12, sh21, sh22, w, z float32
	var i, kx, ky, nsteps int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	sflag = (*sparam)[1-1]
	if *n <= 0 || (sflag+2.0 == 0.0) {
		return
	}
	if *incx == *incy && *incx > 0 {
		//
		nsteps = (*n) * (*incx)
		if sflag < 0.0 {
			sh11 = (*sparam)[1]
			sh12 = (*sparam)[3]
			sh21 = (*sparam)[2]
			sh22 = (*sparam)[4]
			for i = 1; i <= nsteps; i += *incx {
				w = (*sx)[i-1]
				z = (*sy)[i-1]
				(*sx)[i-1] = w*sh11 + z*sh12
				(*sy)[i-1] = w*sh21 + z*sh22
			}
		} else if sflag == 0.0 {
			sh12 = (*sparam)[3]
			sh21 = (*sparam)[2]
			for i = 1; i <= nsteps; i += *incx {
				w = (*sx)[i-1]
				z = (*sy)[i-1]
				(*sx)[i-1] = w + z*sh12
				(*sy)[i-1] = w*sh21 + z
			}
		} else {
			sh11 = (*sparam)[1]
			sh22 = (*sparam)[4]
			for i = 1; i <= nsteps; i += *incx {
				w = (*sx)[i-1]
				z = (*sy)[i-1]
				(*sx)[i-1] = w*sh11 + z
				(*sy)[i-1] = -w + sh22*z
			}
		}
	} else {
		kx = 1
		ky = 1
		if *incx < 0 {
			kx = 1 + (1-(*n))*(*incx)
		}
		if *incy < 0 {
			ky = 1 + (1-(*n))*(*incy)
		}

		if sflag < 0.0 {
			sh11 = (*sparam)[1]
			sh12 = (*sparam)[3]
			sh21 = (*sparam)[2]
			sh22 = (*sparam)[4]
			for i = 1; i <= *n; i++ {
				w = (*sx)[kx-1]
				z = (*sy)[ky-1]
				(*sx)[kx-1] = w*sh11 + z*sh12
				(*sy)[ky-1] = w*sh21 + z*sh22
				kx += *incx
				ky += *incy
			}
		} else if sflag == 0.0 {
			sh12 = (*sparam)[3]
			sh21 = (*sparam)[2]
			for i = 1; i <= *n; i++ {
				w = (*sx)[kx-1]
				z = (*sy)[ky-1]
				(*sx)[kx-1] = w + z*sh12
				(*sy)[ky-1] = w*sh21 + z
				kx += *incx
				ky += *incy
			}
		} else {
			sh11 = (*sparam)[1]
			sh22 = (*sparam)[4]
			for i = 1; i <= *n; i++ {
				w = (*sx)[kx-1]
				z = (*sy)[ky-1]
				(*sx)[kx-1] = w*sh11 + z
				(*sy)[ky-1] = -w + sh22*z
				kx += *incx
				ky += *incy
			}
		}
	}
	return
}
