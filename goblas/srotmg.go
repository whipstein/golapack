package goblas

// Srotmg ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SROTMG(sd1,sd2,sx1,sy1,sparam)
//
//       .. Scalar Arguments ..
//       REAL sd1,sd2,sx1,sy1
//       ..
//       .. Array Arguments ..
//       REAL sparam(5)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
//    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRTsd1*sx1,SQRTsd2*>    SY2)**T.
//    WITH sparam1=sflag, H HAS 1.0 OF THE FOLLOWING FORMS..
//
//    sflag=-1.E0     sflag=0.E0        sflag=1.E0     sflag=-2.E0
//
//      (sh11  sh12)    (1.E0  sh12)    (sh11  1.E0)    (1.E0  0.E0)
//    H=(          )    (          )    (          )    (          )
//      (sh21  sh22),   (sh21  1.E0),   (-1.E0 sh22),   (0.E0  1.E0).
//    LOCATIONS 2-4 OF sparam CONTAIN sh11,sh21,sh12, AND sh22
//    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
//    VALUE OF sparam1 ARE NOT STORED IN sparam.)
//
//    THE VALUES OF gamsq AND rgamsq SET IN THE DATA STATEMENT MAY BE
//    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
//    OF sd1 AND sd2.  ALL ACTUAL SCALING OF DATA IS DONE USING gam.
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in,out] sd1
// \verbatim
//          sd1 is REAL
// \endverbatim
//
// \param[in,out] sd2
// \verbatim
//          sd2 is REAL
// \endverbatim
//
// \param[in,out] sx1
// \verbatim
//          sx1 is REAL
// \endverbatim
//
// \param[in] sy1
// \verbatim
//          sy1 is REAL
// \endverbatim
//
// \param[out] sparam
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
func Srotmg(major *byte, sd1, sd2, sx1, sy1 *float32, sparam *[]float32) {
	var gam, gamsq, rgamsq, sflag, sh11, sh12, sh21, sh22, sp1, sp2, sq1, sq2, stemp, su float32
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	gam, gamsq, rgamsq = 4096, 1.67772e7, 5.96046e-8

	if *sd1 < 0.0 {
		//        GO 0.0-H-D-AND-sx1..
		sflag = -1.0
		sh11 = 0.0
		sh12 = 0.0
		sh21 = 0.0
		sh22 = 0.0

		*sd1 = 0.0
		*sd2 = 0.0
		*sx1 = 0.0
	} else {
		//        CASE-sd1-NONNEGATIVE
		sp2 = (*sd2) * (*sy1)
		if sp2 == 0.0 {
			sflag = -2.0
			(*sparam)[0] = sflag
			return
		}
		//        REGULAR-CASE..
		sp1 = (*sd1) * (*sx1)
		sq2 = sp2 * (*sy1)
		sq1 = sp1 * (*sx1)
		//
		if (absf32(sq1)) > (absf32(sq2)) {
			sh21 = -(*sy1) / (*sx1)
			sh12 = sp2 / sp1
			//
			su = 1.0 - sh12*sh21
			//
			if su > 0.0 {
				sflag = 0.0
				*sd1 = (*sd1) / su
				*sd2 = (*sd2) / su
				*sx1 = (*sx1) * su
			}
		} else {
			if sq2 < 0.0 {
				//              GO 0.0-H-D-AND-sx1..
				sflag = -1.0
				sh11 = 0.0
				sh12 = 0.0
				sh21 = 0.0
				sh22 = 0.0

				*sd1 = 0.0
				*sd2 = 0.0
				*sx1 = 0.0
			} else {
				sflag = 1.0
				sh11 = sp1 / sp2
				sh22 = (*sx1) / (*sy1)
				su = 1.0 + sh11*sh22
				stemp = (*sd2) / su
				*sd2 = (*sd1) / su
				*sd1 = stemp
				*sx1 = (*sy1) * su
			}
		}
		//     PROCESURE..SCALE-CHECK
		if *sd1 != 0.0 {
			for (*sd1 <= rgamsq) || (*sd1 >= gamsq) {
				if sflag == 0.0 {
					sh11 = 1.0
					sh22 = 1.0
					sflag = -1.0
				} else {
					sh21 = -1.0
					sh12 = 1.0
					sflag = -1.0
				}
				if *sd1 <= rgamsq {
					*sd1 = (*sd1) * powf32(gam, 2)
					*sx1 = (*sx1) / gam
					sh11 = sh11 / gam
					sh12 = sh12 / gam
				} else {
					*sd1 = (*sd1) / powf32(gam, 2)
					*sx1 = (*sx1) * gam
					sh11 = sh11 * gam
					sh12 = sh12 * gam
				}
			}
		}
		if (*sd2) != 0.0 {
			for (absf32((*sd2)) <= rgamsq) || (absf32((*sd2)) >= gamsq) {
				if sflag == 0.0 {
					sh11 = 1.0
					sh22 = 1.0
					sflag = -1.0
				} else {
					sh21 = -1.0
					sh12 = 1.0
					sflag = -1.0
				}
				if absf32(*sd2) <= rgamsq {
					*sd2 = (*sd2) * powf32(gam, 2)
					sh21 = sh21 / gam
					sh22 = sh22 / gam
				} else {
					*sd2 = (*sd2) / powf32(gam, 2)
					sh21 = sh21 * gam
					sh22 = sh22 * gam
				}
			}
		}
	}
	if sflag < 0.0 {
		(*sparam)[1] = sh11
		(*sparam)[2] = sh21
		(*sparam)[3] = sh12
		(*sparam)[4] = sh22
	} else if sflag == 0.0 {
		(*sparam)[2] = sh21
		(*sparam)[3] = sh12
	} else {
		(*sparam)[1] = sh11
		(*sparam)[4] = sh22
	}
	(*sparam)[0] = sflag
	return
}
