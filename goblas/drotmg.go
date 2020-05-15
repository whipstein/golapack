package goblas

import (
	"math"
)

// Drotmg ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drotmg(dd1,dd2,dx1,dy1,dparam)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION dd1,dd2,dx1,dy1
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION dparam(5)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
//    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRTdd1*dx1,DSQRTdd2*>    DY2)**T.
//    WITH dparam1=dflag, H HAS one OF THE FOLLOWING FORMS..
//
//    dflag=-1.D0     dflag=0.D0        dflag=1.D0     dflag=-2.D0
//
//      (dh11  dh12)    (1.D0  dh12)    (dh11  1.D0)    (1.D0  0.D0)
//    H=(          )    (          )    (          )    (          )
//      (dh21  dh22),   (dh21  1.D0),   (-1.D0 dh22),   (0.D0  1.D0).
//    LOCATIONS 2-4 OF dparam CONTAIN dh11, dh21, dh12, AND dh22
//    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
//    VALUE OF dparam1 ARE NOT STORED IN dparam.)
//
//    THE VALUES OF gamsq AND rgamsq SET IN THE DATA STATEMENT MAY BE
//    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
//    OF dd1 AND dd2.  ALL ACTUAL SCALING OF DATA IS DONE USING gam.
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in,out] dd1
// \verbatim
//          dd1 is DOUBLE PRECISION
// \endverbatim
//
// \param[in,out] dd2
// \verbatim
//          dd2 is DOUBLE PRECISION
// \endverbatim
//
// \param[in,out] dx1
// \verbatim
//          dx1 is DOUBLE PRECISION
// \endverbatim
//
// \param[in] dy1
// \verbatim
//          dy1 is DOUBLE PRECISION
// \endverbatim
//
// \param[out] dparam
// \verbatim
//          dparam is DOUBLE PRECISION array, dimension (5)
//     dparam1=dflag
//     dparam(2)=dh11
//     dparam(3)=dh21
//     dparam(4)=dh12
//     dparam(5)=dh22
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
//  =====================================================================
func Drotmg(major *byte, dd1 *float64, dd2 *float64, dx1 *float64, dy1 *float64, dparam *[]float64) {
	var dflag, dh11, dh12, dh21, dh22, dp1, dp2, dq1, dq2, dtemp, du, gam, gamsq, rgamsq float64
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	gam, gamsq, rgamsq = 4096, 16777216, 5.9604645e-8
	//     ..
	if *dd1 < 0.0 {
		//        GO zero-H-D-AND-dx1..
		dflag = -1.0

		*dd1 = 0.0
		*dd2 = 0.0
		*dx1 = 0.0
	} else {
		//        CASE-dd1-NONNEGATIVE
		dp2 = (*dd2) * (*dy1)
		if dp2 == 0.0 {
			dflag = -2.0
			(*dparam)[0] = dflag
			return
		}
		//        REGULAR-CASE..
		dp1 = (*dd1) * (*dx1)
		dq2 = dp2 * (*dy1)
		dq1 = dp1 * (*dx1)
		//
		if math.Abs(dq1) > math.Abs(dq2) {
			dh21 = -(*dy1) / (*dx1)
			dh12 = dp2 / dp1
			//
			du = 1.0 - dh12*dh21
			//
			if du > 0.0 {
				dflag = 0.0
				*dd1 = (*dd1) / du
				*dd2 = (*dd2) / du
				*dx1 = (*dx1) * du
			}
		} else {
			if dq2 < 0.0 {
				//              GO zero-H-D-AND-dx1..
				dflag = -1.0

				*dd1 = 0.0
				*dd2 = 0.0
				*dx1 = 0.0
			} else {
				dflag = 1.0
				dh11 = dp1 / dp2
				dh22 = (*dx1) / (*dy1)
				du = 1.0 + dh11*dh22
				dtemp = (*dd2) / du
				*dd2 = (*dd1) / du
				*dd1 = dtemp
				*dx1 = (*dy1) * du
			}
		}
		//     PROCEDURE..scale-CHECK
		if *dd1 != 0.0 {
			for (*dd1 <= rgamsq) || (*dd1 >= gamsq) {
				if dflag == 0.0 {
					dh11 = 1.0
					dh22 = 1.0
					dflag = -1.0
				} else {
					dh21 = -1.0
					dh12 = 1.0
					dflag = -1.0
				}
				if *dd1 <= rgamsq {
					*dd1 *= math.Pow(gam, 2)
					*dx1 /= gam
					dh11 /= gam
					dh12 /= gam
				} else {
					*dd1 = (*dd1) / math.Pow(gam, 2)
					*dx1 = (*dx1) * gam
					dh11 *= gam
					dh12 *= gam
				}
			}
		}
		if *dd2 != 0.0 {
			for (math.Abs(*dd2) <= rgamsq) || (math.Abs(*dd2) >= gamsq) {
				if dflag == 0.0 {
					dh11 = 1.0
					dh22 = 1.0
					dflag = -1.0
				} else {
					dh21 = -1.0
					dh12 = 1.0
					dflag = -1.0
				}
				if math.Abs(*dd2) <= rgamsq {
					*dd2 *= math.Pow(gam, 2)
					dh21 = dh21 / gam
					dh22 = dh22 / gam
				} else {
					*dd2 = (*dd2) / math.Pow(gam, 2)
					dh21 = dh21 * gam
					dh22 = dh22 * gam
				}
			}
		}
	}
	if dflag < 0.0 {
		(*dparam)[1] = dh11
		(*dparam)[2] = dh21
		(*dparam)[3] = dh12
		(*dparam)[4] = dh22
	} else if dflag == 0.0 {
		(*dparam)[2] = dh21
		(*dparam)[3] = dh12
	} else {
		(*dparam)[1] = dh11
		(*dparam)[4] = dh22
	}
	(*dparam)[0] = dflag
	return
}
