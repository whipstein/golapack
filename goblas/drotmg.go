package goblas

import "math"
import 

// \brief \b Drotmg
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drotmg(DD1,DD2,DX1,DY1,DPARAM)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION DD1,DD2,DX1,DY1
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DPARAM(5)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
//    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)//DX1,DSQRT(DD2)//    DY2)////T.
//    WITH DPARAM1=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
//
//    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
//
//      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
//    H=(        )    (        )    (        )    (        )
//      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
//    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
//    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
//    VALUE OF DPARAM1 ARE NOT STORED IN DPARAM.)
//
//    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
//    INEXACT.  THIS IS ok AS THEY ARE ONLY USED FOR TESTING THE SIZE
//    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in,out] DD1
// \verbatim
//          DD1 is DOUBLE PRECISION
// \endverbatim
//
// \param[in,out] DD2
// \verbatim
//          DD2 is DOUBLE PRECISION
// \endverbatim
//
// \param[in,out] DX1
// \verbatim
//          DX1 is DOUBLE PRECISION
// \endverbatim
//
// \param[in] DY1
// \verbatim
//          DY1 is DOUBLE PRECISION
// \endverbatim
//
// \param[out] DPARAM
// \verbatim
//          DPARAM is DOUBLE PRECISION array, dimension (5)
//     DPARAM1=DFLAG
//     DPARAM(2)=DH11
//     DPARAM(3)=DH21
//     DPARAM(4)=DH12
//     DPARAM(5)=DH22
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
func Drotmg(dd1 *float64, dd2 *float64, dx1 *float64, dy1 *float64, dparam *[]float64) {
	dflag := new(float64)
	dh11 := new(float64)
	dh12 := new(float64)
	dh21 := new(float64)
	dh22 := new(float64)
	dp1 := new(float64)
	dp2 := new(float64)
	dq1 := new(float64)
	dq2 := new(float64)
	dtemp := new(float64)
	du := new(float64)
	gam := new(float64)
	gamsq := new(float64)
	one := new(float64)
	rgamsq := new(float64)
	two := new(float64)
	zero := new(float64)
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
	//*     .. Data statements ..
	//*
	(*zero), (*one), (*two) = 0., 1., 2.
	(*gam), (*gamsq), (*rgamsq) = 4096., 16777216., 5.9604645e-8
	//*     ..
	if (*dd1) < (*zero) {
		//*        GO ZERO-H-D-AND-DX1..
		(*dflag) = -(*one)
		(*dh11) = (*zero)
		(*dh12) = (*zero)
		(*dh21) = (*zero)
		(*dh22) = (*zero)
		//*
		(*dd1) = (*zero)
		(*dd2) = (*zero)
		(*dx1) = (*zero)
	} else {
		//*        CASE-DD1-NONNEGATIVE
		(*dp2) = (*dd2) * (*dy1)
		if (*dp2) == (*zero) {
			(*dflag) = -(*two)
			(*dparam)[0] = (*dflag)
			return
		}
		//*        REGULAR-CASE..
		(*dp1) = (*dd1) * (*dx1)
		(*dq2) = (*dp2) * (*dy1)
		(*dq1) = (*dp1) * (*dx1)
		//*
		if ABS((*dq1)) > dabs(dq2) {
			(*dh21) = -(*dy1) / (*dx1)
			(*dh12) = (*dp2) / (*dp1)
			//*
			(*du) = (*one) - (*dh12)*(*dh21)
			//*
			if (*du) > (*zero) {
				(*dflag) = (*zero)
				(*dd1) = (*dd1) / (*du)
				(*dd2) = (*dd2) / (*du)
				(*dx1) = (*dx1) * (*du)
			}
		} else {
			if (*dq2) < (*zero) {
				//*              GO ZERO-H-D-AND-DX1..
				(*dflag) = -(*one)
				(*dh11) = (*zero)
				(*dh12) = (*zero)
				(*dh21) = (*zero)
				(*dh22) = (*zero)
				//*
				(*dd1) = (*zero)
				(*dd2) = (*zero)
				(*dx1) = (*zero)
			} else {
				(*dflag) = (*one)
				(*dh11) = (*dp1) / (*dp2)
				(*dh22) = (*dx1) / (*dy1)
				(*du) = (*one) + (*dh11)*(*dh22)
				(*dtemp) = (*dd2) / (*du)
				(*dd2) = (*dd1) / (*du)
				(*dd1) = (*dtemp)
				(*dx1) = (*dy1) * (*du)
			}
		}
		//*     PROCEDURE..SCALE-CHECK
		if (*dd1) != (*zero) {
			for ((*dd1) <= (*rgamsq)) || ((*dd1) >= (*gamsq)) {
				if (*dflag) == (*zero) {
					(*dh11) = (*one)
					(*dh22) = (*one)
					(*dflag) = -(*one)
				} else {
					(*dh21) = -(*one)
					(*dh12) = (*one)
					(*dflag) = -(*one)
				}
				if (*dd1) <= (*rgamsq) {
					(*dd1) = (*dd1) * math.pow((*gam), 2)
					(*dx1) = (*dx1) / (*gam)
					(*dh11) = (*dh11) / (*gam)
					(*dh12) = (*dh12) / (*gam)
				} else {
					(*dd1) = (*dd1) / math.pow((*gam), 2)
					(*dx1) = (*dx1) * (*gam)
					(*dh11) = (*dh11) * (*gam)
					(*dh12) = (*dh12) * (*gam)
				}
			}
		}
		if (*dd2) != (*zero) {
			for (dabs((*dd2)) <= (*rgamsq)) || (dabs((*dd2)) >= (*gamsq)) {
				if (*dflag) == (*zero) {
					(*dh11) = (*one)
					(*dh22) = (*one)
					(*dflag) = -(*one)
				} else {
					(*dh21) = -(*one)
					(*dh12) = (*one)
					(*dflag) = -(*one)
				}
				if ABS((*dd2)) <= (*rgamsq) {
					(*dd2) = (*dd2) * math.pow((*gam), 2)
					(*dh21) = (*dh21) / (*gam)
					(*dh22) = (*dh22) / (*gam)
				} else {
					(*dd2) = (*dd2) / math.pow((*gam), 2)
					(*dh21) = (*dh21) * (*gam)
					(*dh22) = (*dh22) * (*gam)
				}
			}
		}
	}
	if (*dflag) < (*zero) {
		(*dparam)[1] = (*dh11)
		(*dparam)[2] = (*dh21)
		(*dparam)[3] = (*dh12)
		(*dparam)[4] = (*dh22)
	} else if (*dflag) == (*zero) {
		(*dparam)[2] = (*dh21)
		(*dparam)[3] = (*dh12)
	} else {
		(*dparam)[1] = (*dh11)
		(*dparam)[4] = (*dh22)
	}
	(*dparam)[0] = (*dflag)
	return
}
