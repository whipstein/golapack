package goblas

// \brief \b Drotm
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Drotm(N,DX,incx,DY,incy,DPARAM)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DPARAM(5),DX(//),DY(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
//
//    (DX////T), WHERE ////T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
//    (DY////T)
//
//    DX(LX+I//incx), I = 0 TO N-1, WHERE LX = 1 IF incx .GE. 0, ELSE
//    LX = (-incx)//N, AND SIMILARLY FOR SY USING LY AND incy.
//    WITH DPARAM1=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
//
//    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
//
//      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
//    H=(        )    (        )    (        )    (        )
//      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
//    SEE Drotmg FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
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
// \param[in,out] DX
// \verbatim
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of DX
// \endverbatim
//
// \param[in,out] DY
// \verbatim
//          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of DY
// \endverbatim
//
// \param[in] DPARAM
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
func Drotm(n *int, dx *[]float64, incx *int, dy *[]float64, incy *int, dparam *[]float64) {
	dflag := new(float64)
	dh11 := new(float64)
	dh12 := new(float64)
	dh21 := new(float64)
	dh22 := new(float64)
	two := new(float64)
	w := new(float64)
	z := new(float64)
	zero := new(float64)
	i := new(int)
	kx := new(int)
	ky := new(int)
	nsteps := new(int)
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
	//*     .. Data statements ..
	(*zero), (*two) = 0., 2.
	//*     ..
	//*
	(*dflag) = (*dparam)[0]
	if (*n) <= 0 || ((*dflag)+(*two) == (*zero)) {
		return
	}
	if (*incx) == (*incy) && (*incx) > 0 {
		//*
		(*nsteps) = (*n) * (*incx)
		if (*dflag) < (*zero) {
			(*dh11) = (*dparam)[1]
			(*dh12) = (*dparam)[3]
			(*dh21) = (*dparam)[2]
			(*dh22) = (*dparam)[4]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*dx)[(*i)-1]
				(*z) = (*dy)[(*i)-1]
				(*dx)[(*i)-1] = (*w)*(*dh11) + (*z)*(*dh12)
				(*dy)[(*i)-1] = (*w)*(*dh21) + (*z)*(*dh22)
			}
		} else if (*dflag) == (*zero) {
			(*dh12) = (*dparam)[3]
			(*dh21) = (*dparam)[2]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*dx)[(*i)-1]
				(*z) = (*dy)[(*i)-1]
				(*dx)[(*i)-1] = (*w) + (*z)*(*dh12)
				(*dy)[(*i)-1] = (*w)*(*dh21) + (*z)
			}
		} else {
			(*dh11) = (*dparam)[1]
			(*dh22) = (*dparam)[4]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*dx)[(*i)-1]
				(*z) = (*dy)[(*i)-1]
				(*dx)[(*i)-1] = (*w)*(*dh11) + (*z)
				(*dy)[(*i)-1] = -(*w) + (*dh22)*(*z)
			}
		}
	} else {
		(*kx) = 1
		(*ky) = 1
		if (*incx) < 0 {
			(*kx) = 1 + (1-(*n))*(*incx)
		}
		if (*incy) < 0 {
			(*ky) = 1 + (1-(*n))*(*incy)
		}
		//*
		if (*dflag) < (*zero) {
			(*dh11) = (*dparam)[1]
			(*dh12) = (*dparam)[3]
			(*dh21) = (*dparam)[2]
			(*dh22) = (*dparam)[4]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*dx)[(*kx)-1]
				(*z) = (*dy)[(*ky)-1]
				(*dx)[(*kx)-1] = (*w)*(*dh11) + (*z)*(*dh12)
				(*dy)[(*ky)-1] = (*w)*(*dh21) + (*z)*(*dh22)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		} else if (*dflag) == (*zero) {
			(*dh12) = (*dparam)[3]
			(*dh21) = (*dparam)[2]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*dx)[(*kx)-1]
				(*z) = (*dy)[(*ky)-1]
				(*dx)[(*kx)-1] = (*w) + (*z)*(*dh12)
				(*dy)[(*ky)-1] = (*w)*(*dh21) + (*z)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		} else {
			(*dh11) = (*dparam)[1]
			(*dh22) = (*dparam)[4]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*dx)[(*kx)-1]
				(*z) = (*dy)[(*ky)-1]
				(*dx)[(*kx)-1] = (*w)*(*dh11) + (*z)
				(*dy)[(*ky)-1] = -(*w) + (*dh22)*(*z)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		}
	}
	return
}
