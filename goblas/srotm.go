package goblas

// \brief \b Srotm
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Srotm(N,SX,incx,SY,incy,SPARAM)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,N
//       ..
//       .. Array Arguments ..
//       REAL SPARAM(5),SX(//),SY(//)
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
//    (SX////T), WHERE ////T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
//    (SX////T)
//
//    SX(LX+I//incx), I = 0 TO N-1, WHERE LX = 1 IF incx .GE. 0, ELSE
//    LX = (-incx)//N, AND SIMILARLY FOR SY USING USING LY AND incy.
//    WITH SPARAM1=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
//
//    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
//
//      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
//    H=(        )    (        )    (        )    (        )
//      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
//    SEE  Srotmg FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
//
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
// \param[in,out] SX
// \verbatim
//          SX is REAL array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of SX
// \endverbatim
//
// \param[in,out] SY
// \verbatim
//          SY is REAL array, dimension ( 1 + ( N - 1)//abs( incy))
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of SY
// \endverbatim
//
// \param[in] SPARAM
// \verbatim
//          SPARAM is REAL array, dimension (5)
//     SPARAM1=SFLAG
//     SPARAM(2)=SH11
//     SPARAM(3)=SH21
//     SPARAM(4)=SH12
//     SPARAM(5)=SH22
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
func Srotm(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, sparam *[]float64) {
	sflag := new(float64)
	sh11 := new(float64)
	sh12 := new(float64)
	sh21 := new(float64)
	sh22 := new(float64)
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
	(*sflag) = (*sparam)[0]
	if (*n) <= 0 || ((*sflag)+(*two) == (*zero)) {
		return
	}
	if (*incx) == (*incy) && (*incx) > 0 {
		//*
		(*nsteps) = (*n) * (*incx)
		if (*sflag) < (*zero) {
			(*sh11) = (*sparam)[1]
			(*sh12) = (*sparam)[3]
			(*sh21) = (*sparam)[2]
			(*sh22) = (*sparam)[4]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*sx)[(*i)-1]
				(*z) = (*sy)[(*i)-1]
				(*sx)[(*i)-1] = (*w)*(*sh11) + (*z)*(*sh12)
				(*sy)[(*i)-1] = (*w)*(*sh21) + (*z)*(*sh22)
			}
		} else if (*sflag) == (*zero) {
			(*sh12) = (*sparam)[3]
			(*sh21) = (*sparam)[2]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*sx)[(*i)-1]
				(*z) = (*sy)[(*i)-1]
				(*sx)[(*i)-1] = (*w) + (*z)*(*sh12)
				(*sy)[(*i)-1] = (*w)*(*sh21) + (*z)
			}
		} else {
			(*sh11) = (*sparam)[1]
			(*sh22) = (*sparam)[4]
			for (*i) = 1; (*i) <= (*nsteps); (*i) += (*incx) {
				(*w) = (*sx)[(*i)-1]
				(*z) = (*sy)[(*i)-1]
				(*sx)[(*i)-1] = (*w)*(*sh11) + (*z)
				(*sy)[(*i)-1] = -(*w) + (*sh22)*(*z)
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
		if (*sflag) < (*zero) {
			(*sh11) = (*sparam)[1]
			(*sh12) = (*sparam)[3]
			(*sh21) = (*sparam)[2]
			(*sh22) = (*sparam)[4]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*sx)[(*kx)-1]
				(*z) = (*sy)[(*ky)-1]
				(*sx)[(*kx)-1] = (*w)*(*sh11) + (*z)*(*sh12)
				(*sy)[(*ky)-1] = (*w)*(*sh21) + (*z)*(*sh22)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		} else if (*sflag) == (*zero) {
			(*sh12) = (*sparam)[3]
			(*sh21) = (*sparam)[2]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*sx)[(*kx)-1]
				(*z) = (*sy)[(*ky)-1]
				(*sx)[(*kx)-1] = (*w) + (*z)*(*sh12)
				(*sy)[(*ky)-1] = (*w)*(*sh21) + (*z)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		} else {
			(*sh11) = (*sparam)[1]
			(*sh22) = (*sparam)[4]
			for (*i) = 1; (*i) <= (*n); (*i)++ {
				(*w) = (*sx)[(*kx)-1]
				(*z) = (*sy)[(*ky)-1]
				(*sx)[(*kx)-1] = (*w)*(*sh11) + (*z)
				(*sy)[(*ky)-1] = -(*w) + (*sh22)*(*z)
				(*kx) = (*kx) + (*incx)
				(*ky) = (*ky) + (*incy)
			}
		}
	}
	return
}
