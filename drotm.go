package golapack

// Drotm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
//
//       .. Scalar Arguments ..
//       INTEGER INCX,INCY,N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
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
//    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
//    (DY**T)
//
//    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
//    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
//    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
//
//    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
//
//      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
//    H=(          )    (          )    (          )    (          )
//      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
//    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
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
//          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//         storage spacing between elements of DX
// \endverbatim
//
// \param[in,out] DY
// \verbatim
//          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// \endverbatim
//
// \param[in] INCY
// \verbatim
//          INCY is INTEGER
//         storage spacing between elements of DY
// \endverbatim
//
// \param[in] DPARAM
// \verbatim
//          DPARAM is DOUBLE PRECISION array, dimension (5)
//     DPARAM(1)=DFLAG
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
func Drotm(n *int, dx *[]float64, dxoff, incx *int, dy *[]float64, dyoff, incy *int, dparam *[]float64) {
	var dflag, dh11, dh12, dh21, dh22, w, z float64
	var i, kx, ky, nsteps int
	var zero float64 = 0.0
	var two float64 = 2.0

	dflag = (*dparam)[0]
	if (*n) <= 0 || (dflag+two == zero) {
		return
	}
	if (*incx) == (*incy) && (*incx) > 0 {

		nsteps = (*n) * (*incx)
		if dflag < zero {
			dh11 = (*dparam)[1]
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			dh22 = (*dparam)[4]
			for i = 1; i <= nsteps; i += (*incx) {
				w = (*dx)[i-1+(*dxoff)]
				z = (*dy)[i-1+(*dyoff)]
				(*dx)[i-1+(*dxoff)] = w*dh11 + z*dh12
				(*dy)[i-1+(*dyoff)] = w*dh21 + z*dh22
			}
		} else if dflag == zero {
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			for i = 1; i <= nsteps; i += (*incx) {
				w = (*dx)[i-1+(*dxoff)]
				z = (*dy)[i-1+(*dyoff)]
				(*dx)[i-1+(*dxoff)] = w + z*dh12
				(*dy)[i-1+(*dyoff)] = w*dh21 + z
			}
		} else {
			dh11 = (*dparam)[1]
			dh22 = (*dparam)[4]
			for i = 1; i <= nsteps; i += (*incx) {
				w = (*dx)[i-1+(*dxoff)]
				z = (*dy)[i-1+(*dyoff)]
				(*dx)[i-1+(*dxoff)] = w*dh11 + z
				(*dy)[i-1+(*dyoff)] = -w + dh22*z
			}
		}
	} else {
		kx = 1
		ky = 1
		if (*incx) < 0 {
			kx = 1 + (1-(*n))*(*incx)
		}
		if (*incy) < 0 {
			ky = 1 + (1-(*n))*(*incy)
		}

		if dflag < zero {
			dh11 = (*dparam)[1]
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			dh22 = (*dparam)[4]
			for i = 1; i <= (*n); i++ {
				w = (*dx)[kx-1+(*dxoff)]
				z = (*dy)[ky-1+(*dyoff)]
				(*dx)[kx-1+(*dxoff)] = w*dh11 + z*dh12
				(*dy)[ky-1+(*dyoff)] = w*dh21 + z*dh22
				kx = kx + (*incx)
				ky = ky + (*incy)
			}
		} else if dflag == zero {
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			for i = 1; i <= (*n); i++ {
				w = (*dx)[kx-1+(*dxoff)]
				z = (*dy)[ky-1+(*dyoff)]
				(*dx)[kx-1+(*dxoff)] = w + z*dh12
				(*dy)[ky-1+(*dyoff)] = w*dh21 + z
				kx = kx + (*incx)
				ky = ky + (*incy)
			}
		} else {
			dh11 = (*dparam)[1]
			dh22 = (*dparam)[4]
			for i = 1; i <= (*n); i++ {
				w = (*dx)[kx-1+(*dxoff)]
				z = (*dy)[ky-1+(*dyoff)]
				(*dx)[kx-1+(*dxoff)] = w*dh11 + z
				(*dy)[ky-1+(*dyoff)] = -w + dh22*z
				kx = kx + (*incx)
				ky = ky + (*incy)
			}
		}
	}
}
