package goblas

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
//       SUBROUTINE Drotm(n,dx,incx,dy,incy,dparam)
//
//       .. Scalar Arguments ..
//       INTEGER incx,incy,n
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION dparam(5),dx(*),dy(*)
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
//    (dx**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF dx ARE IN
//    (dy**T)
//
//    dx(LX+i*incx), i = 0 TO n-1, WHERE LX = 1 IF incx .GE. 0, ELSE
//    LX = (-incx)*n, AND SIMILARLY FOR sy USING LY AND incy.
//    WITH dparam1=dflag, H HAS one OF THE FOLLOWING FORMS..
//
//    dflag=-1.D0     dflag=0.D0        dflag=1.D0     dflag=-2.D0
//
//      (dh11  dh12)    (1.D0  dh12)    (dh11  1.D0)    (1.D0  0.D0)
//    H=(          )    (          )    (          )    (          )
//      (dh21  dh22),   (dh21  1.D0),   (-1.D0 dh22),   (0.D0  1.D0).
//    SEE Drotmg FOR a DESCRIPTION OF DATA STORAGE IN dparam.
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
// \param[in,out] dx
// \verbatim
//          dx is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incx ) )
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of dx
// \endverbatim
//
// \param[in,out] dy
// \verbatim
//          dy is DOUBLE PRECISION array, dimension ( 1 + ( n - 1 )*abs( incy ) )
// \endverbatim
//
// \param[in] incy
// \verbatim
//          incy is INTEGER
//         storage spacing between elements of dy
// \endverbatim
//
// \param[in] dparam
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
func Drotm(major *byte, n *int, dx *[]float64, incx *int, dy *[]float64, incy *int, dparam *[]float64) {
	var dflag, dh11, dh12, dh21, dh22, w, z float64
	var i, kx, ky, nsteps int
	//
	//  -- Reference BLAS level1 routine (version 3.8.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	dflag = (*dparam)[1-1]
	if *n <= 0 || (dflag+2.0 == 0.0) {
		return
	}
	if *incx == *incy && *incx > 0 {
		//
		nsteps = (*n) * (*incx)
		if dflag < 0.0 {
			dh11 = (*dparam)[1]
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			dh22 = (*dparam)[4]
			for i = 1; i <= nsteps; i += *incx {
				w = (*dx)[i-1]
				z = (*dy)[i-1]
				(*dx)[i-1] = w*dh11 + z*dh12
				(*dy)[i-1] = w*dh21 + z*dh22
			}
		} else if dflag == 0.0 {
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			for i = 1; i <= nsteps; i += *incx {
				w = (*dx)[i-1]
				z = (*dy)[i-1]
				(*dx)[i-1] = w + z*dh12
				(*dy)[i-1] = w*dh21 + z
			}
		} else {
			dh11 = (*dparam)[1]
			dh22 = (*dparam)[4]
			for i = 1; i <= nsteps; i += *incx {
				w = (*dx)[i-1]
				z = (*dy)[i-1]
				(*dx)[i-1] = w*dh11 + z
				(*dy)[i-1] = -w + dh22*z
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
		//
		if dflag < 0.0 {
			dh11 = (*dparam)[1]
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			dh22 = (*dparam)[4]
			for i = 1; i <= *n; i++ {
				w = (*dx)[kx-1]
				z = (*dy)[ky-1]
				(*dx)[kx-1] = w*dh11 + z*dh12
				(*dy)[ky-1] = w*dh21 + z*dh22
				kx += *incx
				ky += *incy
			}
		} else if dflag == 0.0 {
			dh12 = (*dparam)[3]
			dh21 = (*dparam)[2]
			for i = 1; i <= *n; i++ {
				w = (*dx)[kx-1]
				z = (*dy)[ky-1]
				(*dx)[kx-1] = w + z*dh12
				(*dy)[ky-1] = w*dh21 + z
				kx += *incx
				ky += *incy
			}
		} else {
			dh11 = (*dparam)[1]
			dh22 = (*dparam)[4]
			for i = 1; i <= *n; i++ {
				w = (*dx)[kx-1]
				z = (*dy)[ky-1]
				(*dx)[kx-1] = w*dh11 + z
				(*dy)[ky-1] = -w + dh22*z
				kx += *incx
				ky += *incy
			}
		}
	}
	return
}
