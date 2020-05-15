package goblas

// \brief \b Izamax
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION Izamax(N,ZX,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 ZX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Izamax finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
// \param[in] ZX
// \verbatim
//          ZX is COMPLEX//16 array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of ZX
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
// \ingroup aux_blas
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//     jack dongarra, 1/15/85.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Izamax(n *int, zx *[]complex128, incx *int) (izamaxReturn *int) {
	izamaxreturn := new(int)
	dmax := new(float64)
	i := new(int)
	ix := new(int)
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
	//*     .. External Functions ..
	//*     ..
	(*izamax) = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	(*izamax) = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		(*dmax) = (*Dcabs1(&((*zx)[0])))
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if Dcabs1(&((*zx)[(*i)-1])) > (*dmax) {
				(*izamax) = (*i)
				(*dmax) = (*Dcabs1(&((*zx)[(*i)-1])))
			}
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*ix) = 1
		(*dmax) = (*Dcabs1(&((*zx)[0])))
		(*ix) = (*ix) + (*incx)
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if Dcabs1(&((*zx)[(*ix)-1])) > (*dmax) {
				(*izamax) = (*i)
				(*dmax) = (*Dcabs1(&((*zx)[(*ix)-1])))
			}
			(*ix) = (*ix) + (*incx)
		}
	}
	return
}
