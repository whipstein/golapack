package goblas

// \brief \b Icamax
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       INTEGER FUNCTION Icamax(N,CX,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       COMPLEX CX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Icamax finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
// \param[in] CX
// \verbatim
//          CX is COMPLEX array, dimension ( 1 + ( N - 1)//abs( incx))
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//         storage spacing between elements of CX
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
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array1 declarations changed to array(//)
// \endverbatim
//
//  =====================================================================
func Icamax(n *int, cx *[]complex128, incx *int) (icamaxReturn *int) {
	icamaxreturn := new(int)
	smax := new(float64)
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
	(*icamax) = 0
	if (*n) < 1 || (*incx) <= 0 {
		return
	}
	(*icamax) = 1
	if (*n) == 1 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		(*smax) = (*Saxpy(&((*cx)[0])))
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if Saxpy(&((*cx)[(*i)-1])) > (*smax) {
				(*icamax) = (*i)
				(*smax) = (*Saxpy(&((*cx)[(*i)-1])))
			}
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*ix) = 1
		(*smax) = (*Saxpy(&((*cx)[0])))
		(*ix) = (*ix) + (*incx)
		for (*i) = 2; (*i) <= (*n); (*i)++ {
			if Saxpy(&((*cx)[(*ix)-1])) > (*smax) {
				(*icamax) = (*i)
				(*smax) = (*Saxpy(&((*cx)[(*ix)-1])))
			}
			(*ix) = (*ix) + (*incx)
		}
	}
	return
}
