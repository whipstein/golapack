package goblas

import 
// \brief \b Sscal
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Sscal(N,SA,SX,incx)
//
//       .. Scalar Arguments ..
//       REAL SA
//       INTEGER incx,N
//       ..
//       .. Array Arguments ..
//       REAL SX(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Sscal scales a vector by a constant.
//    uses unrolled loops for increment equal to 1.
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
// \param[in] SA
// \verbatim
//          SA is REAL
//           On entry, SA specifies the scalar alpha.
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
func Sscal(n *int, sa *float64, sx *[]float64, incx *int) {
	i := new(int)
	m := new(int)
	mp1 := new(int)
	nincx := new(int)
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
	if (*n) <= 0 || (*incx) <= 0 {
		return
	}
	if (*incx) == 1 {
		//*
		//*        code for increment equal to 1
		//*
		//*
		//*        clean-up loop
		//*
		(*m) = (MOD((*n), int(5)))
		if (*m) != 0 {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*sx)[(*i)-1] = (*sa) * (*sx)[(*i)-1]
			}
			if (*n) < 5 {
				return
			}
		}
		(*mp1) = (*m) + 1
		for (*i) = (*mp1); (*i) <= (*n); (*i) += 5 {
			(*sx)[(*i)-1] = (*sa) * (*sx)[(*i)-1]
			(*sx)[(*i)+0] = (*sa) * (*sx)[(*i)+0]
			(*sx)[(*i)+1] = (*sa) * (*sx)[(*i)+1]
			(*sx)[(*i)+2] = (*sa) * (*sx)[(*i)+2]
			(*sx)[(*i)+3] = (*sa) * (*sx)[(*i)+3]
		}
	} else {
		//*
		//*        code for increment not equal to 1
		//*
		(*nincx) = (*n) * (*incx)
		for (*i) = 1; (*i) <= (*nincx); (*i) += (*incx) {
			(*sx)[(*i)-1] = (*sa) * (*sx)[(*i)-1]
		}
	}
	return
}
