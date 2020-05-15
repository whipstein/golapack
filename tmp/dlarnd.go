package goblas

import 

// Dlarnd returns a random real number from a uniform or normal
// distribution.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dlarnd( idist, iseed)
//
//       .. Scalar Arguments ..
//       intEGER            idist
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlarnd returns a random real number from a uniform or normal
// distribution.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] idist
// \verbatim
//          idist is intEGER
//          Specifies the distribution of the random numbers:
//          = 1:  uniform (0,1)
//          = 2:  uniform (-1,1)
//          = 3:  normal (0,1)
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is intEGER array, dimension (4)
//          On entry, the seed of the random number generator; the array
//          elements must be between 0 and 4095, and iseed(4) must be
//          odd.
//          On exit, the seed is updated.
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
// \date December 2016
//
// \ingroup double_matgen
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  This routine calls the auxiliary routine Dlaran to generate a random
//  real number from a uniform (0,1) distribution. The Box-Muller method
//  is used to transform numbers from a uniform to a normal distribution.
// \endverbatim
//
//  =====================================================================
func Dlarnd(idist *int, iseed *[]int) (dlarndReturn *float64) {
	dlarndReturn = new(float64)
	one := new(float64)
	two := new(float64)
	twopi := new(float64)
	t1 := new(float64)
	t2 := new(float64)
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*one) = 1.0e+0
	(*two) = 2.0e+0
	(*twopi) = 6.2831853071795864769252867663e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Generate a real random number from a uniform (0,1) distribution
	//
	(*t1) = (*Dlaran((iseed)))
	//
	if (*(idist)) == 1 {
		//
		//        uniform (0,1)
		//
		(*(dlarndReturn)) = (*t1)
	} else if (*(idist)) == 2 {
		//
		//        uniform (-1,1)
		//
		(*(dlarndReturn)) = (*two)*(*t1) - (*one)
	} else if (*(idist)) == 3 {
		//
		//        normal (0,1)
		//
		(*t2) = (*Dlaran((iseed)))
		(*(dlarndReturn)) = SQRt(-(*two)*lOG(t1)) * COS((*twopi)*(*t2))
	}
	return
	//
	//     End of Dlarnd
	//
}
