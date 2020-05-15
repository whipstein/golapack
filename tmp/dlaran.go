package goblas

import 

// Dlaran returns a random real number from a uniform (0,1)
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
//       DOUBLE PRECISION FUNCTION Dlaran( iseed)
//
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
// Dlaran returns a random real number from a uniform (0,1)
// distribution.
// \endverbatim
//
//  Arguments:
//  ==========
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
// \ingroup list_temp
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  This routine uses a multiplicative congruential method with modulus
//  2**48 and multiplier 33952834046453 (see G.S.Fishman,
//  'Multiplicative congruential random number generators with modulus
//  2**b: an exhaustive analysis for b = 32 and a partial analysis for
//  b = 48', Math. Comp. 189, pp 331-344, 1990).
//
//  48-bit integers are stored in 4 integer array elements with 12 bits
//  per element. Hence the routine is portable across machines with
//  integers of 32 bits or more.
// \endverbatim
//
//  =====================================================================
func Dlaran(iseed *[]int) (dlaranReturn *float64) {
	dlaranReturn = new(float64)
	m1 := new(int)
	m2 := new(int)
	m3 := new(int)
	m4 := new(int)
	one := new(float64)
	ipw2 := new(int)
	r := new(float64)
	it1 := new(int)
	it2 := new(int)
	it3 := new(int)
	it4 := new(int)
	rndout := new(float64)
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*m1) = 494
	(*m2) = 322
	(*m3) = 2508
	(*m4) = 2549
	(*one) = 1.0e+0
	(*ipw2) = 4096
	(*r) = (*one) / (*ipw2)
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
Label10:
	;
	//
	//     multiply the seed by the multiplier modulo 2**48
	//
	(*it4) = (*(iseed))[3] * (*m4)
	(*it3) = (*it4) / (*ipw2)
	(*it4) = (*it4) - (*ipw2)*(*it3)
	(*it3) = (*it3) + (*(iseed))[2]*(*m4) + (*(iseed))[3]*(*m3)
	(*it2) = (*it3) / (*ipw2)
	(*it3) = (*it3) - (*ipw2)*(*it2)
	(*it2) = (*it2) + (*(iseed))[1]*(*m4) + (*(iseed))[2]*(*m3) + (*(iseed))[3]*(*m2)
	(*it1) = (*it2) / (*ipw2)
	(*it2) = (*it2) - (*ipw2)*(*it1)
	(*it1) = (*it1) + (*(iseed))[0]*(*m4) + (*(iseed))[1]*(*m3) + (*(iseed))[2]*(*m2) + (*(iseed))[3]*(*m1)
	(*it1) = (MOD((*it1), (*ipw2)))
	//
	//     return updated seed
	//
	(*(iseed))[0] = (*it1)
	(*(iseed))[1] = (*it2)
	(*(iseed))[2] = (*it3)
	(*(iseed))[3] = (*it4)
	//
	//     convert 48-bit integer to a real number in the interval (0,1)
	//
	(*rndout) = (*r) * (DBLE((*it1)) + (*r)*(DBLE((*it2))+(*r)*(DBLE((*it3))+(*r)*(DBLE((*it4))))))
	//
	if (*rndout) == 1.0e+0 {
		//        If a real number has n bits of precision, and the first
		//        n bits of the 48-bit integer above happen to be all 1 (which
		//        will occur about once every 2**n calls), then Dlaran will
		//        be rounded to exactly 1.0.
		//        Since Dlaran is not supposed to return exactly 0.0 or 1.0
		//        (and some callers of Dlaran, such as CLARND, depend on that),
		//        the statistically correct thing to do in this situation is
		//        simply to iterate again.
		//        N.B. the case Dlaran = 0.0 should not be possible.
		//
		goto Label10
	}
	//
	(*(dlaranReturn)) = (*rndout)
	return
	//
	//     End of Dlaran
	//
}
