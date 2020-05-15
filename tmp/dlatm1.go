package goblas

import (
	"math"
)

//    Dlatm1 computes the entries of d(1..N) as specified by
//    mode, cond and irsign. idist and iseed determine the generation
//    of random numbers. Dlatm1 is called by DLAtmR to generate
//    random test matrices for lapACK programs.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlatm1( mode, cond, irsign, idist, iseed, d, n, info)
//
//       .. Scalar Arguments ..
//       intEGER            idist, info, irsign, mode, N
//       DOUBLE PRECISION   cond
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   d(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Dlatm1 computes the entries of d(1..N) as specified by
//    mode, cond and irsign. idist and iseed determine the generation
//    of random numbers. Dlatm1 is called by DLAtmR to generate
//    random test matrices for lapACK programs.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] mode
// \verbatim
//          mode is intEGER
//           On entry describes how D is to be computed:
//           mode = 0 means do not change D.
//           mode = 1 sets d(1)=1 and d(2:N)=1.0/cond
//           mode = 2 sets d(1:N-1)=1 and d(n)=1.0/cond
//           mode = 3 sets d(i)=cond**(-(I-1)/(N-1))
//           mode = 4 sets d(i)=1 - (i-1)/(N-1)*(1 - 1/cond)
//           mode = 5 sets D to random numbers in the range
//                    ( 1/cond, 1) such that their logarithms
//                    are uniformly distributed.
//           mode = 6 set D to random numbers from same distribution
//                    as the rest of the matrix.
//           mode < 0 has the same meaning as ABS(mode), except that
//              the order of the elements of D is reversed.
//           Thus if mode is positive, D has entries ranging from
//              1 to 1/cond, if negative, from 1/cond to 1,
//           Not modified.
// \endverbatim
//
// \param[in] cond
// \verbatim
//          cond is DOUBLE PRECISION
//           On entry, used as described under mode above.
//           If used, it must be >= 1. Not modified.
// \endverbatim
//
// \param[in] irsign
// \verbatim
//          irsign is intEGER
//           On entry, if mode neither -6, 0 nor 6, determines sign of
//           entries of D
//           0 => leave entries of D unchanged
//           1 => multiply each entry of D by 1 or -1 with probability .5
// \endverbatim
//
// \param[in] idist
// \verbatim
//          idist is intEGER
//           On entry, idist specifies the type of distribution to be
//           used to generate a random matrix .
//           1 => UNIFORM( 0, 1)
//           2 => UNIFORM( -1, 1)
//           3 => normaL( 0, 1)
//           Not modified.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is intEGER array, dimension ( 4)
//           On entry iseed specifies the seed of the random number
//           generator. The random number generator uses a
//           linear congruential sequence limited to small
//           integers, and so should produce machine independent
//           random numbers. The values of iseed are changed on
//           exit, and can be used in the next call to Dlatm1
//           to continue the same random number sequence.
//           Changed on exit.
// \endverbatim
//
// \param[in,out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension ( N)
//           Array to be computed according to mode, cond and irsign.
//           May be changed on exit if mode is nonzero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//           Number of entries of D. Not modified.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
//            0  => normal termination
//           -1  => if mode not in range -6 to 6
//           -2  => if mode neither -6, 0 nor 6, and
//                  irsign neither 0 nor 1
//           -3  => if mode neither -6, 0 nor 6 and cond less than 1
//           -4  => if mode equals 6 or -6 and idist not in range 1 to 3
//           -7  => if N negative
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
//  =====================================================================
func Dlatm1(mode *int, cond *float64, irsign *int, idist *int, iseed *[]int, d *[]float64, n *int, info *int) {
	one := new(float64)
	half := new(float64)
	i := new(int)
	alpha := new(float64)
	temp := new(float64)
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
	(*one) = 1.0
	(*half) = 0.5
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Decode and Test the input parameters. Initialize flags & seed.
	//
	(*(info)) = 0
	//
	//     Quick return if possible
	//
	if (*(n)) == 0 {
		return
	}
	//
	//     Set info if an error
	//
	if (*(mode)) < -6 || (*(mode)) > 6 {
		(*(info)) = -1
	} else if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && ((*(irsign)) != 0 && (*(irsign)) != 1) {
		(*(info)) = -2
	} else if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && (*(cond)) < (*one) {
		(*(info)) = -3
	} else if ((*(mode)) == 6 || (*(mode)) == -6) && ((*(idist)) < 1 || (*(idist)) > 3) {
		(*(info)) = -4
	} else if (*(n)) < 0 {
		(*(info)) = -7
	}
	//
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("Dlatm1"); return &y }(), -(*(info)))
		return
	}
	//
	//     Compute D according to cond and mode
	//
	if (*(mode)) != 0 {
		switch ABS((*(mode))) {
		case 1:
			goto Label10
		case 2:
			goto Label30
		case 3:
			goto Label50
		case 4:
			goto Label70
		case 5:
			goto Label90
		case 6:
			goto Label110
		}
		//
		//        one large D value:
		//
	Label10:
		;
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-(1)] = (*one) / (*(cond))
			//Label20:
		}
		(*(d))[0] = (*one)
		goto Label120
		//
		//        one small D value:
		//
	Label30:
		;
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-(1)] = (*one)
			//Label40:
		}
		(*(d))[(*(n))-(1)] = (*one) / (*(cond))
		goto Label120
		//
		//        Exponentially distributed D values:
		//
	Label50:
		;
		(*(d))[0] = (*one)
		if (*(n)) > 1 {
			(*alpha) = math.Pow((*(cond)), (-(*one) / DBLE((*(n))-1)))
			for (*i) = 2; (*i) <= (*(n)); (*i)++ {
				(*(d))[(*i)-(1)] = math.Pow((*alpha), ((*i) - 1))
				//Label60:
			}
		}
		goto Label120
		//
		//        Arithmetically distributed D values:
		//
	Label70:
		;
		(*(d))[0] = (*one)
		if (*(n)) > 1 {
			(*temp) = (*one) / (*(cond))
			(*alpha) = ((*one) - (*temp)) / DBLE((*(n))-1)
			for (*i) = 2; (*i) <= (*(n)); (*i)++ {
				(*(d))[(*i)-(1)] = DBLE((*(n))-(*i))*(*alpha) + (*temp)
				//Label80:
			}
		}
		goto Label120
		//
		//        Randomly distributed D values on ( 1/cond, 1):
		//
	Label90:
		;
		(*alpha) = (*lOG((*one) / (*(cond))))
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-(1)] = (*EXP((*alpha) * Dlaran((iseed))))
			//Label100:
		}
		goto Label120
		//
		//        Randomly distributed D values from idist
		//
	Label110:
		;
		Dlarnv((idist), (iseed), (n), (d))
		//
	Label120:
		;
		//
		//        If mode neither -6 nor 0 nor 6, and irsign = 1, assign
		//        random signs to D
		//
		if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && (*(irsign)) == 1 {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*temp) = (*Dlaran((iseed)))
				if (*temp) > (*half) {
					(*(d))[(*i)-(1)] = -(*(d))[(*i)-(1)]
				}
				//Label130:
			}
		}
		//
		//        Reverse if mode < 0
		//
		if (*(mode)) < 0 {
			for (*i) = 1; (*i) <= (*(n))/2; (*i)++ {
				(*temp) = (*(d))[(*i)-(1)]
				(*(d))[(*i)-(1)] = (*(d))[(*(n))+1-(*i)-(1)]
				(*(d))[(*(n))+1-(*i)-(1)] = (*temp)
				//Label140:
			}
		}
		//
	}
	//
	return
	//
	//     End of Dlatm1
	//
}
