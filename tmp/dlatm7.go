package goblas

import (
	"math"
)

//    dlatm7 computes the entries of D as specified by mode
//    cond and irsign. idist and iseed determine the generation
//    of random numbers. dlatm7 is called by DLAtmT to generate
//    random test matrices.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE dlatm7( mode, cond, irsign, idist, iseed, d, n,
//                          rank, info)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   cond
//       inTEGER            idist, info, irsign, mode, n, rank
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   d(*)
//       inTEGER            iseed( 4)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    dlatm7 computes the entries of D as specified by mode
//    cond and irsign. idist and iseed determine the generation
//    of random numbers. dlatm7 is called by DLAtmT to generate
//    random test matrices.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \verbatim
//  mode   - inTEGER
//           On entry describes how D is to be computed:
//           mode = 0 means do not change D.
//
//           mode = 1 sets d1=1 and d(2:rank)=1.0/cond
//           mode = 2 sets d(1:rank-1)=1 and d(rank)=1.0/cond
//           mode = 3 sets d(i)=cond**(-(I-1)/(rank-1)) I=1:rank
//
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
//
//  cond   - DOUBLE PRECISION
//           On entry, used as described under mode above.
//           If used, it must be >= 1. Not modified.
//
//  irsign - inTEGER
//           On entry, if mode neither -6, 0 nor 6, determines sign of
//           entries of D
//           0 => leave entries of D unchanged
//           1 => multiply each entry of D by 1 or -1 with probability .5
//
//  idist  - CHARACTER*1
//           On entry, idist specifies the type of distribution to be
//           used to generate a random matrix .
//           1 => UNIFORM( 0, 1)
//           2 => UNIFORM( -1, 1)
//           3 => normaL( 0, 1)
//           Not modified.
//
//  iseed  - inTEGER array, dimension ( 4)
//           On entry iseed specifies the seed of the random number
//           generator. The random number generator uses a
//           linear congruential sequence limited to small
//           integers, and so should produce machine independent
//           random numbers. The values of iseed are changed on
//           exit, and can be used in the next call to dlatm7
//           to continue the same random number sequence.
//           Changed on exit.
//
//  D      - DOUBLE PRECISION array, dimension ( Min( M, N))
//           Array to be computed according to mode, cond and irsign.
//           May be changed on exit if mode is nonzero.
//
//  N      - inTEGER
//           Number of entries of D. Not modified.
//
//  rank   - inTEGER
//           The rank of matrix to be generated for modes 1,2,3 only.
//           d( rank+1:N) = 0.
//           Not modified.
//
//  info   - inTEGER
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
func dlatm7(mode *int, cond *float64, irsign *int, idist *int, iseed *[]int, d *[]float64, n *int, rank *int, info *int) {
	one := new(float64)
	zero := new(float64)
	half := new(float64)
	alpha := new(float64)
	temp := new(float64)
	i := new(int)
	//
	//  -- lapACK computational routine (version 3.7.0) --
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
	(*zero) = 0.0
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
		Xerbla(func() *[]byte {y := []byte("dlatm7"); return &y }(), -(*(info)))
		return
	}
	//
	//     Compute D according to cond and mode
	//
	if (*(mode)) != 0 {
		switch ABS((*(mode))) {
		case 1:
			goto Label100
		case 2:
			goto Label130
		case 3:
			goto Label160
		case 4:
			goto Label190
		case 5:
			goto Label210
		case 6:
			goto Label230
		}
		//
		//        one large D value:
		//
	Label100:
		;
		for (*i) = 2; (*i) <= (*(rank)); (*i)++ {
			(*(d))[(*i)-1] = (*one) / (*(cond))
			//Label110:
		}
		for (*i) = (*(rank)) + 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-1] = (*zero)
			//Label120:
		}
		(*(d))[0] = (*one)
		goto Label240
		//
		//        one small D value:
		//
	Label130:
		;
		for (*i) = 1; (*i) <= (*(rank))-1; (*i)++ {
			(*(d))[(*i)-1] = (*one)
			//Label140:
		}
		for (*i) = (*(rank)) + 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-1] = (*zero)
			//Label150:
		}
		(*(d))[(*(rank))-1] = (*one) / (*(cond))
		goto Label240
		//
		//        Exponentially distributed D values:
		//
	Label160:
		;
		(*(d))[0] = (*one)
		if (*(n)) > 1 && (*(rank)) > 1 {
			(*alpha) = math.Pow((*(cond)), (-(*one) / DBLE((*(rank))-1)))
			for (*i) = 2; (*i) <= (*(rank)); (*i)++ {
				(*(d))[(*i)-1] = math.Pow((*alpha), ((*i) - 1))
				//Label170:
			}
			for (*i) = (*(rank)) + 1; (*i) <= (*(n)); (*i)++ {
				(*(d))[(*i)-1] = (*zero)
				//Label180:
			}
		}
		goto Label240
		//
		//        Arithmetically distributed D values:
		//
	Label190:
		;
		(*(d))[0] = (*one)
		if (*(n)) > 1 {
			(*temp) = (*one) / (*(cond))
			(*alpha) = ((*one) - (*temp)) / DBLE((*(n))-1)
			for (*i) = 2; (*i) <= (*(n)); (*i)++ {
				(*(d))[(*i)-1] = DBLE((*(n))-(*i))*(*alpha) + (*temp)
				//Label200:
			}
		}
		goto Label240
		//
		//        Randomly distributed D values on ( 1/cond, 1):
		//
	Label210:
		;
		(*alpha) = (*LOG((*one) / (*(cond))))
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(d))[(*i)-1] = (*EXP((*alpha) * Dlaran((iseed))))
			//Label220:
		}
		goto Label240
		//
		//        Randomly distributed D values from idist
		//
	Label230:
		;
		Dlarnv((idist), (iseed), (n), (d))
		//
	Label240:
		;
		//
		//        If mode neither -6 nor 0 nor 6, and irsign = 1, assign
		//        random signs to D
		//
		if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && (*(irsign)) == 1 {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*temp) = (*Dlaran((iseed)))
				if (*temp) > (*half) {
					(*(d))[(*i)-1] = -(*(d))[(*i)-1]
				}
				//Label250:
			}
		}
		//
		//        Reverse if mode < 0
		//
		if (*(mode)) < 0 {
			for (*i) = 1; (*i) <= (*(n))/2; (*i)++ {
				(*temp) = (*(d))[(*i)-1]
				(*(d))[(*i)-1] = (*(d))[(*(n))+1-(*i)-1]
				(*(d))[(*(n))+1-(*i)-1] = (*temp)
				//Label260:
			}
		}
		//
	}
	//
	return
	//
	//     End of dlatm7
	//
}
