package goblas

import (
	"math"

	"github.com/whipstein/golapack/blas"
)

// Dlattb generates a triangular test matrix in 2-dimensional storage.
// imat and uplo uniquely specify the properties of the test matrix,
// which is returned in the array A.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlattb( imat, uplo, trans, diag, iseed, n, kd, AB,
//                          ldab, b, work, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       inTEGER            imat, info, kd, ldab, N
//       ..
//       .. Array Arguments ..
//       inTEGER            iseed( 4)
//       DOUBLE PRECISION   AB( ldab, *), B(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlattb generates a triangular test matrix in 2-dimensional storage.
// imat and uplo uniquely specify the properties of the test matrix,
// which is returned in the array A.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] imat
// \verbatim
//          imat is inTEGER
//          An integer key describing which matrix to generate for this
//          path.
// \endverbatim
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the matrix A will be upper or lower
//          triangular.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies whether the matrix or its transpose will be used.
//          = 'N':  No transpose
//          = 'T':  Transpose
//          = 'C':  Conjugate transpose (= transpose)
// \endverbatim
//
// \param[out] diag
// \verbatim
//          diag is CHARACTER*1
//          Specifies whether or not the matrix A is unit triangular.
//          = 'N':  Non-unit triangular
//          = 'U':  Unit triangular
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is inTEGER array, dimension (4)
//          The seed vector for the random number generator (used in
//          Dlatms).  Modified on exit.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The order of the matrix to be generated.
// \endverbatim
//
// \param[in] kd
// \verbatim
//          kd is inTEGER
//          The number of superdiagonals or subdiagonals of the banded
//          triangular matrix A.  kd >= 0.
// \endverbatim
//
// \param[out] AB
// \verbatim
//          AB is DOUBLE PRECISION array, dimension (ldab,N)
//          The upper or lower triangular banded matrix a, stored in the
//          first kd+1 rows of AB.  Let j be a column of a, 1<=j<=n.
//          If uplo = 'U', AB(kd+1+i-j,j) = a(i,j) for max(1,j-kd)<=i<=j.
//          If uplo = 'L', AB(1+i-j,j)    = a(i,j) for j<=i<=min(n,j+kd).
// \endverbatim
//
// \param[in] ldab
// \verbatim
//          ldab is inTEGER
//          The leading dimension of the array AB.  ldab >= kd+1.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (n)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (2*n)
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is inTEGER
//          = 0:  successful exit
//          < 0: if info = -k, the k-th argument had an illegal value
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
// \ingroup double_lin
//
//  =====================================================================
func Dlattb(imat *int, uplo *byte, trans *byte, diag *byte, iseed *[]int, n *int, kd *int, ab *[][]float64, ldab *int, b *[]float64, work *[]float64, info *int) {
	one := new(float64)
	two := new(float64)
	zero := new(float64)
	upper := new(bool)
	dist := new(byte)
	packit := new(byte)
	_type := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	ioff := new(int)
	iy := new(int)
	j := new(int)
	jcount := new(int)
	kl := new(int)
	ku := new(int)
	lenj := new(int)
	mode := new(int)
	anorm := new(float64)
	bignum := new(float64)
	bnorm := new(float64)
	bscal := new(float64)
	cndnum := new(float64)
	plus1 := new(float64)
	plus2 := new(float64)
	rexp := new(float64)
	sfac := new(float64)
	smlnum := new(float64)
	star1 := new(float64)
	texp := new(float64)
	tleft := new(float64)
	tnorm := new(float64)
	tscal := new(float64)
	ulp := new(float64)
	unfl := new(float64)
	//
	//  -- lapACK test routine (version 3.7.0) --
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
	(*zero) = 0.0e+0
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
	(*path)[0] = *func() *[]byte {y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y := []byte("TB"); return &y }()
	(*unfl) = (*Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()))
	(*ulp) = Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()) * Dlamch(func() *[]byte {y := []byte("Base"); return &y }())
	(*smlnum) = (*unfl)
	(*bignum) = ((*one) - (*ulp)) / (*smlnum)
	Dlabad(smlnum, bignum)
	if ((*(imat)) >= 6 && (*(imat)) <= 9) || (*(imat)) == 17 {
		(*(diag)) = 'U'
	} else {
		(*(diag)) = 'N'
	}
	(*(info)) = 0
	//
	//     Quick return if N.LE.0.
	//
	if (*(n)) <= 0 {
		return
	}
	//
	//     Call Dlatb4 to set parameters for SLAtmS.
	//
	(*upper) = (*blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()))
	if *upper {
		Dlatb4(path, (imat), (n), (n), _type, kl, ku, anorm, mode, cndnum, dist)
		(*ku) = (*(kd))
		(*ioff) = 1 + MAX(0, (*(kd))-(*(n))+1)
		(*kl) = 0
		(*packit) = 'Q'
	} else {
		Dlatb4(path, -(*(imat)), (n), (n), _type, kl, ku, anorm, mode, cndnum, dist)
		(*kl) = (*(kd))
		(*ioff) = 1
		(*ku) = 0
		(*packit) = 'B'
	}
	//
	//     imat <= 5:  Non-unit triangular matrix
	//
	if (*(imat)) <= 5 {
		Dlatms((n), (n), dist, (iseed), _type, (b), mode, cndnum, anorm, kl, ku, packit, &((*(ab))[(*ioff)-1][0]), (ldab), (work), (info))
		//
		//     imat > 5:  Unit triangular matrix
		//     The diagonal is deliberately set to something other than 1.
		//
		//     imat = 6:  Matrix is the identity
		//
	} else if (*(imat)) == 6 {
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = MAX(1, (*(kd))+2-(*j)); (*i) <= (*(kd)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label10:
				}
				(*(ab))[(*(kd))+0][(*j)-1] = (*j)
				//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(ab))[0][(*j)-1] = (*j)
				for (*i) = 2; (*i) <= (Min((*(kd))+1, (*(n))-(*j)+1)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label30:
				}
				//Label40:
			}
		}
		//
		//     imat > 6:  Non-trivial unit triangular matrix
		//
		//     A unit triangular matrix T with condition cndnum is formed.
		//     In this version, T only has bandwidth 2, the rest of it is zero.
		//
	} else if (*(imat)) <= 9 {
		(*tnorm) = (SQRt((*cndnum)))
		//
		//        Initialize AB to zero.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = MAX(1, (*(kd))+2-(*j)); (*i) <= (*(kd)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label50:
				}
				(*(ab))[(*(kd))+0][(*j)-1] = (DBLE((*j)))
				//Label60:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 2; (*i) <= (Min((*(kd))+1, (*(n))-(*j)+1)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label70:
				}
				(*(ab))[0][(*j)-1] = (DBLE((*j)))
				//Label80:
			}
		}
		//
		//        Special case:  T is tridiagonal.  Set every other offdiagonal
		//        so that the matrix has norm tnorm+1.
		//
		if (*(kd)) == 1 {
			if *upper {
				(*(ab))[0][1] = (*SIGN(tnorm, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*lenj) = ((*(n)) - 3) / 2
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, (work))
				for (*j) = 1; (*j) <= (*lenj); (*j)++ {
					(*(ab))[0][2*((*j)+1)-1] = (*tnorm) * (*(work))[(*j)-1]
					//Label90:
				}
			} else {
				(*(ab))[1][0] = (*SIGN(tnorm, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*lenj) = ((*(n)) - 3) / 2
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, (work))
				for (*j) = 1; (*j) <= (*lenj); (*j)++ {
					(*(ab))[1][2*(*j)+0] = (*tnorm) * (*(work))[(*j)-1]
					//Label100:
				}
			}
		} else if (*(kd)) > 1 {
			//
			//           Form a unit triangular matrix T with condition cndnum.  T is
			//           given by
			//                   | 1   +   *                      |
			//                   |     1   +                      |
			//               T = |         1   +   *              |
			//                   |             1   +              |
			//                   |                 1   +   *      |
			//                   |                     1   +      |
			//                   |                          . . . |
			//        Each element marked with a '*' is formed by taking the product
			//        of the adjacent elements marked with '+'.  The '*'s can be
			//        chosen freely, and the '+'s are chosen so that the inverse of
			//        T will have elements of the same magnitude as T.
			//
			//        The two offdiagonals of T are stored in work.
			//
			(*star1) = (*SIGN(tnorm, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
			(*sfac) = (SQRt((*tnorm)))
			(*plus1) = (*SIGN(sfac, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
			for (*j) = 1; (*j) <= (*(n)); (*j) += 2 {
				(*plus2) = (*star1) / (*plus1)
				(*(work))[(*j)-1] = (*plus1)
				(*(work))[(*(n))+(*j)-1] = (*star1)
				if (*j)+1 <= (*(n)) {
					(*(work))[(*j)+0] = (*plus2)
					(*(work))[(*(n))+(*j)+0] = (*zero)
					(*plus1) = (*star1) / (*plus2)
					//
					//                 Generate a new *-value with norm between sqrt(tnorm)
					//                 and tnorm.
					//
					(*rexp) = (*Dlarnd(func() *int {y := 2; return &y }(), (iseed)))
					if (*rexp) < (*zero) {
						(*star1) = -math.Pow((*sfac), ((*one) - (*rexp)))
					} else {
						(*star1) = math.Pow((*sfac), ((*one) + (*rexp)))
					}
				}
				//Label110:
			}
			//
			//           Copy the tridiagonal T to AB.
			//
			if *upper {
				Dcopy((*(n))-1, (work), func() *int {y := 1; return &y }(), &((*(ab))[(*(kd))-1][1]), (ldab))
				Dcopy((*(n))-2, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(ab))[(*(kd))-0][2]), (ldab))
			} else {
				Dcopy((*(n))-1, (work), func() *int {y := 1; return &y }(), &((*(ab))[1][0]), (ldab))
				Dcopy((*(n))-2, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(ab))[2][0]), (ldab))
			}
		}
		//
		//     imat > 9:  Pathological test cases.  These triangular matrices
		//     are badly scaled or badly conditioned, so when used in solving a
		//     triangular system they may cause overflow in the solution vector.
		//
	} else if (*(imat)) == 10 {
		//
		//        Type 10:  Generate a triangular matrix with elements between
		//        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
		//        Make the right hand side large so that it requires scaling.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j), (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]))
				(*(ab))[(*(kd))+0][(*j)-1] = (*SIGN(two, &((*(ab))[(*(kd))+0][(*j)-1])))
				//Label120:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j)+1, (*(kd))+1))
				if (*lenj) > 0 {
					Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[0][(*j)-1]))
				}
				(*(ab))[0][(*j)-1] = (*SIGN(two, &((*(ab))[0][(*j)-1])))
				//Label130:
			}
		}
		//
		//        Set the right hand side so that the largest value is bignum.
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		(*iy) = (*Idamax((n), (b), func() *int {y := 1; return &y }()))
		(*bnorm) = (ABS(((*(b))[(*iy)-1])))
		(*bscal) = (*bignum) / MAX((*one), (*bnorm))
		Dscal((n), bscal, (b), func() *int {y := 1; return &y }())
		//
	} else if (*(imat)) == 11 {
		//
		//        Type 11:  Make the first diagonal element in the solve small to
		//        cause immediate overflow when dividing by t(j,j).
		//        In type 11, the offdiagonal elements are small (cnorm(j) < 1).
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		(*tscal) = (*one) / DBLE((*(kd))+1)
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j), (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]))
				Dscal((*lenj)-1, tscal, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]), func() *int {y := 1; return &y }())
				(*(ab))[(*(kd))+0][(*j)-1] = (*SIGN(one, &((*(ab))[(*(kd))+0][(*j)-1])))
				//Label140:
			}
			(*(ab))[(*(kd))+0][(*(n))-1] = (*smlnum) * (*(ab))[(*(kd))+0][(*(n))-1]
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j)+1, (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[0][(*j)-1]))
				if (*lenj) > 1 {
					Dscal((*lenj)-1, tscal, &((*(ab))[1][(*j)-1]), func() *int {y := 1; return &y }())
				}
				(*(ab))[0][(*j)-1] = (*SIGN(one, &((*(ab))[0][(*j)-1])))
				//Label150:
			}
			(*(ab))[0][0] = (*smlnum) * (*(ab))[0][0]
		}
		//
	} else if (*(imat)) == 12 {
		//
		//        Type 12:  Make the first diagonal element in the solve small to
		//        cause immediate overflow when dividing by t(j,j).
		//        In type 12, the offdiagonal elements are O1 (cnorm(j) > 1).
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j), (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]))
				(*(ab))[(*(kd))+0][(*j)-1] = (*SIGN(one, &((*(ab))[(*(kd))+0][(*j)-1])))
				//Label160:
			}
			(*(ab))[(*(kd))+0][(*(n))-1] = (*smlnum) * (*(ab))[(*(kd))+0][(*(n))-1]
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j)+1, (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[0][(*j)-1]))
				(*(ab))[0][(*j)-1] = (*SIGN(one, &((*(ab))[0][(*j)-1])))
				//Label170:
			}
			(*(ab))[0][0] = (*smlnum) * (*(ab))[0][0]
		}
		//
	} else if (*(imat)) == 13 {
		//
		//        Type 13:  T is diagonal with small numbers on the diagonal to
		//        make the growth factor underflow, but a small right hand side
		//        chosen so that the solution does not overflow.
		//
		if *upper {
			(*jcount) = 1
			for (*j) = (*(n)); (*j) <= 1; (*j) += -1 {
				for (*i) = MAX(1, (*(kd))+1-((*j)-1)); (*i) <= (*(kd)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label180:
				}
				if (*jcount) <= 2 {
					(*(ab))[(*(kd))+0][(*j)-1] = (*smlnum)
				} else {
					(*(ab))[(*(kd))+0][(*j)-1] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				//Label190:
			}
		} else {
			(*jcount) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 2; (*i) <= (Min((*(n))-(*j)+1, (*(kd))+1)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label200:
				}
				if (*jcount) <= 2 {
					(*(ab))[0][(*j)-1] = (*smlnum)
				} else {
					(*(ab))[0][(*j)-1] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				//Label210:
			}
		}
		//
		//        Set the right hand side alternately zero and small.
		//
		if *upper {
			(*(b))[0] = (*zero)
			for (*i) = (*(n)); (*i) <= 2; (*i) += -2 {
				(*(b))[(*i)-1] = (*zero)
				(*(b))[(*i)-0] = (*smlnum)
				//Label220:
			}
		} else {
			(*(b))[(*(n))-1] = (*zero)
			for (*i) = 1; (*i) <= (*(n))-1; (*i) += 2 {
				(*(b))[(*i)-1] = (*zero)
				(*(b))[(*i)+0] = (*smlnum)
				//Label230:
			}
		}
		//
	} else if (*(imat)) == 14 {
		//
		//        Type 14:  Make the diagonal elements small to cause gradual
		//        overflow when dividing by t(j,j).  To control the amount of
		//        scaling needed, the matrix is bidiagonal.
		//
		(*texp) = (*one) / DBLE((*(kd))+1)
		(*tscal) = math.Pow((*smlnum), (*texp))
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = MAX(1, (*(kd))+2-(*j)); (*i) <= (*(kd)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label240:
				}
				if (*j) > 1 && (*(kd)) > 0 {
					(*(ab))[(*(kd))-1][(*j)-1] = -(*one)
				}
				(*(ab))[(*(kd))+0][(*j)-1] = (*tscal)
				//Label250:
			}
			(*(b))[(*(n))-1] = (*one)
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 3; (*i) <= (Min((*(n))-(*j)+1, (*(kd))+1)); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = (*zero)
					//Label260:
				}
				if (*j) < (*(n)) && (*(kd)) > 0 {
					(*(ab))[1][(*j)-1] = -(*one)
				}
				(*(ab))[0][(*j)-1] = (*tscal)
				//Label270:
			}
			(*(b))[0] = (*one)
		}
		//
	} else if (*(imat)) == 15 {
		//
		//        Type 15:  one zero diagonal element.
		//
		(*iy) = (*(n))/2 + 1
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j), (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]))
				if (*j) != (*iy) {
					(*(ab))[(*(kd))+0][(*j)-1] = (*SIGN(two, &((*(ab))[(*(kd))+0][(*j)-1])))
				} else {
					(*(ab))[(*(kd))+0][(*j)-1] = (*zero)
				}
				//Label280:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j)+1, (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[0][(*j)-1]))
				if (*j) != (*iy) {
					(*(ab))[0][(*j)-1] = (*SIGN(two, &((*(ab))[0][(*j)-1])))
				} else {
					(*(ab))[0][(*j)-1] = (*zero)
				}
				//Label290:
			}
		}
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		Dscal((n), two, (b), func() *int {y := 1; return &y }())
		//
	} else if (*(imat)) == 16 {
		//
		//        Type 16:  Make the offdiagonal elements large to cause overflow
		//        when adding a column of T.  In the non-transposed case, the
		//        matrix is _constructed to cause overflow when adding a column in
		//        every other step.
		//
		(*tscal) = (*unfl) / (*ulp)
		(*tscal) = ((*one) - (*ulp)) / (*tscal)
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*(kd))+1; (*i)++ {
				(*(ab))[(*i)-1][(*j)-1] = (*zero)
				//Label300:
			}
			//Label310:
		}
		(*texp) = (*one)
		if (*(kd)) > 0 {
			if *upper {
				for (*j) = (*(n)); (*j) <= 1; (*j) += -(*(kd)) {
					for (*i) = (*j); (*i) <= (MAX(1, (*j)-(*(kd))+1)); (*i) += -2 {
						(*(ab))[1+((*j)-(*i))-1][(*i)-1] = -(*tscal) / DBLE((*(kd))+2)
						(*(ab))[(*(kd))+0][(*i)-1] = (*one)
						(*(b))[(*i)-1] = (*texp) * ((*one) - (*ulp))
						if (*i) > (MAX(1, (*j)-(*(kd))+1)) {
							(*(ab))[2+((*j)-(*i))-1][(*i)-0] = -((*tscal) / DBLE((*(kd))+2)) / DBLE((*(kd))+3)
							(*(ab))[(*(kd))+0][(*i)-0] = (*one)
							(*(b))[(*i)-0] = (*texp) * DBLE(((*(kd))+1)*((*(kd))+1)+(*(kd)))
						}
						(*texp) = (*texp) * (*two)
						//Label320:
					}
					(*(b))[MAX(1, (*j)-(*(kd))+1)-1] = (DBLE((*(kd))+2) / DBLE((*(kd))+3)) * (*tscal)
					//Label330:
				}
			} else {
				for (*j) = 1; (*j) <= (*(n)); (*j) += (*(kd)) {
					(*texp) = (*one)
					(*lenj) = (Min((*(kd))+1, (*(n))-(*j)+1))
					for (*i) = (*j); (*i) <= (Min((*(n)), (*j)+(*(kd))-1)); (*i) += 2 {
						(*(ab))[(*lenj)-((*i)-(*j))-1][(*j)-1] = -(*tscal) / DBLE((*(kd))+2)
						(*(ab))[0][(*j)-1] = (*one)
						(*(b))[(*j)-1] = (*texp) * ((*one) - (*ulp))
						if (*i) < (Min((*(n)), (*j)+(*(kd))-1)) {
							(*(ab))[(*lenj)-((*i)-(*j)+1)-1][(*i)+0] = -((*tscal) / DBLE((*(kd))+2)) / DBLE((*(kd))+3)
							(*(ab))[0][(*i)+0] = (*one)
							(*(b))[(*i)+0] = (*texp) * DBLE(((*(kd))+1)*((*(kd))+1)+(*(kd)))
						}
						(*texp) = (*texp) * (*two)
						//Label340:
					}
					(*(b))[Min((*(n)), (*j)+(*(kd))-1)-1] = (DBLE((*(kd))+2) / DBLE((*(kd))+3)) * (*tscal)
					//Label350:
				}
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(ab))[0][(*j)-1] = (*one)
				(*(b))[(*j)-1] = (DBLE((*j)))
				//Label360:
			}
		}
		//
	} else if (*(imat)) == 17 {
		//
		//        Type 17:  Generate a unit triangular matrix with elements
		//        between -1 and 1, and make the right hand side large so that it
		//        requires scaling.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j)-1, (*(kd))))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+1-(*lenj)-1][(*j)-1]))
				(*(ab))[(*(kd))+0][(*j)-1] = (DBLE((*j)))
				//Label370:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j), (*(kd))))
				if (*lenj) > 0 {
					Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[1][(*j)-1]))
				}
				(*(ab))[0][(*j)-1] = (DBLE((*j)))
				//Label380:
			}
		}
		//
		//        Set the right hand side so that the largest value is bignum.
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		(*iy) = (*Idamax((n), (b), func() *int {y := 1; return &y }()))
		(*bnorm) = (ABS(((*(b))[(*iy)-1])))
		(*bscal) = (*bignum) / MAX((*one), (*bnorm))
		Dscal((n), bscal, (b), func() *int {y := 1; return &y }())
		//
	} else if (*(imat)) == 18 {
		//
		//        Type 18:  Generate a triangular matrix with elements between
		//        bignum/kd and bignum so that at least one of the column
		//        norms will exceed bignum.
		//
		(*tleft) = (*bignum) / MAX((*one), DBLE((*(kd))))
		(*tscal) = (*bignum) * (DBLE((*(kd))) / DBLE((*(kd))+1))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*j), (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[(*(kd))+2-(*lenj)-1][(*j)-1]))
				for (*i) = (*(kd)) + 2 - (*lenj); (*i) <= (*(kd))+1; (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = SIGN(tleft, &((*(ab))[(*i)-1][(*j)-1])) + (*tscal)*(*(ab))[(*i)-1][(*j)-1]
					//Label390:
				}
				//Label400:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*lenj) = (Min((*(n))-(*j)+1, (*(kd))+1))
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), lenj, &((*(ab))[0][(*j)-1]))
				for (*i) = 1; (*i) <= (*lenj); (*i)++ {
					(*(ab))[(*i)-1][(*j)-1] = SIGN(tleft, &((*(ab))[(*i)-1][(*j)-1])) + (*tscal)*(*(ab))[(*i)-1][(*j)-1]
					//Label410:
				}
				//Label420:
			}
		}
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		Dscal((n), two, (b), func() *int {y := 1; return &y }())
	}
	//
	//     Flip the matrix if the transpose will be used.
	//
	if !blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		if *upper {
			for (*j) = 1; (*j) <= (*(n))/2; (*j)++ {
				(*lenj) = (Min((*(n))-2*(*j)+1, (*(kd))+1))
				Dswap(lenj, &((*(ab))[(*(kd))+0][(*j)-1]), (*(ldab))-1, &((*(ab))[(*(kd))+2-(*lenj)-1][(*(n))-(*j)+0]), -1)
				//Label430:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n))/2; (*j)++ {
				(*lenj) = (Min((*(n))-2*(*j)+1, (*(kd))+1))
				Dswap(lenj, &((*(ab))[0][(*j)-1]), func() *int {y := 1; return &y }(), &((*(ab))[(*lenj)-1][(*(n))-(*j)+2-(*lenj)-1]), -(*(ldab))+1)
				//Label440:
			}
		}
	}
	//
	return
	//
	//     End of Dlattb
	//
}
