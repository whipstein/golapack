package goblas

import (
	"math"

	"github.com/whipstein/golapack/blas"
)

// Dlattr generates a triangular test matrix.
// imat and uplo uniquely specify the properties of the test
// matrix, which is returned in the array A.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlattr( imat, uplo, trans, diag, iseed, n, a, lda, b,
//                          work, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       intEGER            imat, info, lda, N
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), B(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlattr generates a triangular test matrix.
// imat and uplo uniquely specify the properties of the test
// matrix, which is returned in the array A.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] imat
// \verbatim
//          imat is intEGER
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
//          = 'C':  Conjugate transpose (= Transpose)
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
//          iseed is intEGER array, dimension (4)
//          The seed vector for the random number generator (used in
//          Dlatms).  Modified on exit.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The order of the matrix to be generated.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The triangular matrix A.  If uplo = 'U', the leading n by n
//          upper triangular part of the array A contains the upper
//          triangular matrix, and the strictly lower triangular part of
//          A is not referenced.  If uplo = 'L', the leading n by n lower
//          triangular part of the array A contains the lower triangular
//          matrix, and the strictly upper triangular part of A is not
//          referenced.  If diag = 'U', the diagonal elements of A are
//          set so that a(k,k) = k for 1 <= k <= n.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.  lda >= max(1,N).
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (n)
//          The right hand side vector, if imat > 10.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (3*n)
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
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
func Dlattr(imat *int, uplo *byte, trans *byte, diag *byte, iseed *[]int, n *int, a *[][]float64, lda *int, b *[]float64, work *[]float64, info *int) {
	one := new(float64)
	two := new(float64)
	zero := new(float64)
	upper := new(bool)
	dist := new(byte)
	_type := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	iy := new(int)
	j := new(int)
	jcount := new(int)
	kl := new(int)
	ku := new(int)
	mode := new(int)
	anorm := new(float64)
	bignum := new(float64)
	bnorm := new(float64)
	bscal := new(float64)
	c := new(float64)
	cndnum := new(float64)
	plus1 := new(float64)
	plus2 := new(float64)
	ra := new(float64)
	rb := new(float64)
	rexp := new(float64)
	s := new(float64)
	sfac := new(float64)
	smlnum := new(float64)
	star1 := new(float64)
	texp := new(float64)
	tleft := new(float64)
	tscal := new(float64)
	ulp := new(float64)
	unfl := new(float64)
	x := new(float64)
	y := new(float64)
	z := new(float64)
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
	(*path)[1] = *func() *[]byte {y := []byte("TR"); return &y }()
	(*unfl) = (*Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()))
	(*ulp) = Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()) * Dlamch(func() *[]byte {y := []byte("Base"); return &y }())
	(*smlnum) = (*unfl)
	(*bignum) = ((*one) - (*ulp)) / (*smlnum)
	Dlabad(smlnum, bignum)
	if ((*(imat)) >= 7 && (*(imat)) <= 10) || (*(imat)) == 18 {
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
	} else {
		Dlatb4(path, -(*(imat)), (n), (n), _type, kl, ku, anorm, mode, cndnum, dist)
	}
	//
	//     imat <= 6:  Non-unit triangular matrix
	//
	if (*(imat)) <= 6 {
		Dlatms((n), (n), dist, (iseed), _type, (b), mode, cndnum, anorm, kl, ku, func() *[]byte {y := []byte("No packing"); return &y }(), (a), (lda), (work), (info))
		//
		//     imat > 6:  Unit triangular matrix
		//     The diagonal is deliberately set to something other than 1.
		//
		//     imat = 7:  Matrix is the identity
		//
	} else if (*(imat)) == 7 {
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label10:
				}
				(*(a))[(*j)-1][(*j)-1] = (*j)
				//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*j)-1][(*j)-1] = (*j)
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label30:
				}
				//Label40:
			}
		}
		//
		//     imat > 7:  Non-trivial unit triangular matrix
		//
		//     Generate a unit triangular matrix T with condition cndnum by
		//     forming a triangular matrix with known singular values and
		//     filling in the zero entries with Givens rotations.
		//
	} else if (*(imat)) <= 10 {
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label50:
				}
				(*(a))[(*j)-1][(*j)-1] = (*j)
				//Label60:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*j)-1][(*j)-1] = (*j)
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label70:
				}
				//Label80:
			}
		}
		//
		//        Since the trace of a unit triangular matrix is 1, the product
		//        of its singular values must be 1.  Let s = sqrt(cndnum),
		//        x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
		//        The following triangular matrix has singular values s, 1, 1,
		//        ..., 1, 1/s:
		//
		//        1  y  y  y  ...  y  y  z
		//           1  0  0  ...  0  0  y
		//              1  0  ...  0  0  y
		//                 .  ...  .  .  .
		//                     .   .  .  .
		//                         1  0  y
		//                            1  y
		//                               1
		//
		//        To fill in the zeros, we first multiply by a matrix with small
		//        condition number of the form
		//
		//        1  0  0  0  0  ...
		//           1  +  *  0  0  ...
		//              1  +  0  0  0
		//                 1  +  *  0  0
		//                    1  +  0  0
		//                       ...
		//                          1  +  0
		//                             1  0
		//                                1
		//
		//        Each element marked with a '*' is formed by taking the product
		//        of the adjacent elements marked with '+'.  The '*'s can be
		//        chosen freely, and the '+'s are chosen so that the inverse of
		//        T will have elements of the same magnitude as T.  If the *'s in
		//        both T and inv (t) have small magnitude, T is well conditioned.
		//        The two offdiagonals of T are stored in work.
		//
		//        The product of these two matrices has the form
		//
		//        1  y  y  y  y  y  .  y  y  z
		//           1  +  *  0  0  .  0  0  y
		//              1  +  0  0  .  0  0  y
		//                 1  +  *  .  .  .  .
		//                    1  +  .  .  .  .
		//                       .  .  .  .  .
		//                          .  .  .  .
		//                             1  +  y
		//                                1  y
		//                                   1
		//
		//        Now we multiply by Givens rotations, using the fact that
		//
		//             [c   s][1   w][-c  -s] = [1  -w]
		//             [-s   c][0   1][s  -c]   [0   1]
		//        and
		//             [-c  -s][1   0][c   s] = [1   0]
		//             [s  -c][w   1][-s   c]   [-w   1]
		//
		//        where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).
		//
		(*star1) = 0.25
		(*sfac) = 0.5
		(*plus1) = (*sfac)
		for (*j) = 1; (*j) <= (*(n)); (*j) += 2 {
			(*plus2) = (*star1) / (*plus1)
			(*(work))[(*j)-1] = (*plus1)
			(*(work))[(*(n))+(*j)-1] = (*star1)
			if (*j)+1 <= (*(n)) {
				(*(work))[(*j)+0] = (*plus2)
				(*(work))[(*(n))+(*j)+0] = (*zero)
				(*plus1) = (*star1) / (*plus2)
				(*rexp) = (*Dlarnd(func() *int {y := 2; return &y }(), (iseed)))
				(*star1) = (*star1) * (math.Pow((*sfac), (*rexp)))
				if (*rexp) < (*zero) {
					(*star1) = -math.Pow((*sfac), ((*one) - (*rexp)))
				} else {
					(*star1) = math.Pow((*sfac), ((*one) + (*rexp)))
				}
			}
			//Label90:
		}
		//
		(*x) = SQRt((*cndnum)) - 1/SQRt((*cndnum))
		if (*(n)) > 2 {
			(*y) = SQRt(2./((*(n))-2)) * (*x)
		} else {
			(*y) = (*zero)
		}
		(*z) = (*x) * (*x)
		//
		if *upper {
			if (*(n)) > 3 {
				Dcopy((*(n))-3, (work), func() *int {y := 1; return &y }(), &((*(a))[1][2]), (*(lda))+1)
				if (*(n)) > 4 {
					Dcopy((*(n))-4, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(a))[1][3]), (*(lda))+1)
				}
			}
			for (*j) = 2; (*j) <= (*(n))-1; (*j)++ {
				(*(a))[0][(*j)-1] = (*y)
				(*(a))[(*j)-1][(*(n))-1] = (*y)
				//Label100:
			}
			(*(a))[0][(*(n))-1] = (*z)
		} else {
			if (*(n)) > 3 {
				Dcopy((*(n))-3, (work), func() *int {y := 1; return &y }(), &((*(a))[2][1]), (*(lda))+1)
				if (*(n)) > 4 {
					Dcopy((*(n))-4, &((*(work))[(*(n))+0]), func() *int {y := 1; return &y }(), &((*(a))[3][1]), (*(lda))+1)
				}
			}
			for (*j) = 2; (*j) <= (*(n))-1; (*j)++ {
				(*(a))[(*j)-1][0] = (*y)
				(*(a))[(*(n))-1][(*j)-1] = (*y)
				//Label110:
			}
			(*(a))[(*(n))-1][0] = (*z)
		}
		//
		//        Fill in the zeros using Givens rotations.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
				(*ra) = (*(a))[(*j)-1][(*j)+0]
				(*rb) = 2.0
				Drotg(ra, rb, c, s)
				//
				//              Multiply by[c  s; -s  c] on the left.
				//
				if (*(n)) > (*j)+1 {
					Drot((*(n))-(*j)-1, &((*(a))[(*j)-1][(*j)+1]), (lda), &((*(a))[(*j)+0][(*j)+1]), (lda), c, s)
				}
				//
				//              Multiply by[-c -s;  s -c] on the right.
				//
				if (*j) > 1 {
					Drot((*j)-1, &((*(a))[0][(*j)+0]), func() *int {y := 1; return &y }(), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }(), -(*c), -(*s))
				}
				//
				//              Negate a(J,J+1).
				//
				(*(a))[(*j)-1][(*j)+0] = -(*(a))[(*j)-1][(*j)+0]
				//Label120:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
				(*ra) = (*(a))[(*j)+0][(*j)-1]
				(*rb) = 2.0
				Drotg(ra, rb, c, s)
				//
				//              Multiply by[c -s;  s  c] on the right.
				//
				if (*(n)) > (*j)+1 {
					Drot((*(n))-(*j)-1, &((*(a))[(*j)+1][(*j)+0]), func() *int {y := 1; return &y }(), &((*(a))[(*j)+1][(*j)-1]), func() *int {y := 1; return &y }(), c, -(*s))
				}
				//
				//              Multiply by[-c  s; -s -c] on the left.
				//
				if (*j) > 1 {
					Drot((*j)-1, &((*(a))[(*j)-1][0]), (lda), &((*(a))[(*j)+0][0]), (lda), -(*c), s)
				}
				//
				//              Negate a(J+1,j).
				//
				(*(a))[(*j)+0][(*j)-1] = -(*(a))[(*j)+0][(*j)-1]
				//Label130:
			}
		}
		//
		//     imat > 10:  Pathological test cases.  These triangular matrices
		//     are badly scaled or badly conditioned, so when used in solving a
		//     triangular system they may cause overflow in the solution vector.
		//
	} else if (*(imat)) == 11 {
		//
		//        Type 11:  Generate a triangular matrix with elements between
		//        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
		//        Make the right hand side large so that it requires scaling.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[0][(*j)-1]))
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(two, &((*(a))[(*j)-1][(*j)-1])))
				//Label140:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*j)-1][(*j)-1]))
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(two, &((*(a))[(*j)-1][(*j)-1])))
				//Label150:
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
	} else if (*(imat)) == 12 {
		//
		//        Type 12:  Make the first diagonal element in the solve small to
		//        cause immediate overflow when dividing by t(j,j).
		//        In type 12, the offdiagonal elements are small (cnorm(j) < 1).
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		(*tscal) = (*one) / MAX((*one), DBLE((*(n))-1))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[0][(*j)-1]))
				Dscal((*j)-1, tscal, &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(one, &((*(a))[(*j)-1][(*j)-1])))
				//Label160:
			}
			(*(a))[(*(n))-1][(*(n))-1] = (*smlnum) * (*(a))[(*(n))-1][(*(n))-1]
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*j)-1][(*j)-1]))
				if (*(n)) > (*j) {
					Dscal((*(n))-(*j), tscal, &((*(a))[(*j)+0][(*j)-1]), func() *int {y := 1; return &y }())
				}
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(one, &((*(a))[(*j)-1][(*j)-1])))
				//Label170:
			}
			(*(a))[0][0] = (*smlnum) * (*(a))[0][0]
		}
		//
	} else if (*(imat)) == 13 {
		//
		//        Type 13:  Make the first diagonal element in the solve small to
		//        cause immediate overflow when dividing by t(j,j).
		//        In type 13, the offdiagonal elements are O1 (cnorm(j) > 1).
		//
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[0][(*j)-1]))
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(one, &((*(a))[(*j)-1][(*j)-1])))
				//Label180:
			}
			(*(a))[(*(n))-1][(*(n))-1] = (*smlnum) * (*(a))[(*(n))-1][(*(n))-1]
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*j)-1][(*j)-1]))
				(*(a))[(*j)-1][(*j)-1] = (*SIGN(one, &((*(a))[(*j)-1][(*j)-1])))
				//Label190:
			}
			(*(a))[0][0] = (*smlnum) * (*(a))[0][0]
		}
		//
	} else if (*(imat)) == 14 {
		//
		//        Type 14:  T is diagonal with small numbers on the diagonal to
		//        make the growth factor underflow, but a small right hand side
		//        chosen so that the solution does not overflow.
		//
		if *upper {
			(*jcount) = 1
			for (*j) = (*(n)); (*j) <= 1; (*j) += -1 {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label200:
				}
				if (*jcount) <= 2 {
					(*(a))[(*j)-1][(*j)-1] = (*smlnum)
				} else {
					(*(a))[(*j)-1][(*j)-1] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				//Label210:
			}
		} else {
			(*jcount) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label220:
				}
				if (*jcount) <= 2 {
					(*(a))[(*j)-1][(*j)-1] = (*smlnum)
				} else {
					(*(a))[(*j)-1][(*j)-1] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				//Label230:
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
				//Label240:
			}
		} else {
			(*(b))[(*(n))-1] = (*zero)
			for (*i) = 1; (*i) <= (*(n))-1; (*i) += 2 {
				(*(b))[(*i)-1] = (*zero)
				(*(b))[(*i)+0] = (*smlnum)
				//Label250:
			}
		}
		//
	} else if (*(imat)) == 15 {
		//
		//        Type 15:  Make the diagonal elements small to cause gradual
		//        overflow when dividing by t(j,j).  To control the amount of
		//        scaling needed, the matrix is bidiagonal.
		//
		(*texp) = (*one) / MAX((*one), DBLE((*(n))-1))
		(*tscal) = math.Pow((*smlnum), (*texp))
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-2; (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = 0.
					//Label260:
				}
				if (*j) > 1 {
					(*(a))[(*j)-0][(*j)-1] = -(*one)
				}
				(*(a))[(*j)-1][(*j)-1] = (*tscal)
				//Label270:
			}
			(*(b))[(*(n))-1] = (*one)
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) + 2; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = 0.
					//Label280:
				}
				if (*j) < (*(n)) {
					(*(a))[(*j)+0][(*j)-1] = -(*one)
				}
				(*(a))[(*j)-1][(*j)-1] = (*tscal)
				//Label290:
			}
			(*(b))[0] = (*one)
		}
		//
	} else if (*(imat)) == 16 {
		//
		//        Type 16:  one zero diagonal element.
		//
		(*iy) = (*(n))/2 + 1
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[0][(*j)-1]))
				if (*j) != (*iy) {
					(*(a))[(*j)-1][(*j)-1] = (*SIGN(two, &((*(a))[(*j)-1][(*j)-1])))
				} else {
					(*(a))[(*j)-1][(*j)-1] = (*zero)
				}
				//Label300:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*j)-1][(*j)-1]))
				if (*j) != (*iy) {
					(*(a))[(*j)-1][(*j)-1] = (*SIGN(two, &((*(a))[(*j)-1][(*j)-1])))
				} else {
					(*(a))[(*j)-1][(*j)-1] = (*zero)
				}
				//Label310:
			}
		}
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		Dscal((n), two, (b), func() *int {y := 1; return &y }())
		//
	} else if (*(imat)) == 17 {
		//
		//        Type 17:  Make the offdiagonal elements large to cause overflow
		//        when adding a column of T.  In the non-transposed case, the
		//        matrix is _constructed to cause overflow when adding a column in
		//        every other step.
		//
		(*tscal) = (*unfl) / (*ulp)
		(*tscal) = ((*one) - (*ulp)) / (*tscal)
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = 1; (*i) <= (*(n)); (*i)++ {
				(*(a))[(*i)-1][(*j)-1] = 0.
				//Label320:
			}
			//Label330:
		}
		(*texp) = (*one)
		if *upper {
			for (*j) = (*(n)); (*j) <= 2; (*j) += -2 {
				(*(a))[0][(*j)-1] = -(*tscal) / DBLE((*(n))+1)
				(*(a))[(*j)-1][(*j)-1] = (*one)
				(*(b))[(*j)-1] = (*texp) * ((*one) - (*ulp))
				(*(a))[0][(*j)-0] = -((*tscal) / DBLE((*(n))+1)) / DBLE((*(n))+2)
				(*(a))[(*j)-0][(*j)-0] = (*one)
				(*(b))[(*j)-0] = (*texp) * DBLE((*(n))*(*(n))+(*(n))-1)
				(*texp) = (*texp) * 2.
				//Label340:
			}
			(*(b))[0] = (DBLE((*(n))+1) / DBLE((*(n))+2)) * (*tscal)
		} else {
			for (*j) = 1; (*j) <= (*(n))-1; (*j) += 2 {
				(*(a))[(*(n))-1][(*j)-1] = -(*tscal) / DBLE((*(n))+1)
				(*(a))[(*j)-1][(*j)-1] = (*one)
				(*(b))[(*j)-1] = (*texp) * ((*one) - (*ulp))
				(*(a))[(*(n))-1][(*j)+0] = -((*tscal) / DBLE((*(n))+1)) / DBLE((*(n))+2)
				(*(a))[(*j)+0][(*j)+0] = (*one)
				(*(b))[(*j)+0] = (*texp) * DBLE((*(n))*(*(n))+(*(n))-1)
				(*texp) = (*texp) * 2.
				//Label350:
			}
			(*(b))[(*(n))-1] = (DBLE((*(n))+1) / DBLE((*(n))+2)) * (*tscal)
		}
		//
	} else if (*(imat)) == 18 {
		//
		//        Type 18:  Generate a unit triangular matrix with elements
		//        between -1 and 1, and make the right hand side large so that it
		//        requires scaling.
		//
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*j)-1, &((*(a))[0][(*j)-1]))
				(*(a))[(*j)-1][(*j)-1] = (*zero)
				//Label360:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				if (*j) < (*(n)) {
					Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j), &((*(a))[(*j)+0][(*j)-1]))
				}
				(*(a))[(*j)-1][(*j)-1] = (*zero)
				//Label370:
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
	} else if (*(imat)) == 19 {
		//
		//        Type 19:  Generate a triangular matrix with elements between
		//        bignum/(n-1) and bignum so that at least one of the column
		//        norms will exceed bignum.
		//        1/3/91:  Dlatrs no longer can handle this case
		//
		(*tleft) = (*bignum) / MAX((*one), DBLE((*(n))-1))
		(*tscal) = (*bignum) * (DBLE((*(n))-1) / MAX((*one), DBLE((*(n)))))
		if *upper {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[0][(*j)-1]))
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = SIGN(tleft, &((*(a))[(*i)-1][(*j)-1])) + (*tscal)*(*(a))[(*i)-1][(*j)-1]
					//Label380:
				}
				//Label390:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*j)-1][(*j)-1]))
				for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = SIGN(tleft, &((*(a))[(*i)-1][(*j)-1])) + (*tscal)*(*(a))[(*i)-1][(*j)-1]
					//Label400:
				}
				//Label410:
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
				Dswap((*(n))-2*(*j)+1, &((*(a))[(*j)-1][(*j)-1]), (lda), &((*(a))[(*j)+0][(*(n))-(*j)+0]), -1)
				//Label420:
			}
		} else {
			for (*j) = 1; (*j) <= (*(n))/2; (*j)++ {
				Dswap((*(n))-2*(*j)+1, &((*(a))[(*j)-1][(*j)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*(n))-(*j)+0][(*j)+0]), -(*(lda)))
				//Label430:
			}
		}
	}
	//
	return
	//
	//     End of Dlattr
	//
}
