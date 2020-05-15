package goblas

import (
	"math"

	"github.com/whipstein/golapack/blas"
)

// Dlattp generates a triangular test matrix in packed storage.
// imat and uplo uniquely specify the properties of the test
// matrix, which is returned in the array AP.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlattp( imat, uplo, trans, diag, iseed, n, a, b, work,
//                          info)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       intEGER            imat, info, N
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a(*), B(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlattp generates a triangular test matrix in packed storage.
// imat and uplo uniquely specify the properties of the test
// matrix, which is returned in the array AP.
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
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The upper or lower triangular matrix a, packed columnwise in
//          a linear array.  The j-th column of A is stored in the array
//          AP as follows:
//          if uplo = 'U', AP((j-1)*j/2 + i) = a(i,j) for 1<=i<=j;
//          if uplo = 'L',
//             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = a(i,j) for j<=i<=n.
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
func Dlattp(imat *int, uplo *byte, trans *byte, diag *byte, iseed *[]int, n *int, a *[]float64, b *[]float64, work *[]float64, info *int) {
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
	iy := new(int)
	j := new(int)
	jc := new(int)
	jcnext := new(int)
	jcount := new(int)
	jj := new(int)
	jl := new(int)
	jr := new(int)
	jx := new(int)
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
	stemp := new(float64)
	t := new(float64)
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
	(*path)[1] = *func() *[]byte {y := []byte("TP"); return &y }()
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
		(*packit) = 'C'
	} else {
		Dlatb4(path, -(*(imat)), (n), (n), _type, kl, ku, anorm, mode, cndnum, dist)
		(*packit) = 'R'
	}
	//
	//     imat <= 6:  Non-unit triangular matrix
	//
	if (*(imat)) <= 6 {
		Dlatms((n), (n), dist, (iseed), _type, (b), mode, cndnum, anorm, kl, ku, packit, (a), (n), (work), (info))
		//
		//     imat > 6:  Unit triangular matrix
		//     The diagonal is deliberately set to something other than 1.
		//
		//     imat = 7:  Matrix is the identity
		//
	} else if (*(imat)) == 7 {
		if *upper {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*jc)+(*i)-0] = (*zero)
					//Label10:
				}
				(*(a))[(*jc)+(*j)-0] = (*j)
				(*jc) = (*jc) + (*j)
				//Label20:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*jc)-1] = (*j)
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*jc)+(*i)-(*j)-1] = (*zero)
					//Label30:
				}
				(*jc) = (*jc) + (*(n)) - (*j) + 1
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
			(*jc) = 0
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*jc)+(*i)-1] = (*zero)
					//Label50:
				}
				(*(a))[(*jc)+(*j)-1] = (*j)
				(*jc) = (*jc) + (*j)
				//Label60:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*jc)-1] = (*j)
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*jc)+(*i)-(*j)-1] = (*zero)
					//Label70:
				}
				(*jc) = (*jc) + (*(n)) - (*j) + 1
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
		(*x) = SQRt((*cndnum)) - (*one)/SQRt((*cndnum))
		if (*(n)) > 2 {
			(*y) = SQRt((*two)/DBLE((*(n))-2)) * (*x)
		} else {
			(*y) = (*zero)
		}
		(*z) = (*x) * (*x)
		//
		if *upper {
			//
			//           Set the upper triangle of A with a unit triangular matrix
			//           of known condition number.
			//
			(*jc) = 1
			for (*j) = 2; (*j) <= (*(n)); (*j)++ {
				(*(a))[(*jc)+0] = (*y)
				if (*j) > 2 {
					(*(a))[(*jc)+(*j)-0] = (*(work))[(*j)-1]
				}
				if (*j) > 3 {
					(*(a))[(*jc)+(*j)-1] = (*(work))[(*(n))+(*j)-2]
				}
				(*jc) = (*jc) + (*j)
				//Label100:
			}
			(*jc) = (*jc) - (*(n))
			(*(a))[(*jc)+0] = (*z)
			for (*j) = 2; (*j) <= (*(n))-1; (*j)++ {
				(*(a))[(*jc)+(*j)-1] = (*y)
				//Label110:
			}
		} else {
			//
			//           Set the lower triangle of A with a unit triangular matrix
			//           of known condition number.
			//
			for (*i) = 2; (*i) <= (*(n))-1; (*i)++ {
				(*(a))[(*i)-1] = (*y)
				//Label120:
			}
			(*(a))[(*(n))-1] = (*z)
			(*jc) = (*(n)) + 1
			for (*j) = 2; (*j) <= (*(n))-1; (*j)++ {
				(*(a))[(*jc)+0] = (*(work))[(*j)-0]
				if (*j) < (*(n))-1 {
					(*(a))[(*jc)+1] = (*(work))[(*(n))+(*j)-0]
				}
				(*(a))[(*jc)+(*(n))-(*j)-1] = (*y)
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label130:
			}
		}
		//
		//        Fill in the zeros using Givens rotations
		//
		if *upper {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
				(*jcnext) = (*jc) + (*j)
				(*ra) = (*(a))[(*jcnext)+(*j)-0]
				(*rb) = (*two)
				Drotg(ra, rb, c, s)
				//
				//              Multiply by[c  s; -s  c] on the left.
				//
				if (*(n)) > (*j)+1 {
					(*jx) = (*jcnext) + (*j)
					for (*i) = (*j) + 2; (*i) <= (*(n)); (*i)++ {
						(*stemp) = (*c)*(*(a))[(*jx)+(*j)-1] + (*s)*(*(a))[(*jx)+(*j)+0]
						(*(a))[(*jx)+(*j)+0] = -(*s)*(*(a))[(*jx)+(*j)-1] + (*c)*(*(a))[(*jx)+(*j)+0]
						(*(a))[(*jx)+(*j)-1] = (*stemp)
						(*jx) = (*jx) + (*i)
						//Label140:
					}
				}
				//
				//              Multiply by[-c -s;  s -c] on the right.
				//
				if (*j) > 1 {
					Drot((*j)-1, &((*(a))[(*jcnext)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*jc)-1]), func() *int {y := 1; return &y }(), -(*c), -(*s))
				}
				//
				//              Negate a(J,J+1).
				//
				(*(a))[(*jcnext)+(*j)-0] = -(*(a))[(*jcnext)+(*j)-0]
				(*jc) = (*jcnext)
				//Label150:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n))-1; (*j)++ {
				(*jcnext) = (*jc) + (*(n)) - (*j) + 1
				(*ra) = (*(a))[(*jc)+0]
				(*rb) = (*two)
				Drotg(ra, rb, c, s)
				//
				//              Multiply by[c -s;  s  c] on the right.
				//
				if (*(n)) > (*j)+1 {
					Drot((*(n))-(*j)-1, &((*(a))[(*jcnext)+0]), func() *int {y := 1; return &y }(), &((*(a))[(*jc)+1]), func() *int {y := 1; return &y }(), c, -(*s))
				}
				//
				//              Multiply by[-c  s; -s -c] on the left.
				//
				if (*j) > 1 {
					(*jx) = 1
					for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
						(*stemp) = -(*c)*(*(a))[(*jx)+(*j)-(*i)-1] + (*s)*(*(a))[(*jx)+(*j)-(*i)+0]
						(*(a))[(*jx)+(*j)-(*i)+0] = -(*s)*(*(a))[(*jx)+(*j)-(*i)-1] - (*c)*(*(a))[(*jx)+(*j)-(*i)+0]
						(*(a))[(*jx)+(*j)-(*i)-1] = (*stemp)
						(*jx) = (*jx) + (*(n)) - (*i) + 1
						//Label160:
					}
				}
				//
				//              Negate a(J+1,j).
				//
				(*(a))[(*jc)+0] = -(*(a))[(*jc)+0]
				(*jc) = (*jcnext)
				//Label170:
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[(*jc)-1]))
				(*(a))[(*jc)+(*j)-0] = (*SIGN(two, &((*(a))[(*jc)+(*j)-0])))
				(*jc) = (*jc) + (*j)
				//Label180:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*jc)-1]))
				(*(a))[(*jc)-1] = (*SIGN(two, &((*(a))[(*jc)-1])))
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label190:
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*j)-1, &((*(a))[(*jc)-1]))
				Dscal((*j)-1, tscal, &((*(a))[(*jc)-1]), func() *int {y := 1; return &y }())
				(*(a))[(*jc)+(*j)-0] = (*SIGN(one, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*jc) = (*jc) + (*j)
				//Label200:
			}
			(*(a))[(*(n))*((*(n))+1)/1] = (*smlnum)
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j), &((*(a))[(*jc)+0]))
				Dscal((*(n))-(*j), tscal, &((*(a))[(*jc)+0]), func() *int {y := 1; return &y }())
				(*(a))[(*jc)-1] = (*SIGN(one, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label210:
			}
			(*(a))[0] = (*smlnum)
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*j)-1, &((*(a))[(*jc)-1]))
				(*(a))[(*jc)+(*j)-0] = (*SIGN(one, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*jc) = (*jc) + (*j)
				//Label220:
			}
			(*(a))[(*(n))*((*(n))+1)/1] = (*smlnum)
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j), &((*(a))[(*jc)+0]))
				(*(a))[(*jc)-1] = (*SIGN(one, Dlarnd(func() *int {y := 2; return &y }(), (iseed))))
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label230:
			}
			(*(a))[0] = (*smlnum)
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
			(*jc) = ((*(n))-1)*(*(n))/2 + 1
			for (*j) = (*(n)); (*j) <= 1; (*j) += -1 {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*jc)+(*i)-0] = (*zero)
					//Label240:
				}
				if (*jcount) <= 2 {
					(*(a))[(*jc)+(*j)-0] = (*smlnum)
				} else {
					(*(a))[(*jc)+(*j)-0] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				(*jc) = (*jc) - (*j) + 1
				//Label250:
			}
		} else {
			(*jcount) = 1
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) + 1; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*jc)+(*i)-(*j)-1] = (*zero)
					//Label260:
				}
				if (*jcount) <= 2 {
					(*(a))[(*jc)-1] = (*smlnum)
				} else {
					(*(a))[(*jc)-1] = (*one)
				}
				(*jcount) = (*jcount) + 1
				if (*jcount) > 4 {
					(*jcount) = 1
				}
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label270:
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
				//Label280:
			}
		} else {
			(*(b))[(*(n))-1] = (*zero)
			for (*i) = 1; (*i) <= (*(n))-1; (*i) += 2 {
				(*(b))[(*i)-1] = (*zero)
				(*(b))[(*i)+0] = (*smlnum)
				//Label290:
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-2; (*i)++ {
					(*(a))[(*jc)+(*i)-0] = (*zero)
					//Label300:
				}
				if (*j) > 1 {
					(*(a))[(*jc)+(*j)-1] = -(*one)
				}
				(*(a))[(*jc)+(*j)-0] = (*tscal)
				(*jc) = (*jc) + (*j)
				//Label310:
			}
			(*(b))[(*(n))-1] = (*one)
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) + 2; (*i) <= (*(n)); (*i)++ {
					(*(a))[(*jc)+(*i)-(*j)-1] = (*zero)
					//Label320:
				}
				if (*j) < (*(n)) {
					(*(a))[(*jc)+0] = -(*one)
				}
				(*(a))[(*jc)-1] = (*tscal)
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label330:
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[(*jc)-1]))
				if (*j) != (*iy) {
					(*(a))[(*jc)+(*j)-0] = (*SIGN(two, &((*(a))[(*jc)+(*j)-0])))
				} else {
					(*(a))[(*jc)+(*j)-0] = (*zero)
				}
				(*jc) = (*jc) + (*j)
				//Label340:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*jc)-1]))
				if (*j) != (*iy) {
					(*(a))[(*jc)-1] = (*SIGN(two, &((*(a))[(*jc)-1])))
				} else {
					(*(a))[(*jc)-1] = (*zero)
				}
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label350:
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
		for (*j) = 1; (*j) <= (*(n))*((*(n))+1)/2; (*j)++ {
			(*(a))[(*j)-1] = (*zero)
			//Label360:
		}
		(*texp) = (*one)
		if *upper {
			(*jc) = ((*(n))-1)*(*(n))/2 + 1
			for (*j) = (*(n)); (*j) <= 2; (*j) += -2 {
				(*(a))[(*jc)-1] = -(*tscal) / DBLE((*(n))+1)
				(*(a))[(*jc)+(*j)-0] = (*one)
				(*(b))[(*j)-1] = (*texp) * ((*one) - (*ulp))
				(*jc) = (*jc) - (*j) + 1
				(*(a))[(*jc)-1] = -((*tscal) / DBLE((*(n))+1)) / DBLE((*(n))+2)
				(*(a))[(*jc)+(*j)-1] = (*one)
				(*(b))[(*j)-0] = (*texp) * DBLE((*(n))*(*(n))+(*(n))-1)
				(*texp) = (*texp) * (*two)
				(*jc) = (*jc) - (*j) + 2
				//Label370:
			}
			(*(b))[0] = (DBLE((*(n))+1) / DBLE((*(n))+2)) * (*tscal)
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n))-1; (*j) += 2 {
				(*(a))[(*jc)+(*(n))-(*j)-1] = -(*tscal) / DBLE((*(n))+1)
				(*(a))[(*jc)-1] = (*one)
				(*(b))[(*j)-1] = (*texp) * ((*one) - (*ulp))
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				(*(a))[(*jc)+(*(n))-(*j)-0] = -((*tscal) / DBLE((*(n))+1)) / DBLE((*(n))+2)
				(*(a))[(*jc)-1] = (*one)
				(*(b))[(*j)+0] = (*texp) * DBLE((*(n))*(*(n))+(*(n))-1)
				(*texp) = (*texp) * (*two)
				(*jc) = (*jc) + (*(n)) - (*j)
				//Label380:
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
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*j)-1, &((*(a))[(*jc)-1]))
				(*(a))[(*jc)+(*j)-0] = (*zero)
				(*jc) = (*jc) + (*j)
				//Label390:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				if (*j) < (*(n)) {
					Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j), &((*(a))[(*jc)+0]))
				}
				(*(a))[(*jc)-1] = (*zero)
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label400:
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
		//
		(*tleft) = (*bignum) / MAX((*one), DBLE((*(n))-1))
		(*tscal) = (*bignum) * (DBLE((*(n))-1) / MAX((*one), DBLE((*(n)))))
		if *upper {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), j, &((*(a))[(*jc)-1]))
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*(a))[(*jc)+(*i)-0] = SIGN(tleft, &((*(a))[(*jc)+(*i)-0])) + (*tscal)*(*(a))[(*jc)+(*i)-0]
					//Label410:
				}
				(*jc) = (*jc) + (*j)
				//Label420:
			}
		} else {
			(*jc) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y }(), (iseed), (*(n))-(*j)+1, &((*(a))[(*jc)-1]))
				for (*i) = (*j); (*i) <= (*(n)); (*i)++ {
					(*(a))[(*jc)+(*i)-(*j)-1] = SIGN(tleft, &((*(a))[(*jc)+(*i)-(*j)-1])) + (*tscal)*(*(a))[(*jc)+(*i)-(*j)-1]
					//Label430:
				}
				(*jc) = (*jc) + (*(n)) - (*j) + 1
				//Label440:
			}
		}
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), (b))
		Dscal((n), two, (b), func() *int {y := 1; return &y }())
	}
	//
	//     Flip the matrix across its counter-diagonal if the transpose will
	//     be used.
	//
	if !blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		if *upper {
			(*jj) = 1
			(*jr) = (*(n)) * ((*(n)) + 1) / 2
			for (*j) = 1; (*j) <= (*(n))/2; (*j)++ {
				(*jl) = (*jj)
				for (*i) = (*j); (*i) <= (*(n))-(*j); (*i)++ {
					(*t) = (*(a))[(*jr)-(*i)+(*j)-1]
					(*(a))[(*jr)-(*i)+(*j)-1] = (*(a))[(*jl)-1]
					(*(a))[(*jl)-1] = (*t)
					(*jl) = (*jl) + (*i)
					//Label450:
				}
				(*jj) = (*jj) + (*j) + 1
				(*jr) = (*jr) - ((*(n)) - (*j) + 1)
				//Label460:
			}
		} else {
			(*jl) = 1
			(*jj) = (*(n)) * ((*(n)) + 1) / 2
			for (*j) = 1; (*j) <= (*(n))/2; (*j)++ {
				(*jr) = (*jj)
				for (*i) = (*j); (*i) <= (*(n))-(*j); (*i)++ {
					(*t) = (*(a))[(*jl)+(*i)-(*j)-1]
					(*(a))[(*jl)+(*i)-(*j)-1] = (*(a))[(*jr)-1]
					(*(a))[(*jr)-1] = (*t)
					(*jr) = (*jr) - (*i)
					//Label470:
				}
				(*jl) = (*jl) + (*(n)) - (*j) + 1
				(*jj) = (*jj) - (*j) - 1
				//Label480:
			}
		}
	}
	//
	return
	//
	//     End of Dlattp
	//
}
