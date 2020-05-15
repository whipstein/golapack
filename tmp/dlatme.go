package goblas

import (
	"github.com/whipstein/golapack/blas"
)

//    dlatme generates random non-symmetric square matrices with
//    specified eigenvalues for testing lapACK programs.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE dlatme( n, dist, iseed, d, mode, cond, dmax, EI,
//         rsign,
//                          upper, SIM, DS, modes, conds, kl, ku, anorm,
//         a,
//                          lda, work, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          dist, rsign, SIM, upper
//       inTEGER            info, kl, ku, lda, mode, modes, N
//       DOUBLE PRECISION   anorm, cond, conds, dmax
//       ..
//       .. Array Arguments ..
//       CHARACTER          EI(*)
//       inTEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), d(*), DS(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    dlatme generates random non-symmetric square matrices with
//    specified eigenvalues for testing lapACK programs.
//
//    dlatme operates by applying the following sequence of
//    operations:
//
//    1. Set the diagonal to d, where D may be input or
//         computed according to mode, cond, dmax, and rsign
//         as described below.
//
//    2. If complex conjugate pairs are desired (mode=0 and EI1='R',
//         or mode=5), certain pairs of adjacent elements of D are
//         interpreted as the real and complex parts of a complex
//         conjugate pair; A thus becomes block diagonal, with 1x1
//         and 2x2 blocks.
//
//    3. If upper='T', the upper triangle of A is set to random values
//         out of distribution dist.
//
//    4. If SIM='T', A is multiplied on the left by a random matrix
//         x, whose singular values are specified by DS, modes, and
//         conds, and on the right by X inverse.
//
//    5. If kl < N-1, the lower bandwidth is reduced to kl using
//         Householder transformations.  If ku < N-1, the upper
//         bandwidth is reduced to ku.
//
//    6. If anorm is not negative, the matrix is scaled to have
//         maximum-element-norm anorm.
//
//    (Note: since the matrix cannot be reduced beyond Hessenberg form,
//     no packing options are available.)
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] N
// \verbatim
//          N is inTEGER
//           The number of columns (or rows) of A. Not modified.
// \endverbatim
//
// \param[in] dist
// \verbatim
//          dist is CHARACTER*1
//           On entry, dist specifies the type of distribution to be used
//           to generate the random eigen-/singular values, and for the
//           upper triangle (see upper).
//           'U' => UNIFORM( 0, 1)  ( 'U' for uniform)
//           'S' => UNIFORM( -1, 1) ( 'S' for symmetric)
//           'N' => normaL( 0, 1)   ( 'N' for normal)
//           Not modified.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is inTEGER array, dimension ( 4)
//           On entry iseed specifies the seed of the random number
//           generator. They should lie between 0 and 4095 inclusive,
//           and iseed(4) should be odd. The random number generator
//           uses a linear congruential sequence limited to small
//           integers, and so should produce machine independent
//           random numbers. The values of iseed are changed on
//           exit, and can be used in the next call to dlatme
//           to continue the same random number sequence.
//           Changed on exit.
// \endverbatim
//
// \param[in,out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension ( N)
//           This array is used to specify the eigenvalues of A.  If
//           mode=0, then D is assumed to contain the eigenvalues (but
//           see the description of EI), otherwise they will be
//           computed according to mode, cond, dmax, and rsign and
//           placed in D.
//           Modified if mode is nonzero.
// \endverbatim
//
// \param[in] mode
// \verbatim
//          mode is inTEGER
//           On entry this describes how the eigenvalues are to
//           be specified:
//           mode = 0 means use D (with EI) as input
//           mode = 1 sets d1=1 and d(2:N)=1.0/cond
//           mode = 2 sets d(1:N-1)=1 and d(n)=1.0/cond
//           mode = 3 sets d(i)=cond**(-(I-1)/(N-1))
//           mode = 4 sets d(i)=1 - (i-1)/(N-1)*(1 - 1/cond)
//           mode = 5 sets D to random numbers in the range
//                    ( 1/cond, 1) such that their logarithms
//                    are uniformly distributed.  Each odd-even pair
//                    of elements will be either used as two real
//                    eigenvalues or as the real and imaginary part
//                    of a complex conjugate pair of eigenvalues;
//                    the choice of which is done is random, with
//                    50-50 probability, for each pair.
//           mode = 6 set D to random numbers from same distribution
//                    as the rest of the matrix.
//           mode < 0 has the same meaning as ABS(mode), except that
//              the order of the elements of D is reversed.
//           Thus if mode is between 1 and 4, D has entries ranging
//              from 1 to 1/cond, if between -1 and -4, D has entries
//              ranging from 1/cond to 1,
//           Not modified.
// \endverbatim
//
// \param[in] cond
// \verbatim
//          cond is DOUBLE PRECISION
//           On entry, this is used as described under mode above.
//           If used, it must be >= 1. Not modified.
// \endverbatim
//
// \param[in] dmax
// \verbatim
//          dmax is DOUBLE PRECISION
//           If mode is neither -6, 0 nor 6, the contents of d, as
//           computed according to mode and cond, will be scaled by
//           dmax / max(abs(d(i))).  Note that dmax need not be
//           positive: if dmax is negative (or zero), D will be
//           scaled by a negative number (or zero).
//           Not modified.
// \endverbatim
//
// \param[in] EI
// \verbatim
//          ei is CHARACTER*1 array, dimension ( N)
//           If mode is 0, and EI1 is not ' ' (space character),
//           this array specifies which elements of D (on input) are
//           real eigenvalues and which are the real and imaginary parts
//           of a complex conjugate pair of eigenvalues.  The elements
//           of ei may then only have the values 'R' and 'I'.  If
//           EI(j)='R' and EI(j+1)='I', then the j-th eigenvalue is
//           CMPLX( d(j), d(j+1)), and the (j+1)-th is the complex
//           conjugate thereof.  If EI(j)=EI(j+1)='R', then the j-th
//           eigenvalue is d(j) (i.e., real).  EI1 may not be 'I',
//           nor may two adjacent elements of ei both have the value 'I'.
//           If mode is not 0, then ei is ignored.  If mode is 0 and
//           EI1=' ', then the eigenvalues will all be real.
//           Not modified.
// \endverbatim
//
// \param[in] rsign
// \verbatim
//          rsign is CHARACTER*1
//           If mode is not 0, 6, or -6, and rsign='T', then the
//           elements of d, as computed according to mode and cond, will
//           be multiplied by a random sign (+1 or -1).  If rsign='F',
//           they will not be.  rsign may only have the values 'T' or
//           'F'.
//           Not modified.
// \endverbatim
//
// \param[in] upper
// \verbatim
//          upper is CHARACTER*1
//           If upper='T', then the elements of A above the diagonal
//           (and above the 2x2 diagonal blocks, if A has complex
//           eigenvalues) will be set to random numbers out of dist.
//           If upper='F', they will not.  upper may only have the
//           values 'T' or 'F'.
//           Not modified.
// \endverbatim
//
// \param[in] SIM
// \verbatim
//          SIM is CHARACTER*1
//           If SIM='T', then A will be operated on by a "similarity
//           transform", i.e., multiplied on the left by a matrix X and
//           on the right by X inverse.  X = U S V, where U and V are
//           random unitary matrices and S is a (diagonal) matrix of
//           singular values specified by DS, modes, and conds.  If
//           SIM='F', then A will not be transformed.
//           Not modified.
// \endverbatim
//
// \param[in,out] DS
// \verbatim
//          ds is DOUBLE PRECISION array, dimension ( N)
//           This array is used to specify the singular values of x,
//           in the same way that D specifies the eigenvalues of A.
//           If mode=0, the ds contains the singular values, which
//           may not be zero.
//           Modified if mode is nonzero.
// \endverbatim
//
// \param[in] modes
// \verbatim
//          modes is inTEGER
// \endverbatim
//
// \param[in] conds
// \verbatim
//          conds is DOUBLE PRECISION
//           Same as mode and cond, but for specifying the diagonal
//           of S.  modes=-6 and +6 are not allowed (since they would
//           result in randomly ill-conditioned eigenvalues.)
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is inTEGER
//           This specifies the lower bandwidth of the  matrix.  kl=1
//           specifies upper Hessenberg form.  If kl is at least N-1,
//           then A will have full lower bandwidth.  kl must be at
//           least 1.
//           Not modified.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is inTEGER
//           This specifies the upper bandwidth of the  matrix.  ku=1
//           specifies lower Hessenberg form.  If ku is at least N-1,
//           then A will have full upper bandwidth; if ku and kl
//           are both at least N-1, then A will be dense.  Only one of
//           ku and kl may be less than N-1.  ku must be at least 1.
//           Not modified.
// \endverbatim
//
// \param[in] anorm
// \verbatim
//          anorm is DOUBLE PRECISION
//           If anorm is not negative, then A will be scaled by a non-
//           negative real number to make the maximum-element-norm of A
//           to be anorm.
//           Not modified.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension ( lda, N)
//           On exit A is the desired test matrix.
//           Modified.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//           lda specifies the first dimension of A as declared in the
//           calling program.  lda must be at least N.
//           Not modified.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension ( 3*n)
//           workspace.
//           Modified.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is inTEGER
//           Error code.  On exit, info will be set to one of the
//           following values:
//             0 => normal return
//            -1 => N negative
//            -2 => dist illegal string
//            -5 => mode not in range -6 to 6
//            -6 => cond less than 1.0, and mode neither -6, 0 nor 6
//            -8 => EI1 is not ' ' or 'R', EI(j) is not 'R' or 'I', or
//                  two adjacent elements of ei are 'I'.
//            -9 => rsign is not 'T' or 'F'
//           -10 => upper is not 'T' or 'F'
//           -11 => SIM   is not 'T' or 'F'
//           -12 => modes=0 and ds has a zero singular value.
//           -13 => modes is not in the range -5 to 5.
//           -14 => modes is nonzero and conds is less than 1.
//           -15 => kl is less than 1.
//           -16 => ku is less than 1, or kl and ku are both less than
//                  N-1.
//           -19 => lda is less than N.
//            1  => Error return from Dlatm1 (computing D)
//            2  => Cannot scale to dmax (max. eigenvalue is 0)
//            3  => Error return from Dlatm1 (computing DS)
//            4  => Error return from Dlarge
//            5  => zero singular value from Dlatm1.
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
func dlatme(n *int, dist *byte, iseed *[]int, d *[]float64, mode *int, cond *float64, dmax *float64, ei *[]byte, rsign *byte, upper *byte, sim *byte, ds *[]float64, modes *int, conds *float64, kl *int, ku *int, anorm *float64, a *[][]float64, lda *int, work *[]float64, info *int) {
	zero := new(float64)
	one := new(float64)
	half := new(float64)
	badei := new(bool)
	bads := new(bool)
	useei := new(bool)
	i := new(int)
	ic := new(int)
	icols := new(int)
	idist := new(int)
	iinfo := new(int)
	ir := new(int)
	irows := new(int)
	irsign := new(int)
	isim := new(int)
	iupper := new(int)
	j := new(int)
	jc := new(int)
	jcr := new(int)
	jr := new(int)
	alpha := new(float64)
	tau := new(float64)
	temp := new(float64)
	xnorms := new(float64)
	tempa := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
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
	(*zero) = 0.0
	(*one) = 1.0
	(*half) = 1.0 / 2.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     1)      Decode and Test the input parameters.
	//             Initialize flags & seed.
	//
	(*(info)) = 0
	//
	//     Quick return if possible
	//
	if (*(n)) == 0 {
		return
	}
	//
	//     Decode dist
	//
	if blas.Lsame((dist), func() *byte {y := byte('U'); return &y }()) {
		(*idist) = 1
	} else if blas.Lsame((dist), func() *byte {y := byte('S'); return &y }()) {
		(*idist) = 2
	} else if blas.Lsame((dist), func() *byte {y := byte('N'); return &y }()) {
		(*idist) = 3
	} else {
		(*idist) = -1
	}
	//
	//     Check EI
	//
	(*useei) = true
	(*badei) = false
	if (*blas.Lsame(&((*(ei))[0]), func() *byte {y := byte(' '); return &y }())) || (*(mode)) != 0 {
		(*useei) = false
	} else {
		if blas.Lsame(&((*(ei))[0]), func() *byte {y := byte('R'); return &y }()) {
			for (*j) = 2; (*j) <= (*(n)); (*j)++ {
				if blas.Lsame(&((*(ei))[(*j)-1]), func() *byte {y := byte('I'); return &y }()) {
					if blas.Lsame(&((*(ei))[(*j)-0]), func() *byte {y := byte('I'); return &y }()) {
						(*badei) = true
					}
				} else {
					if !blas.Lsame(&((*(ei))[(*j)-1]), func() *byte {y := byte('R'); return &y }()) {
						(*badei) = true
					}
				}
				//Label10:
			}
		} else {
			(*badei) = true
		}
	}
	//
	//     Decode rsign
	//
	if blas.Lsame((rsign), func() *byte {y := byte('T'); return &y }()) {
		(*irsign) = 1
	} else if blas.Lsame((rsign), func() *byte {y := byte('F'); return &y }()) {
		(*irsign) = 0
	} else {
		(*irsign) = -1
	}
	//
	//     Decode upper
	//
	if blas.Lsame((upper), func() *byte {y := byte('T'); return &y }()) {
		(*iupper) = 1
	} else if blas.Lsame((upper), func() *byte {y := byte('F'); return &y }()) {
		(*iupper) = 0
	} else {
		(*iupper) = -1
	}
	//
	//     Decode SIM
	//
	if blas.Lsame((SIM), func() *byte {y := byte('T'); return &y }()) {
		(*isim) = 1
	} else if blas.Lsame((SIM), func() *byte {y := byte('F'); return &y }()) {
		(*isim) = 0
	} else {
		(*isim) = -1
	}
	//
	//     Check DS, if modes=0 and isim=1
	//
	(*bads) = false
	if (*(modes)) == 0 && (*isim) == 1 {
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			if (*(ds))[(*j)-1] == (*zero) {
				(*bads) = true
			}
			//Label20:
		}
	}
	//
	//     Set info if an error
	//
	if (*(n)) < 0 {
		(*(info)) = -1
	} else if (*idist) == -1 {
		(*(info)) = -2
	} else if (ABS((*(mode)))) > 6 {
		(*(info)) = -5
	} else if ((*(mode)) != 0 && ABS((*(mode))) != 6) && (*(cond)) < (*one) {
		(*(info)) = -6
	} else if *BADei {
		(*(info)) = -8
	} else if (*irsign) == -1 {
		(*(info)) = -9
	} else if (*iupper) == -1 {
		(*(info)) = -10
	} else if (*isim) == -1 {
		(*(info)) = -11
	} else if *bads {
		(*(info)) = -12
	} else if (*isim) == 1 && ABS((*(modes))) > 5 {
		(*(info)) = -13
	} else if (*isim) == 1 && (*(modes)) != 0 && (*(conds)) < (*one) {
		(*(info)) = -14
	} else if (*(kl)) < 1 {
		(*(info)) = -15
	} else if (*(ku)) < 1 || ((*(ku)) < (*(n))-1 && (*(kl)) < (*(n))-1) {
		(*(info)) = -16
	} else if (*(lda)) < (MAX(1, (*(n)))) {
		(*(info)) = -19
	}
	//
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("dlatme"); return &y }(), -(*(info)))
		return
	}
	//
	//     Initialize random number generator
	//
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*(iseed))[(*i)-1] = (MOD(ABS(((*(iseed))[(*i)-1])), int(4096)))
		//Label30:
	}
	//
	if (MOD(((*(iseed))[3]), int(2))) != 1 {
		(*(iseed))[3] = (*(iseed))[3] + 1
	}
	//
	//     2)      Set up diagonal of A
	//
	//             Compute D according to cond and mode
	//
	Dlatm1((mode), (cond), irsign, idist, (iseed), (d), (n), iinfo)
	if (*iinfo) != 0 {
		(*(info)) = 1
		return
	}
	if (*(mode)) != 0 && ABS((*(mode))) != 6 {
		//
		//        Scale by dmax
		//
		(*temp) = (ABS(((*(d))[0])))
		for (*i) = 2; (*i) <= (*(n)); (*i)++ {
			(*temp) = (MAX((*temp), ABS(((*(d))[(*i)-1]))))
			//Label40:
		}
		//
		if (*temp) > (*zero) {
			(*alpha) = (*(dmax)) / (*temp)
		} else if (*(dmax)) != (*zero) {
			(*(info)) = 2
			return
		} else {
			(*alpha) = (*zero)
		}
		//
		Dscal((n), alpha, (d), func() *int {y := 1; return &y }())
		//
	}
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (n), (n), zero, zero, (a), (lda))
	Dcopy((n), (d), func() *int {y := 1; return &y }(), (a), (*(lda))+1)
	//
	//     Set up complex conjugate pairs
	//
	if (*(mode)) == 0 {
		if *useei {
			for (*j) = 2; (*j) <= (*(n)); (*j)++ {
				if blas.Lsame(&((*(ei))[(*j)-1]), func() *byte {y := byte('I'); return &y }()) {
					(*(a))[(*j)-0][(*j)-1] = (*(a))[(*j)-1][(*j)-1]
					(*(a))[(*j)-1][(*j)-0] = -(*(a))[(*j)-1][(*j)-1]
					(*(a))[(*j)-1][(*j)-1] = (*(a))[(*j)-0][(*j)-0]
				}
				//Label50:
			}
		}
		//
	} else if (ABS((*(mode)))) == 5 {
		//
		for (*j) = 2; (*j) <= (*(n)); (*j) += 2 {
			if (*Dlaran((iseed))) > (*half) {
				(*(a))[(*j)-0][(*j)-1] = (*(a))[(*j)-1][(*j)-1]
				(*(a))[(*j)-1][(*j)-0] = -(*(a))[(*j)-1][(*j)-1]
				(*(a))[(*j)-1][(*j)-1] = (*(a))[(*j)-0][(*j)-0]
			}
			//Label60:
		}
	}
	//
	//     3)      If upper='T', set upper triangle of A to random numbers.
	//             (but don't modify the corners of 2x2 blocks.)
	//
	if (*iupper) != 0 {
		for (*jc) = 2; (*jc) <= (*(n)); (*jc)++ {
			if (*(a))[(*jc)-0][(*jc)-1] != (*zero) {
				(*jr) = (*jc) - 2
			} else {
				(*jr) = (*jc) - 1
			}
			Dlarnv(idist, (iseed), jr, &((*(a))[0][(*jc)-1]))
			//Label70:
		}
	}
	//
	//     4)      If SIM='T', apply similarity transformation.
	//
	//                                -1
	//             Transform is  X A X , where X = U S V, thus
	//
	//             it is  U S V A V' (1/S) U'
	//
	if (*isim) != 0 {
		//
		//        Compute S (singular values of the eigenvector matrix)
		//        according to conds and modes
		//
		Dlatm1((modes), (conds), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), (iseed), (ds), (n), iinfo)
		if (*iinfo) != 0 {
			(*(info)) = 3
			return
		}
		//
		//        Multiply by V and V'
		//
		Dlarge((n), (a), (lda), (iseed), (work), iinfo)
		if (*iinfo) != 0 {
			(*(info)) = 4
			return
		}
		//
		//        Multiply by S and (1/S)
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			Dscal((n), &((*(ds))[(*j)-1]), &((*(a))[(*j)-1][0]), (lda))
			if (*(ds))[(*j)-1] != (*zero) {
				Dscal((n), (*one)/(*(ds))[(*j)-1], &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
			} else {
				(*(info)) = 5
				return
			}
			//Label80:
		}
		//
		//        Multiply by U and U'
		//
		Dlarge((n), (a), (lda), (iseed), (work), iinfo)
		if (*iinfo) != 0 {
			(*(info)) = 4
			return
		}
	}
	//
	//     5)      Reduce the bandwidth.
	//
	if (*(kl)) < (*(n))-1 {
		//
		//        Reduce bandwidth -- kill column
		//
		for (*jcr) = (*(kl)) + 1; (*jcr) <= (*(n))-1; (*jcr)++ {
			(*ic) = (*jcr) - (*(kl))
			(*irows) = (*(n)) + 1 - (*jcr)
			(*icols) = (*(n)) + (*(kl)) - (*jcr)
			//
			Dcopy(irows, &((*(a))[(*jcr)-1][(*ic)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
			(*xnorms) = (*(work))[0]
			dlarfG(irows, xnorms, &((*(work))[1]), func() *int {y := 1; return &y }(), tau)
			(*(work))[0] = (*one)
			//
			Dgemv(func() *byte {y := byte('T'); return &y }(), irows, ICOLS, one, &((*(a))[(*jcr)-1][(*ic)+0]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*irows)+0]), func() *int {y := 1; return &y }())
			Dger(irows, ICOLS, -(*tau), (work), func() *int {y := 1; return &y }(), &((*(work))[(*irows)+0]), func() *int {y := 1; return &y }(), &((*(a))[(*jcr)-1][(*ic)+0]), (lda))
			//
			Dgemv(func() *byte {y := byte('N'); return &y }(), (n), irows, one, &((*(a))[0][(*jcr)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*irows)+0]), func() *int {y := 1; return &y }())
			Dger((n), irows, -(*tau), &((*(work))[(*irows)+0]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[0][(*jcr)-1]), (lda))
			//
			(*(a))[(*jcr)-1][(*ic)-1] = (*xnorms)
			Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (*irows)-1, func() *int {y := 1; return &y }(), zero, zero, &((*(a))[(*jcr)+0][(*ic)-1]), (lda))
			//Label90:
		}
	} else if (*(ku)) < (*(n))-1 {
		//
		//        Reduce upper bandwidth -- kill a row at a time.
		//
		for (*jcr) = (*(ku)) + 1; (*jcr) <= (*(n))-1; (*jcr)++ {
			(*ir) = (*jcr) - (*(ku))
			(*irows) = (*(n)) + (*(ku)) - (*jcr)
			(*icols) = (*(n)) + 1 - (*jcr)
			//
			Dcopy(ICOLS, &((*(a))[(*ir)-1][(*jcr)-1]), (lda), (work), func() *int {y := 1; return &y }())
			(*xnorms) = (*(work))[0]
			dlarfG(ICOLS, xnorms, &((*(work))[1]), func() *int {y := 1; return &y }(), tau)
			(*(work))[0] = (*one)
			//
			Dgemv(func() *byte {y := byte('N'); return &y }(), irows, ICOLS, one, &((*(a))[(*ir)+0][(*jcr)-1]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*icols)+0]), func() *int {y := 1; return &y }())
			Dger(irows, ICOLS, -(*tau), &((*(work))[(*icols)+0]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }(), &((*(a))[(*ir)+0][(*jcr)-1]), (lda))
			//
			Dgemv(func() *byte {y := byte('C'); return &y }(), ICOLS, (n), one, &((*(a))[(*jcr)-1][0]), (lda), (work), func() *int {y := 1; return &y }(), zero, &((*(work))[(*icols)+0]), func() *int {y := 1; return &y }())
			Dger(ICOLS, (n), -(*tau), (work), func() *int {y := 1; return &y }(), &((*(work))[(*icols)+0]), func() *int {y := 1; return &y }(), &((*(a))[(*jcr)-1][0]), (lda))
			//
			(*(a))[(*ir)-1][(*jcr)-1] = (*xnorms)
			Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), func() *int {y := 1; return &y }(), (*icols)-1, zero, zero, &((*(a))[(*ir)-1][(*jcr)+0]), (lda))
			//Label100:
		}
	}
	//
	//     Scale the matrix to have norm anorm
	//
	if (*(anorm)) >= (*zero) {
		(*temp) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (n), (n), (a), (lda), tempa))
		if (*temp) > (*zero) {
			(*alpha) = (*(anorm)) / (*temp)
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				Dscal((n), alpha, &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
				//Label110:
			}
		}
	}
	//
	return
	//
	//     End of dlatme
	//
}
