package goblas

import (
	"github.com/whipstein/golapack/blas"
)

//    Dlatms generates random matrices with specified singular values
//    (or symmetric/hermitian with specified eigenvalues)
//    for testing lapACK programs.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlatms( m, n, dist, iseed, sym, d, mode, cond, dmax,
//                          kl, ku, PACK, a, lda, work, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          dist, PACK, sym
//       intEGER            info, kl, ku, lda, m, mode, N
//       DOUBLE PRECISION   cond, dmax
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), d(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Dlatms generates random matrices with specified singular values
//    (or symmetric/hermitian with specified eigenvalues)
//    for testing lapACK programs.
//
//    Dlatms operates by applying the following sequence of
//    operations:
//
//      Set the diagonal to d, where D may be input or
//         computed according to mode, cond, dmax, and sym
//         as described below.
//
//      Generate a matrix with the appropriate band structure, by one
//         of two methods:
//
//      Method A:
//          Generate a dense M x N matrix by multiplying D on the left
//              and the right by random unitary matrices, then:
//
//          Reduce the bandwidth according to kl and ku, using
//          Householder transformations.
//
//      Method B:
//          Convert the bandwidth-0 (i.e., diagonal) matrix to a
//              bandwidth-1 matrix using Givens rotations, "chasing"
//              out-of-band elements back, much as in QR; then
//              convert the bandwidth-1 to a bandwidth-2 matrix, etc.
//              Note that for reasonably small bandwidths (relative to
//              M and N) this requires less storage, as a dense matrix
//              is not generated.  Also, for symmetric matrices, only
//              one triangle is generated.
//
//      Method A is chosen if the bandwidth is a large fraction of the
//          order of the matrix, and lda is at least M (so a dense
//          matrix can be stored.)  Method B is chosen if the bandwidth
//          is small (< 1/2 N for symmetric, < .3 N+M for
//          non-symmetric), or lda is less than M and not less than the
//          bandwidth.
//
//      Pack the matrix if desired. Options specified by PACK are:
//         no packing
//         zero out upper half (if symmetric)
//         zero out lower half (if symmetric)
//         store the upper half columnwise (if symmetric or upper
//               triangular)
//         store the lower half columnwise (if symmetric or lower
//               triangular)
//         store the lower triangle in banded format (if symmetric
//               or lower triangular)
//         store the upper triangle in banded format (if symmetric
//               or upper triangular)
//         store the entire matrix in banded format
//      If Method B is chosen, and band format is specified, then the
//         matrix will be generated in the band format, so no repacking
//         will be necessary.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//           The number of rows of A. Not modified.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//           The number of columns of A. Not modified.
// \endverbatim
//
// \param[in] dist
// \verbatim
//          dist is CHARACTER*1
//           On entry, dist specifies the type of distribution to be used
//           to generate the random eigen-/singular values.
//           'U' => UNIFORM( 0, 1)  ( 'U' for uniform)
//           'S' => UNIFORM( -1, 1) ( 'S' for symmetric)
//           'N' => normaL( 0, 1)   ( 'N' for normal)
//           Not modified.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is intEGER array, dimension ( 4)
//           On entry iseed specifies the seed of the random number
//           generator. They should lie between 0 and 4095 inclusive,
//           and iseed(4) should be odd. The random number generator
//           uses a linear congruential sequence limited to small
//           integers, and so should produce machine independent
//           random numbers. The values of iseed are changed on
//           exit, and can be used in the next call to Dlatms
//           to continue the same random number sequence.
//           Changed on exit.
// \endverbatim
//
// \param[in] sym
// \verbatim
//          sym is CHARACTER*1
//           If sym='S' or 'H', the generated matrix is symmetric, with
//             eigenvalues specified by d, cond, mode, and dmax; they
//             may be positive, negative, or zero.
//           If sym='P', the generated matrix is symmetric, with
//             eigenvalues (= singular values) specified by d, cond,
//             mode, and dmax; they will not be negative.
//           If sym='N', the generated matrix is nonsymmetric, with
//             singular values specified by d, cond, mode, and dmax;
//             they will not be negative.
//           Not modified.
// \endverbatim
//
// \param[in,out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension ( Min( M, N))
//           This array is used to specify the singular values or
//           eigenvalues of A (see sym, above.)  If mode=0, then D is
//           assumed to contain the singular/eigenvalues, otherwise
//           they will be computed according to mode, cond, and dmax,
//           and placed in D.
//           Modified if mode is nonzero.
// \endverbatim
//
// \param[in] mode
// \verbatim
//          mode is intEGER
//           On entry this describes how the singular/eigenvalues are to
//           be specified:
//           mode = 0 means use D as input
//           mode = 1 sets d1=1 and d(2:N)=1.0/cond
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
//           If sym='S' or 'H', and mode is neither 0, 6, nor -6, then
//              the elements of D will also be multiplied by a random
//              sign (i.e., +1 or -1.)
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
//           dmax / max(abs(d(i))); thus, the maximum absolute eigen- or
//           singular value (which is to say the norm) will be ABS(dmax).
//           Note that dmax need not be positive: if dmax is negative
//           (or zero), D will be scaled by a negative number (or zero).
//           Not modified.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is intEGER
//           This specifies the lower bandwidth of the  matrix. For
//           example, kl=0 implies upper triangular, kl=1 implies upper
//           Hessenberg, and kl being at least M-1 means that the matrix
//           has full lower bandwidth.  kl must equal ku if the matrix
//           is symmetric.
//           Not modified.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//           This specifies the upper bandwidth of the  matrix. For
//           example, ku=0 implies lower triangular, ku=1 implies lower
//           Hessenberg, and ku being at least N-1 means that the matrix
//           has full upper bandwidth.  kl must equal ku if the matrix
//           is symmetric.
//           Not modified.
// \endverbatim
//
// \param[in] PACK
// \verbatim
//          PACK is CHARACTER*1
//           This specifies packing of matrix as follows:
//           'N' => no packing
//           'U' => zero out all subdiagonal entries (if symmetric)
//           'L' => zero out all superdiagonal entries (if symmetric)
//           'C' => store the upper triangle columnwise
//                  (only if the matrix is symmetric or upper triangular)
//           'R' => store the lower triangle columnwise
//                  (only if the matrix is symmetric or lower triangular)
//           'B' => store the lower triangle in band storage scheme
//                  (only if matrix symmetric or lower triangular)
//           'Q' => store the upper triangle in band storage scheme
//                  (only if matrix symmetric or upper triangular)
//           'Z' => store the entire matrix in band storage scheme
//                      (pivoting can be provided for by using this
//                      option to store A in the trailing rows of
//                      the allocated storage)
//
//           Using these options, the various lapACK packed and banded
//           storage schemes can be obtained:
//           GB               - use 'Z'
//           PB, SB or TB     - use 'B' or 'Q'
//           PP, SP or TP     - use 'C' or 'R'
//
//           If two calls to Dlatms differ only in the PACK parameter,
//           they will generate mathematically equivalent matrices.
//           Not modified.
// \endverbatim
//
// \param[in,out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension ( lda, N)
//           On exit A is the desired test matrix.  A is first generated
//           in full (unpacked) form, and then packed, if so specified
//           by PACK.  Thus, the first M elements of the first N
//           columns will always be modified.  If PACK specifies a
//           packed or banded storage scheme, all lda elements of the
//           first N columns will be modified; the elements of the
//           array which do not correspond to elements of the generated
//           matrix are set to zero.
//           Modified.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//           lda specifies the first dimension of A as declared in the
//           calling program.  If PACK='N', 'U', 'L', 'C', or 'R', then
//           lda must be at least M.  If PACK='B' or 'Q', then lda must
//           be at least Min( kl, M-1) (which is equal to Min(ku,N-1)).
//           If PACK='Z', lda must be large enough to hold the packed
//           array: Min( ku, N-1) + Min( kl, M-1) + 1.
//           Not modified.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension ( 3*MAX( N, M))
//           workspace.
//           Modified.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
//           Error code.  On exit, info will be set to one of the
//           following values:
//             0 => normal return
//            -1 => M negative or unequal to N and sym='S', 'H', or 'P'
//            -2 => N negative
//            -3 => dist illegal string
//            -5 => sym illegal string
//            -7 => mode not in range -6 to 6
//            -8 => cond less than 1.0, and mode neither -6, 0 nor 6
//           -10 => kl negative
//           -11 => ku negative, or sym='S' or 'H' and ku not equal to kl
//           -12 => PACK illegal string, or PACK='U' or 'L', and sym='N';
//                  or PACK='C' or 'Q' and sym='N' and kl is not zero;
//                  or PACK='R' or 'B' and sym='N' and ku is not zero;
//                  or PACK='U', 'L', 'C', 'R', 'B', or 'Q', and M is not
//                  N.
//           -14 => lda is less than m, or PACK='Z' and lda is less than
//                  Min(ku,N-1) + Min(kl,M-1) + 1.
//            1  => Error return from Dlatm1
//            2  => Cannot scale to dmax (max. sing. value is 0)
//            3  => Error return from Dlagge or SLAGSY
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
func Dlatms(m *int, n *int, dist *byte, iseed *[]int, sym *byte, d *[]float64, mode *int, cond *float64, dmax *float64, kl *int, ku *int, pack *byte, a *[][]float64, lda *int, work *[]float64, info *int) {
	zero := new(float64)
	one := new(float64)
	twopi := new(float64)
	givens := new(bool)
	ilextr := new(bool)
	iltemp := new(bool)
	topdwn := new(bool)
	i := new(int)
	ic := new(int)
	icol := new(int)
	idist := new(int)
	iendch := new(int)
	iinfo := new(int)
	il := new(int)
	ilda := new(int)
	ioffg := new(int)
	ioffst := new(int)
	ipack := new(int)
	ipackg := new(int)
	ir := new(int)
	ir1 := new(int)
	ir2 := new(int)
	irow := new(int)
	irsign := new(int)
	iskew := new(int)
	isym := new(int)
	isympk := new(int)
	j := new(int)
	jc := new(int)
	jch := new(int)
	jkl := new(int)
	jku := new(int)
	jr := new(int)
	k := new(int)
	llb := new(int)
	minlda := new(int)
	mnmin := new(int)
	mr := new(int)
	nc := new(int)
	uub := new(int)
	alpha := new(float64)
	angle := new(float64)
	c := new(float64)
	dummy := new(float64)
	extra := new(float64)
	s := new(float64)
	temp := new(float64)
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
	(*twopi) = 6.2831853071795864769252867663e+0
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
	//     1)      Decode and Test the input parameters.
	//             Initialize flags & seed.
	//
	(*(info)) = 0
	//
	//     Quick return if possible
	//
	if (*(m)) == 0 || (*(n)) == 0 {
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
	//     Decode sym
	//
	if blas.Lsame((sym), func() *byte {y := byte('N'); return &y }()) {
		(*isym) = 1
		(*irsign) = 0
	} else if blas.Lsame((sym), func() *byte {y := byte('P'); return &y }()) {
		(*isym) = 2
		(*irsign) = 0
	} else if blas.Lsame((sym), func() *byte {y := byte('S'); return &y }()) {
		(*isym) = 2
		(*irsign) = 1
	} else if blas.Lsame((sym), func() *byte {y := byte('H'); return &y }()) {
		(*isym) = 2
		(*irsign) = 1
	} else {
		(*isym) = -1
	}
	//
	//     Decode PACK
	//
	(*isympk) = 0
	if blas.Lsame((pack), func() *byte {y := byte('N'); return &y }()) {
		(*ipack) = 0
	} else if blas.Lsame((pack), func() *byte {y := byte('U'); return &y }()) {
		(*ipack) = 1
		(*isympk) = 1
	} else if blas.Lsame((pack), func() *byte {y := byte('L'); return &y }()) {
		(*ipack) = 2
		(*isympk) = 1
	} else if blas.Lsame((pack), func() *byte {y := byte('C'); return &y }()) {
		(*ipack) = 3
		(*isympk) = 2
	} else if blas.Lsame((pack), func() *byte {y := byte('R'); return &y }()) {
		(*ipack) = 4
		(*isympk) = 3
	} else if blas.Lsame((pack), func() *byte {y := byte('B'); return &y }()) {
		(*ipack) = 5
		(*isympk) = 3
	} else if blas.Lsame((pack), func() *byte {y := byte('Q'); return &y }()) {
		(*ipack) = 6
		(*isympk) = 2
	} else if blas.Lsame((pack), func() *byte {y := byte('Z'); return &y }()) {
		(*ipack) = 7
	} else {
		(*ipack) = -1
	}
	//
	//     Set certain internal parameters
	//
	(*mnmin) = (Min((*(m)), (*(n))))
	(*llb) = (Min((*(kl)), (*(m))-1))
	(*uub) = (Min((*(ku)), (*(n))-1))
	(*mr) = (Min((*(m)), (*(n))+(*llb)))
	(*nc) = (Min((*(n)), (*(m))+(*uub)))
	//
	if (*ipack) == 5 || (*ipack) == 6 {
		(*minlda) = (*uub) + 1
	} else if (*ipack) == 7 {
		(*minlda) = (*llb) + (*uub) + 1
	} else {
		(*minlda) = (*(m))
	}
	//
	//     Use Givens rotation method if bandwidth small enough,
	//     or if lda is too small to store the matrix unpacked.
	//
	(*givens) = false
	if (*isym) == 1 {
		if (DBLE((*llb) + (*uub))) < 0.3*DBLE(MAX(1, (*mr)+(*nc))) {
			(*givens) = true
		}
	} else {
		if 2*(*llb) < (*(m)) {
			(*givens) = true
		}
	}
	if (*(lda)) < (*(m)) && (*(lda)) >= (*minlda) {
		(*givens) = true
	}
	//
	//     Set info if an error
	//
	if (*(m)) < 0 {
		(*(info)) = -1
	} else if (*(m)) != (*(n)) && (*isym) != 1 {
		(*(info)) = -1
	} else if (*(n)) < 0 {
		(*(info)) = -2
	} else if (*idist) == -1 {
		(*(info)) = -3
	} else if (*isym) == -1 {
		(*(info)) = -5
	} else if (ABS((*(mode)))) > 6 {
		(*(info)) = -7
	} else if ((*(mode)) != 0 && ABS((*(mode))) != 6) && (*(cond)) < (*one) {
		(*(info)) = -8
	} else if (*(kl)) < 0 {
		(*(info)) = -10
	} else if (*(ku)) < 0 || ((*isym) != 1 && (*(kl)) != (*(ku))) {
		(*(info)) = -11
	} else if (*ipack) == -1 || ((*isympk) == 1 && (*isym) == 1) || ((*isympk) == 2 && (*isym) == 1 && (*(kl)) > 0) || ((*isympk) == 3 && (*isym) == 1 && (*(ku)) > 0) || ((*isympk) != 0 && (*(m)) != (*(n))) {
		(*(info)) = -12
	} else if (*(lda)) < (MAX(1, (*minlda))) {
		(*(info)) = -14
	}
	//
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("Dlatms"); return &y }(), -(*(info)))
		return
	}
	//
	//     Initialize random number generator
	//
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*(iseed))[(*i)-1] = (MOD(ABS(((*(iseed))[(*i)-1])), int(4096)))
		//Label10:
	}
	//
	if (MOD(((*(iseed))[3]), int(2))) != 1 {
		(*(iseed))[3] = (*(iseed))[3] + 1
	}
	//
	//     2)      Set up D  if indicated.
	//
	//             Compute D according to cond and mode
	//
	Dlatm1((mode), (cond), irsign, idist, (iseed), (d), mnmin, iinfo)
	if (*iinfo) != 0 {
		(*(info)) = 1
		return
	}
	//
	//     Choose Top-Down if D is (apparently) increasing,
	//     Bottom-Up if D is (apparently) decreasing.
	//
	if (ABS(((*(d))[0]))) <= (ABS(((*(d))[(*mnmin)-1]))) {
		(*topdwn) = true
	} else {
		(*topdwn) = false
	}
	//
	if (*(mode)) != 0 && ABS((*(mode))) != 6 {
		//
		//        Scale by dmax
		//
		(*temp) = (ABS(((*(d))[0])))
		for (*i) = 2; (*i) <= (*mnmin); (*i)++ {
			(*temp) = (MAX((*temp), ABS(((*(d))[(*i)-1]))))
			//Label20:
		}
		//
		if (*temp) > (*zero) {
			(*alpha) = (*(dmax)) / (*temp)
		} else {
			(*(info)) = 2
			return
		}
		//
		Dscal(mnmin, alpha, (d), func() *int {y := 1; return &y }())
		//
	}
	//
	//     3)      Generate Banded Matrix using Givens rotations.
	//             Also the special case of uub=llb=0
	//
	//               Compute Addressing _constants to cover all
	//               storage formats.  Whether GE, SY, GB, or SB,
	//               upper or lower triangle or both,
	//               the (i,j)-th element is in
	//               a( i - iskew*j + ioffst, j)
	//
	if (*ipack) > 4 {
		(*ilda) = (*(lda)) - 1
		(*iskew) = 1
		if (*ipack) > 5 {
			(*ioffst) = (*uub) + 1
		} else {
			(*ioffst) = 1
		}
	} else {
		(*ilda) = (*(lda))
		(*iskew) = 0
		(*ioffst) = 0
	}
	//
	//     ipackg is the format that the matrix is generated in. If this is
	//     different from ipACK, then the matrix must be repacked at the
	//     end.  It also signals how to compute the norm, for scaling.
	//
	(*ipackg) = 0
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (lda), (n), zero, zero, (a), (lda))
	//
	//     diagonal Matrix -- We are done, unless it
	//     is to be stored SP/PP/TP (PACK='R' or 'C')
	//
	if (*llb) == 0 && (*uub) == 0 {
		Dcopy(mnmin, (d), func() *int {y := 1; return &y }(), &((*(a))[1-(*iskew)+(*ioffst)-1][0]), (*ilda)+1)
		if (*ipack) <= 2 || (*ipack) >= 5 {
			(*ipackg) = (*ipack)
		}
		//
	} else if *givens {
		//
		//        Check whether to use Givens rotations,
		//        Householder transformations, or nothing.
		//
		if (*isym) == 1 {
			//
			//           Non-symmetric -- A = U D V
			//
			if (*ipack) > 4 {
				(*ipackg) = (*ipack)
			} else {
				(*ipackg) = 0
			}
			//
			Dcopy(mnmin, (d), func() *int {y := 1; return &y }(), &((*(a))[1-(*iskew)+(*ioffst)-1][0]), (*ilda)+1)
			//
			if *topdwn {
				(*jkl) = 0
				for (*jku) = 1; (*jku) <= (*uub); (*jku)++ {
					//
					//                 Transform from bandwidth jkl, jku-1 to jkl, jku
					//
					//                 Last row actually rotated is M
					//                 Last column actually rotated is Min( M+jku, N)
					//
					for (*jr) = 1; (*jr) <= Min((*(m))+(*jku), (*(n)))+(*jkl)-1; (*jr)++ {
						(*extra) = (*zero)
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = (*Sin(angle))
						(*icol) = (MAX(1, (*jr)-(*jkl)))
						if (*jr) < (*(m)) {
							(*il) = Min((*(n)), (*jr)+(*jku)) + 1 - (*icol)
							Dlarot(&(true), (*jr) > (*jkl), &(false), il, c, s, &((*(a))[(*jr)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, extra, dummy)
						}
						//
						//                    Chase "extra" back up
						//
						(*ir) = (*jr)
						(*ic) = (*icol)
						for (*jch) = (*jr) - (*jkl); (*jch) <= 1; (*jch) += -(*jkl) - (*jku) {
							if (*ir) < (*(m)) {
								dlartg(&((*(a))[(*ir)+1-(*iskew)*((*ic)+1)+(*ioffst)-1][(*ic)+0]), extra, c, s, dummy)
							}
							(*irow) = (MAX(1, (*jch)-(*jku)))
							(*il) = (*ir) + 2 - (*irow)
							(*temp) = (*zero)
							(*iltemp) = (*jch) > (*jku)
							Dlarot(&(false), iltemp, &(true), il, c, -(*s), &((*(a))[(*irow)-(*iskew)*(*ic)+(*ioffst)-1][(*ic)-1]), ilda, temp, extra)
							if *iltemp {
								dlartg(&((*(a))[(*irow)+1-(*iskew)*((*ic)+1)+(*ioffst)-1][(*ic)+0]), temp, c, s, dummy)
								(*icol) = (MAX(1, (*jch)-(*jku)-(*jkl)))
								(*il) = (*ic) + 2 - (*icol)
								(*extra) = (*zero)
								Dlarot(&(true), (*jch) > (*jku)+(*jkl), &(true), il, c, -(*s), &((*(a))[(*irow)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, extra, temp)
								(*ic) = (*icol)
								(*ir) = (*irow)
							}
							//Label30:
						}
						//Label40:
					}
					//Label50:
				}
				//
				(*jku) = (*uub)
				for (*jkl) = 1; (*jkl) <= (*llb); (*jkl)++ {
					//
					//                 Transform from bandwidth jkl-1, jku to jkl, jku
					//
					for (*jc) = 1; (*jc) <= Min((*(n))+(*jkl), (*(m)))+(*jku)-1; (*jc)++ {
						(*extra) = (*zero)
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = (*Sin(angle))
						(*irow) = (MAX(1, (*jc)-(*jku)))
						if (*jc) < (*(n)) {
							(*il) = Min((*(m)), (*jc)+(*jkl)) + 1 - (*irow)
							Dlarot(&(false), (*jc) > (*jku), &(false), il, c, s, &((*(a))[(*irow)-(*iskew)*(*jc)+(*ioffst)-1][(*jc)-1]), ilda, extra, dummy)
						}
						//
						//                    Chase "extra" back up
						//
						(*ic) = (*jc)
						(*ir) = (*irow)
						for (*jch) = (*jc) - (*jku); (*jch) <= 1; (*jch) += -(*jkl) - (*jku) {
							if (*ic) < (*(n)) {
								dlartg(&((*(a))[(*ir)+1-(*iskew)*((*ic)+1)+(*ioffst)-1][(*ic)+0]), extra, c, s, dummy)
							}
							(*icol) = (MAX(1, (*jch)-(*jkl)))
							(*il) = (*ic) + 2 - (*icol)
							(*temp) = (*zero)
							(*iltemp) = (*jch) > (*jkl)
							Dlarot(&(true), iltemp, &(true), il, c, -(*s), &((*(a))[(*ir)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, temp, extra)
							if *iltemp {
								dlartg(&((*(a))[(*ir)+1-(*iskew)*((*icol)+1)+(*ioffst)-1][(*icol)+0]), temp, c, s, dummy)
								(*irow) = (MAX(1, (*jch)-(*jkl)-(*jku)))
								(*il) = (*ir) + 2 - (*irow)
								(*extra) = (*zero)
								Dlarot(&(false), (*jch) > (*jkl)+(*jku), &(true), il, c, -(*s), &((*(a))[(*irow)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, extra, temp)
								(*ic) = (*icol)
								(*ir) = (*irow)
							}
							//Label60:
						}
						//Label70:
					}
					//Label80:
				}
				//
			} else {
				//
				//              Bottom-Up -- Start at the bottom right.
				//
				(*jkl) = 0
				for (*jku) = 1; (*jku) <= (*uub); (*jku)++ {
					//
					//                 Transform from bandwidth jkl, jku-1 to jkl, jku
					//
					//                 First row actually rotated is M
					//                 First column actually rotated is Min( M+jku, N)
					//
					(*iendch) = Min((*(m)), (*(n))+(*jkl)) - 1
					for (*jc) = Min((*(m))+(*jku), (*(n))) - 1; (*jc) <= 1-(*jkl); (*jc) += -1 {
						(*extra) = (*zero)
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = (*Sin(angle))
						(*irow) = (MAX(1, (*jc)-(*jku)+1))
						if (*jc) > 0 {
							(*il) = Min((*(m)), (*jc)+(*jkl)+1) + 1 - (*irow)
							Dlarot(&(false), &(false), (*jc)+(*jkl) < (*(m)), il, c, s, &((*(a))[(*irow)-(*iskew)*(*jc)+(*ioffst)-1][(*jc)-1]), ilda, dummy, extra)
						}
						//
						//                    Chase "extra" back down
						//
						(*ic) = (*jc)
						for (*jch) = (*jc) + (*jkl); (*jch) <= (*iendch); (*jch) += (*jkl) + (*jku) {
							(*ilextr) = (*ic) > 0
							if *ilextr {
								dlartg(&((*(a))[(*jch)-(*iskew)*(*ic)+(*ioffst)-1][(*ic)-1]), extra, c, s, dummy)
							}
							(*ic) = (MAX(1, (*ic)))
							(*icol) = (Min((*(n))-1, (*jch)+(*jku)))
							(*iltemp) = (*jch)+(*jku) < (*(n))
							(*temp) = (*zero)
							Dlarot(&(true), ilextr, iltemp, (*icol)+2-(*ic), c, s, &((*(a))[(*jch)-(*iskew)*(*ic)+(*ioffst)-1][(*ic)-1]), ilda, extra, temp)
							if *iltemp {
								dlartg(&((*(a))[(*jch)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), temp, c, s, dummy)
								(*il) = Min((*iendch), (*jch)+(*jkl)+(*jku)) + 2 - (*jch)
								(*extra) = (*zero)
								Dlarot(&(false), &(true), (*jch)+(*jkl)+(*jku) <= (*iendch), il, c, s, &((*(a))[(*jch)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, temp, extra)
								(*ic) = (*icol)
							}
							//Label90:
						}
						//Label100:
					}
					//Label110:
				}
				//
				(*jku) = (*uub)
				for (*jkl) = 1; (*jkl) <= (*llb); (*jkl)++ {
					//
					//                 Transform from bandwidth jkl-1, jku to jkl, jku
					//
					//                 First row actually rotated is Min( N+jkl, M)
					//                 First column actually rotated is N
					//
					(*iendch) = Min((*(n)), (*(m))+(*jku)) - 1
					for (*jr) = Min((*(n))+(*jkl), (*(m))) - 1; (*jr) <= 1-(*jku); (*jr) += -1 {
						(*extra) = (*zero)
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = (*Sin(angle))
						(*icol) = (MAX(1, (*jr)-(*jkl)+1))
						if (*jr) > 0 {
							(*il) = Min((*(n)), (*jr)+(*jku)+1) + 1 - (*icol)
							Dlarot(&(true), &(false), (*jr)+(*jku) < (*(n)), il, c, s, &((*(a))[(*jr)-(*iskew)*(*icol)+(*ioffst)-1][(*icol)-1]), ilda, dummy, extra)
						}
						//
						//                    Chase "extra" back down
						//
						(*ir) = (*jr)
						for (*jch) = (*jr) + (*jku); (*jch) <= (*iendch); (*jch) += (*jkl) + (*jku) {
							(*ilextr) = (*ir) > 0
							if *ilextr {
								dlartg(&((*(a))[(*ir)-(*iskew)*(*jch)+(*ioffst)-1][(*jch)-1]), extra, c, s, dummy)
							}
							(*ir) = (MAX(1, (*ir)))
							(*irow) = (Min((*(m))-1, (*jch)+(*jkl)))
							(*iltemp) = (*jch)+(*jkl) < (*(m))
							(*temp) = (*zero)
							Dlarot(&(false), ilextr, iltemp, (*irow)+2-(*ir), c, s, &((*(a))[(*ir)-(*iskew)*(*jch)+(*ioffst)-1][(*jch)-1]), ilda, extra, temp)
							if *iltemp {
								dlartg(&((*(a))[(*irow)-(*iskew)*(*jch)+(*ioffst)-1][(*jch)-1]), temp, c, s, dummy)
								(*il) = Min((*iendch), (*jch)+(*jkl)+(*jku)) + 2 - (*jch)
								(*extra) = (*zero)
								Dlarot(&(true), &(true), (*jch)+(*jkl)+(*jku) <= (*iendch), il, c, s, &((*(a))[(*irow)-(*iskew)*(*jch)+(*ioffst)-1][(*jch)-1]), ilda, temp, extra)
								(*ir) = (*irow)
							}
							//Label120:
						}
						//Label130:
					}
					//Label140:
				}
			}
			//
		} else {
			//
			//           Symmetric -- A = U D U'
			//
			(*ipackg) = (*ipack)
			(*ioffg) = (*ioffst)
			//
			if *topdwn {
				//
				//              Top-Down -- Generate Upper triangle only
				//
				if (*ipack) >= 5 {
					(*ipackg) = 6
					(*ioffg) = (*uub) + 1
				} else {
					(*ipackg) = 1
				}
				Dcopy(mnmin, (d), func() *int {y := 1; return &y }(), &((*(a))[1-(*iskew)+(*ioffg)-1][0]), (*ilda)+1)
				//
				for (*k) = 1; (*k) <= (*uub); (*k)++ {
					for (*jc) = 1; (*jc) <= (*(n))-1; (*jc)++ {
						(*irow) = (MAX(1, (*jc)-(*k)))
						(*il) = (Min((*jc)+1, (*k)+2))
						(*extra) = (*zero)
						(*temp) = (*(a))[(*jc)-(*iskew)*((*jc)+1)+(*ioffg)-1][(*jc)+0]
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = (*Sin(angle))
						Dlarot(&(false), (*jc) > (*k), &(true), il, c, s, &((*(a))[(*irow)-(*iskew)*(*jc)+(*ioffg)-1][(*jc)-1]), ilda, extra, temp)
						Dlarot(&(true), &(true), &(false), Min((*k), (*(n))-(*jc))+1, c, s, &((*(a))[(1-(*iskew))*(*jc)+(*ioffg)-1][(*jc)-1]), ilda, temp, dummy)
						//
						//                    Chase extra back up the matrix
						//
						(*icol) = (*jc)
						for (*jch) = (*jc) - (*k); (*jch) <= 1; (*jch) += -(*k) {
							dlartg(&((*(a))[(*jch)+1-(*iskew)*((*icol)+1)+(*ioffg)-1][(*icol)+0]), extra, c, s, dummy)
							(*temp) = (*(a))[(*jch)-(*iskew)*((*jch)+1)+(*ioffg)-1][(*jch)+0]
							Dlarot(&(true), &(true), &(true), (*k)+2, c, -(*s), &((*(a))[(1-(*iskew))*(*jch)+(*ioffg)-1][(*jch)-1]), ilda, temp, extra)
							(*irow) = (MAX(1, (*jch)-(*k)))
							(*il) = (Min((*jch)+1, (*k)+2))
							(*extra) = (*zero)
							Dlarot(&(false), (*jch) > (*k), &(true), il, c, -(*s), &((*(a))[(*irow)-(*iskew)*(*jch)+(*ioffg)-1][(*jch)-1]), ilda, extra, temp)
							(*icol) = (*jch)
							//Label150:
						}
						//Label160:
					}
					//Label170:
				}
				//
				//              If we need lower triangle, copy from upper. Note that
				//              the order of copying is chosen to work for 'q' -> 'b'
				//
				if (*ipack) != (*ipackg) && (*ipack) != 3 {
					for (*jc) = 1; (*jc) <= (*(n)); (*jc)++ {
						(*irow) = (*ioffst) - (*iskew)*(*jc)
						for (*jr) = (*jc); (*jr) <= (Min((*(n)), (*jc)+(*uub))); (*jr)++ {
							(*(a))[(*jr)+(*irow)-1][(*jc)-1] = (*(a))[(*jc)-(*iskew)*(*jr)+(*ioffg)-1][(*jr)-1]
							//Label180:
						}
						//Label190:
					}
					if (*ipack) == 5 {
						for (*jc) = (*(n)) - (*uub) + 1; (*jc) <= (*(n)); (*jc)++ {
							for (*jr) = (*(n)) + 2 - (*jc); (*jr) <= (*uub)+1; (*jr)++ {
								(*(a))[(*jr)-1][(*jc)-1] = (*zero)
								//Label200:
							}
							//Label210:
						}
					}
					if (*ipackg) == 6 {
						(*ipackg) = (*ipack)
					} else {
						(*ipackg) = 0
					}
				}
			} else {
				//
				//              Bottom-Up -- Generate Lower triangle only
				//
				if (*ipack) >= 5 {
					(*ipackg) = 5
					if (*ipack) == 6 {
						(*ioffg) = 1
					}
				} else {
					(*ipackg) = 2
				}
				Dcopy(mnmin, (d), func() *int {y := 1; return &y }(), &((*(a))[1-(*iskew)+(*ioffg)-1][0]), (*ilda)+1)
				//
				for (*k) = 1; (*k) <= (*uub); (*k)++ {
					for (*jc) = (*(n)) - 1; (*jc) <= 1; (*jc) += -1 {
						(*il) = (Min((*(n))+1-(*jc), (*k)+2))
						(*extra) = (*zero)
						(*temp) = (*(a))[1+(1-(*iskew))*(*jc)+(*ioffg)-1][(*jc)-1]
						(*angle) = (*twopi) * Dlarnd(func() *int {y := 1; return &y }(), (iseed))
						(*c) = (*COS(angle))
						(*s) = -Sin(angle)
						Dlarot(&(false), &(true), (*(n))-(*jc) > (*k), il, c, s, &((*(a))[(1-(*iskew))*(*jc)+(*ioffg)-1][(*jc)-1]), ilda, temp, extra)
						(*icol) = (MAX(1, (*jc)-(*k)+1))
						Dlarot(&(true), &(false), &(true), (*jc)+2-(*icol), c, s, &((*(a))[(*jc)-(*iskew)*(*icol)+(*ioffg)-1][(*icol)-1]), ilda, dummy, temp)
						//
						//                    Chase extra back down the matrix
						//
						(*icol) = (*jc)
						for (*jch) = (*jc) + (*k); (*jch) <= (*(n))-1; (*jch) += (*k) {
							dlartg(&((*(a))[(*jch)-(*iskew)*(*icol)+(*ioffg)-1][(*icol)-1]), extra, c, s, dummy)
							(*temp) = (*(a))[1+(1-(*iskew))*(*jch)+(*ioffg)-1][(*jch)-1]
							Dlarot(&(true), &(true), &(true), (*k)+2, c, s, &((*(a))[(*jch)-(*iskew)*(*icol)+(*ioffg)-1][(*icol)-1]), ilda, extra, temp)
							(*il) = (Min((*(n))+1-(*jch), (*k)+2))
							(*extra) = (*zero)
							Dlarot(&(false), &(true), (*(n))-(*jch) > (*k), il, c, s, &((*(a))[(1-(*iskew))*(*jch)+(*ioffg)-1][(*jch)-1]), ilda, temp, extra)
							(*icol) = (*jch)
							//Label220:
						}
						//Label230:
					}
					//Label240:
				}
				//
				//              If we need upper triangle, copy from lower. Note that
				//              the order of copying is chosen to work for 'b' -> 'q'
				//
				if (*ipack) != (*ipackg) && (*ipack) != 4 {
					for (*jc) = (*(n)); (*jc) <= 1; (*jc) += -1 {
						(*irow) = (*ioffst) - (*iskew)*(*jc)
						for (*jr) = (*jc); (*jr) <= (MAX(1, (*jc)-(*uub))); (*jr) += -1 {
							(*(a))[(*jr)+(*irow)-1][(*jc)-1] = (*(a))[(*jc)-(*iskew)*(*jr)+(*ioffg)-1][(*jr)-1]
							//Label250:
						}
						//Label260:
					}
					if (*ipack) == 6 {
						for (*jc) = 1; (*jc) <= (*uub); (*jc)++ {
							for (*jr) = 1; (*jr) <= (*uub)+1-(*jc); (*jr)++ {
								(*(a))[(*jr)-1][(*jc)-1] = (*zero)
								//Label270:
							}
							//Label280:
						}
					}
					if (*ipackg) == 5 {
						(*ipackg) = (*ipack)
					} else {
						(*ipackg) = 0
					}
				}
			}
		}
		//
	} else {
		//
		//        4)      Generate Banded Matrix by first
		//                Rotating by random Unitary matrices,
		//                then reducing the bandwidth using Householder
		//                transformations.
		//
		//                Note: we should get here only if lda .ge. N
		//
		if (*isym) == 1 {
			//
			//           Non-symmetric -- A = U D V
			//
			Dlagge(MR, nc, llb, uub, (d), (a), (lda), (iseed), (work), iinfo)
		} else {
			//
			//           Symmetric -- A = U D U'
			//
			Dlagsy((m), llb, (d), (a), (lda), (iseed), (work), iinfo)
			//
		}
		if (*iinfo) != 0 {
			(*(info)) = 3
			return
		}
	}
	//
	//     5)      Pack the matrix
	//
	if (*ipack) != (*ipackg) {
		if (*ipack) == 1 {
			//
			//           'U' -- Upper triangular, not packed
			//
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				for (*i) = (*j) + 1; (*i) <= (*(m)); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label290:
				}
				//Label300:
			}
			//
		} else if (*ipack) == 2 {
			//
			//           'L' -- Lower triangular, not packed
			//
			for (*j) = 2; (*j) <= (*(m)); (*j)++ {
				for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*zero)
					//Label310:
				}
				//Label320:
			}
			//
		} else if (*ipack) == 3 {
			//
			//           'C' -- Upper triangle packed Columnwise.
			//
			(*icol) = 1
			(*irow) = 0
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*irow) = (*irow) + 1
					if (*irow) > (*(lda)) {
						(*irow) = 1
						(*icol) = (*icol) + 1
					}
					(*(a))[(*irow)-1][(*icol)-1] = (*(a))[(*i)-1][(*j)-1]
					//Label330:
				}
				//Label340:
			}
			//
		} else if (*ipack) == 4 {
			//
			//           'R' -- Lower triangle packed Columnwise.
			//
			(*icol) = 1
			(*irow) = 0
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				for (*i) = (*j); (*i) <= (*(m)); (*i)++ {
					(*irow) = (*irow) + 1
					if (*irow) > (*(lda)) {
						(*irow) = 1
						(*icol) = (*icol) + 1
					}
					(*(a))[(*irow)-1][(*icol)-1] = (*(a))[(*i)-1][(*j)-1]
					//Label350:
				}
				//Label360:
			}
			//
		} else if (*ipack) >= 5 {
			//
			//           'B' -- The lower triangle is packed as a band matrix.
			//           'Q' -- The upper triangle is packed as a band matrix.
			//           'Z' -- The whole matrix is packed as a band matrix.
			//
			if (*ipack) == 5 {
				(*uub) = 0
			}
			if (*ipack) == 6 {
				(*llb) = 0
			}
			//
			for (*j) = 1; (*j) <= (*uub); (*j)++ {
				for (*i) = Min((*j)+(*llb), (*(m))); (*i) <= 1; (*i) += -1 {
					(*(a))[(*i)-(*j)+(*uub)+0][(*j)-1] = (*(a))[(*i)-1][(*j)-1]
					//Label370:
				}
				//Label380:
			}
			//
			for (*j) = (*uub) + 2; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) - (*uub); (*i) <= (Min((*j)+(*llb), (*(m)))); (*i)++ {
					(*(a))[(*i)-(*j)+(*uub)+0][(*j)-1] = (*(a))[(*i)-1][(*j)-1]
					//Label390:
				}
				//Label400:
			}
		}
		//
		//        If packed, zero out extraneous elements.
		//
		//        Symmetric/Triangular Packed --
		//        zero out everything after a(irow,ICOL)
		//
		if (*ipack) == 3 || (*ipack) == 4 {
			for (*jc) = (*icol); (*jc) <= (*(m)); (*jc)++ {
				for (*jr) = (*irow) + 1; (*jr) <= (*(lda)); (*jr)++ {
					(*(a))[(*jr)-1][(*jc)-1] = (*zero)
					//Label410:
				}
				(*irow) = 0
				//Label420:
			}
			//
		} else if (*ipack) >= 5 {
			//
			//           Packed Band --
			//              1st row is now in a( uub+2-j, j), zero above it
			//              m-th row is now in a( M+uub-j,j), zero below it
			//              last non-zero diagonal is now in a( uub+llb+1,j),
			//                 zero below it, too.
			//
			(*ir1) = (*uub) + (*llb) + 2
			(*ir2) = (*uub) + (*(m)) + 2
			for (*jc) = 1; (*jc) <= (*(n)); (*jc)++ {
				for (*jr) = 1; (*jr) <= (*uub)+1-(*jc); (*jr)++ {
					(*(a))[(*jr)-1][(*jc)-1] = (*zero)
					//Label430:
				}
				for (*jr) = MAX(1, Min((*ir1), (*ir2)-(*jc))); (*jr) <= (*(lda)); (*jr)++ {
					(*(a))[(*jr)-1][(*jc)-1] = (*zero)
					//Label440:
				}
				//Label450:
			}
		}
	}
	//
	return
	//
	//     End of Dlatms
	//
}
