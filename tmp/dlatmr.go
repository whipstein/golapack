package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// \brief \b DLAtmR
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE DLATmr( m, n, dist, iseed, sym, d, mode, cond, dmax,
//                          rsign, GRADE, DL, modeL, condL, DR, modeR,
//                          condR, pivTNG, ipivOT, kl, ku, sparse, anorm,
//                          PACK, a, lda, iwork, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          dist, GRADE, PACK, pivTNG, rsign, sym
//       intEGER            info, kl, ku, lda, m, mode, modeL, modeR, N
//       DOUBLE PRECISION   anorm, cond, condL, condR, dmax, sparse
//       ..
//       .. Array Arguments ..
//       intEGER            ipivOt(*), iseed( 4), iwork(*)
//       DOUBLE PRECISION   a( lda, *), d(*), DL(*), dr(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    DLAtmR generates random matrices of various types for testing
//    lapACK programs.
//
//    DLAtmR operates by applying the following sequence of
//    operations:
//
//      Generate a matrix A with random entries of distribution dist
//         which is symmetric if sym='S', and nonsymmetric
//         if sym='N'.
//
//      Set the diagonal to d, where D may be input or
//         computed according to mode, cond, dmax and rsign
//         as described below.
//
//      Grade the matrix, if desired, from the left and/or right
//         as specified by GRADE. The inputs DL, modeL, condL, DR,
//         modeR and condR also determine the grading as described
//         below.
//
//      Permute, if desired, the rows and/or columns as specified by
//         pivTNG and ipivOT.
//
//      Set random entries to zero, if desired, to get a random sparse
//         matrix as specified by sparse.
//
//      Make A a band matrix, if desired, by zeroing out the matrix
//         outside a band of lower bandwidth kl and upper bandwidth ku.
//
//      Scale a, if desired, to have maximum entry anorm.
//
//      Pack the matrix if desired. Options specified by PACK are:
//         no packing
//         zero out upper half (if symmetric)
//         zero out lower half (if symmetric)
//         store the upper half columnwise (if symmetric or
//             square upper triangular)
//         store the lower half columnwise (if symmetric or
//             square lower triangular)
//             same as upper half rowwise if symmetric
//         store the lower triangle in banded format (if symmetric)
//         store the upper triangle in banded format (if symmetric)
//         store the entire matrix in banded format
//
//    Note: If two calls to DLAtmR differ only in the PACK parameter,
//          they will generate mathematically equivalent matrices.
//
//          If two calls to DLAtmR both have full bandwidth (kl = M-1
//          and ku = N-1), and differ only in the pivTNG and PACK
//          parameters, then the matrices generated will differ only
//          in the order of the rows and/or columns, and otherwise
//          contain the same data. This consistency cannot be and
//          is not maintained with less than full bandwidth.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//           Number of rows of A. Not modified.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//           Number of columns of A. Not modified.
// \endverbatim
//
// \param[in] dist
// \verbatim
//          dist is CHARACTER*1
//           On entry, dist specifies the type of distribution to be used
//           to generate a random matrix .
//           'U' => UNIFORM( 0, 1)  ( 'U' for uniform)
//           'S' => UNIFORM( -1, 1) ( 'S' for symmetric)
//           'N' => normaL( 0, 1)   ( 'N' for normal)
//           Not modified.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is intEGER array, dimension (4)
//           On entry iseed specifies the seed of the random number
//           generator. They should lie between 0 and 4095 inclusive,
//           and iseed(4) should be odd. The random number generator
//           uses a linear congruential sequence limited to small
//           integers, and so should produce machine independent
//           random numbers. The values of iseed are changed on
//           exit, and can be used in the next call to DLAtmR
//           to continue the same random number sequence.
//           Changed on exit.
// \endverbatim
//
// \param[in] sym
// \verbatim
//          sym is CHARACTER*1
//           If sym='S' or 'H', generated matrix is symmetric.
//           If sym='N', generated matrix is nonsymmetric.
//           Not modified.
// \endverbatim
//
// \param[in,out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (min(m,N))
//           On entry this array specifies the diagonal entries
//           of the diagonal of A.  D may either be specified
//           on entry, or set according to mode and cond as described
//           below. May be changed on exit if mode is nonzero.
// \endverbatim
//
// \param[in] mode
// \verbatim
//          mode is intEGER
//           On entry describes how D is to be used:
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
// \param[in] dmax
// \verbatim
//          dmax is DOUBLE PRECISION
//           If mode neither -6, 0 nor 6, the diagonal is scaled by
//           dmax / max(abs(d(i))), so that maximum absolute entry
//           of diagonal is ABS(dmax). If dmax is negative (or zero),
//           diagonal will be scaled by a negative number (or zero).
// \endverbatim
//
// \param[in] rsign
// \verbatim
//          rsign is CHARACTER*1
//           If mode neither -6, 0 nor 6, specifies sign of diagonal
//           as follows:
//           'T' => diagonal entries are multiplied by 1 or -1
//                  with probability .5
//           'F' => diagonal unchanged
//           Not modified.
// \endverbatim
//
// \param[in] GRADE
// \verbatim
//          GRADE is CHARACTER*1
//           Specifies grading of matrix as follows:
//           'N'  => no grading
//           'L'  => matrix premultiplied by diag( DL)
//                   (only if matrix nonsymmetric)
//           'R'  => matrix postmultiplied by diag( DR)
//                   (only if matrix nonsymmetric)
//           'B'  => matrix premultiplied by diag( DL) and
//                         postmultiplied by diag( DR)
//                   (only if matrix nonsymmetric)
//           'S' or 'H'  => matrix premultiplied by diag( DL) and
//                          postmultiplied by diag( DL)
//                          ('S' for symmetric, or 'H' for Hermitian)
//           'E'  => matrix premultiplied by diag( DL) and
//                         postmultiplied by inv( diag( DL))
//                         ( 'E' for eigenvalue invariance)
//                   (only if matrix nonsymmetric)
//                   Note: if GRADE='E', then M must equal N.
//           Not modified.
// \endverbatim
//
// \param[in,out] DL
// \verbatim
//          DL is DOUBLE PRECISION array, dimension (m)
//           If modeL=0, then on entry this array specifies the diagonal
//           entries of a diagonal matrix used as described under GRADE
//           above. If modeL is not zero, then DL will be set according
//           to modeL and condL, analogous to the way D is set according
//           to mode and cond (except there is no dmax parameter for DL).
//           If GRADE='E', then DL cannot have zero entries.
//           Not referenced if GRADE = 'N' or 'R'. Changed on exit.
// \endverbatim
//
// \param[in] modeL
// \verbatim
//          modeL is intEGER
//           This specifies how the diagonal array DL is to be computed,
//           just as mode specifies how D is to be computed.
//           Not modified.
// \endverbatim
//
// \param[in] condL
// \verbatim
//          condL is DOUBLE PRECISION
//           When modeL is not zero, this specifies the condition number
//           of the computed DL.  Not modified.
// \endverbatim
//
// \param[in,out] DR
// \verbatim
//          DR is DOUBLE PRECISION array, dimension (n)
//           If modeR=0, then on entry this array specifies the diagonal
//           entries of a diagonal matrix used as described under GRADE
//           above. If modeR is not zero, then DR will be set according
//           to modeR and condR, analogous to the way D is set according
//           to mode and cond (except there is no dmax parameter for DR).
//           Not referenced if GRADE = 'N', 'L', 'H', 'S' or 'E'.
//           Changed on exit.
// \endverbatim
//
// \param[in] modeR
// \verbatim
//          modeR is intEGER
//           This specifies how the diagonal array DR is to be computed,
//           just as mode specifies how D is to be computed.
//           Not modified.
// \endverbatim
//
// \param[in] condR
// \verbatim
//          condR is DOUBLE PRECISION
//           When modeR is not zero, this specifies the condition number
//           of the computed DR.  Not modified.
// \endverbatim
//
// \param[in] pivTNG
// \verbatim
//          pivTNG is CHARACTER*1
//           On entry specifies pivoting permutations as follows:
//           'N' or ' ' => none.
//           'L' => left or row pivoting (matrix must be nonsymmetric).
//           'R' => right or column pivoting (matrix must be
//                  nonsymmetric).
//           'B' or 'F' => both or full pivoting, i.e., on both sides.
//                         In this case, M must equal N
//
//           If two calls to DLAtmR both have full bandwidth (kl = M-1
//           and ku = N-1), and differ only in the pivTNG and PACK
//           parameters, then the matrices generated will differ only
//           in the order of the rows and/or columns, and otherwise
//           contain the same data. This consistency cannot be
//           maintained with less than full bandwidth.
// \endverbatim
//
// \param[in] ipivOT
// \verbatim
//          ipivOT is intEGER array, dimension (N or M)
//           This array specifies the permutation used.  After the
//           basic matrix is generated, the rows, columns, or both
//           are permuted.   If, say, row pivoting is selected, DLAtmR
//           starts with the *last* row and interchanges the M-th and
//           ipivOt(m)-th rows, then moves to the next-to-last row,
//           interchanging the (M-1)-th and the ipivOt(M-1)-th rows,
//           and so on.  In terms of "2-cycles", the permutation is
//           (1 ipivOt1) (2 ipivOt(2)) ... (M ipivOt(m))
//           where the rightmost cycle is applied first.  This is the
//           *inverse* of the effect of pivoting in LinPACK.  The idea
//           is that factoring (with pivoting) an identity matrix
//           which has been inverse-pivoted in this way should
//           result in a pivot vector identical to ipivOT.
//           Not referenced if pivTNG = 'N'. Not modified.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is intEGER
//           On entry specifies the lower bandwidth of the  matrix. For
//           example, kl=0 implies upper triangular, kl=1 implies upper
//           Hessenberg, and kl at least M-1 implies the matrix is not
//           banded. Must equal ku if matrix is symmetric.
//           Not modified.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//           On entry specifies the upper bandwidth of the  matrix. For
//           example, ku=0 implies lower triangular, ku=1 implies lower
//           Hessenberg, and ku at least N-1 implies the matrix is not
//           banded. Must equal kl if matrix is symmetric.
//           Not modified.
// \endverbatim
//
// \param[in] sparse
// \verbatim
//          sparse is DOUBLE PRECISION
//           On entry specifies the sparsity of the matrix if a sparse
//           matrix is to be generated. sparse should lie between
//           0 and 1. To generate a sparse matrix, for each matrix entry
//           a uniform ( 0, 1) random number x is generated and
//           compared to sparse; if x is larger the matrix entry
//           is unchanged and if x is smaller the entry is set
//           to zero. Thus on the average a fraction sparse of the
//           entries will be set to zero.
//           Not modified.
// \endverbatim
//
// \param[in] anorm
// \verbatim
//          anorm is DOUBLE PRECISION
//           On entry specifies maximum entry of output matrix
//           (output matrix will by multiplied by a _constant so that
//           its largest absolute entry equal anorm)
//           if anorm is nonnegative. If anorm is negative no scaling
//           is done. Not modified.
// \endverbatim
//
// \param[in] PACK
// \verbatim
//          PACK is CHARACTER*1
//           On entry specifies packing of matrix as follows:
//           'N' => no packing
//           'U' => zero out all subdiagonal entries (if symmetric)
//           'L' => zero out all superdiagonal entries (if symmetric)
//           'C' => store the upper triangle columnwise
//                  (only if matrix symmetric or square upper triangular)
//           'R' => store the lower triangle columnwise
//                  (only if matrix symmetric or square lower triangular)
//                  (same as upper half rowwise if symmetric)
//           'B' => store the lower triangle in band storage scheme
//                  (only if matrix symmetric)
//           'Q' => store the upper triangle in band storage scheme
//                  (only if matrix symmetric)
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
//           If two calls to DLAtmR differ only in the PACK parameter,
//           they will generate mathematically equivalent matrices.
//           Not modified.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//           On exit A is the desired test matrix. Only those
//           entries of A which are significant on output
//           will be referenced (even if A is in packed or band
//           storage format). The 'unoccupied corners' of A in
//           band format will be zeroed out.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//           on entry lda specifies the first dimension of A as
//           declared in the calling program.
//           If PACK='N', 'U' or 'L', lda must be at least max ( 1, M).
//           If PACK='C' or 'R', lda must be at least 1.
//           If PACK='B', or 'Q', lda must be Min ( ku+1, N)
//           If PACK='Z', lda must be at least kuu+klL+1, where
//           kuu = Min ( ku, N-1) and klL = Min ( kl, M-1)
//           Not modified.
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is intEGER array, dimension ( N or M)
//           workspace. Not referenced if pivTNG = 'N'. Changed on exit.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
//           Error parameter on exit:
//             0 => normal return
//            -1 => M negative or unequal to N and sym='S' or 'H'
//            -2 => N negative
//            -3 => dist illegal string
//            -5 => sym illegal string
//            -7 => mode not in range -6 to 6
//            -8 => cond less than 1.0, and mode neither -6, 0 nor 6
//           -10 => mode neither -6, 0 nor 6 and rsign illegal string
//           -11 => GRADE illegal string, or GRADE='E' and
//                  M not equal to n, or GRADE='L', 'R', 'B' or 'E' and
//                  sym = 'S' or 'H'
//           -12 => GRADE = 'E' and DL contains zero
//           -13 => modeL not in range -6 to 6 and GRADE= 'L', 'B', 'H',
//                  'S' or 'E'
//           -14 => condL less than 1.0, GRADE='L', 'B', 'H', 'S' or 'E',
//                  and modeL neither -6, 0 nor 6
//           -16 => modeR not in range -6 to 6 and GRADE= 'R' or 'B'
//           -17 => condR less than 1.0, GRADE='R' or 'B', and
//                  modeR neither -6, 0 nor 6
//           -18 => pivTNG illegal string, or pivTNG='B' or 'F' and
//                  M not equal to n, or pivTNG='L' or 'R' and sym='S'
//                  or 'H'
//           -19 => ipivOT contains out of range number and
//                  pivTNG not equal to 'N'
//           -20 => kl negative
//           -21 => ku negative, or sym='S' or 'H' and ku not equal to kl
//           -22 => sparse not in range 0. to 1.
//           -24 => PACK illegal string, or PACK='U', 'L', 'B' or 'Q'
//                  and sym='N', or PACK='C' and sym='N' and either kl
//                  not equal to 0 or N not equal to m, or PACK='R' and
//                  sym='N', and either ku not equal to 0 or N not equal
//                  to M
//           -26 => lda too small
//             1 => Error return from Dlatm1 (computing D)
//             2 => Cannot scale diagonal to dmax (max. entry is 0)
//             3 => Error return from Dlatm1 (computing DL)
//             4 => Error return from Dlatm1 (computing DR)
//             5 => anorm is positive, but matrix _constructed prior to
//                  attempting to scale it to have norm anorm, is zero
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
func DLATmr(m *int, n *int, dist *byte, iseed *[]int, sym *byte, d *[]float64, mode *int, cond *float64, dmax *float64, rsign *byte, GRADE *byte, dl *[]float64, modeL *int, condL *float64, dr *[]float64, moder *int, condr *float64, pivTNG *byte, ipivOT *[]int, kl *int, ku *int, sparse *float64, anorm *float64, pack *byte, a *[][]float64, lda *int, iwork *[]int, info *int) {
	zero := new(float64)
	one := new(float64)
	BADPVt := new(bool)
	Dzero := new(bool)
	FULBNd := new(bool)
	i := new(int)
	idist := new(int)
	IGRAde := new(int)
	Iisub := new(int)
	ipack := new(int)
	ipvtng := new(int)
	irsign := new(int)
	isub := new(int)
	isym := new(int)
	j := new(int)
	Jjsub := new(int)
	jsub := new(int)
	k := new(int)
	kll := new(int)
	kuu := new(int)
	mnmin := new(int)
	mnSUb := new(int)
	MXSUb := new(int)
	NPVts := new(int)
	alpha := new(float64)
	Onorm := new(float64)
	temp := new(float64)
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
	if blas.Lsame((sym), func() *byte {y := byte('S'); return &y }()) {
		(*isym) = 0
	} else if blas.Lsame((sym), func() *byte {y := byte('N'); return &y }()) {
		(*isym) = 1
	} else if blas.Lsame((sym), func() *byte {y := byte('H'); return &y }()) {
		(*isym) = 0
	} else {
		(*isym) = -1
	}
	//
	//     Decode rsign
	//
	if blas.Lsame((rsign), func() *byte {y := byte('F'); return &y }()) {
		(*irsign) = 0
	} else if blas.Lsame((rsign), func() *byte {y := byte('T'); return &y }()) {
		(*irsign) = 1
	} else {
		(*irsign) = -1
	}
	//
	//     Decode pivTNG
	//
	if blas.Lsame((pivTNG), func() *byte {y := byte('N'); return &y }()) {
		(*iPVTNG) = 0
	} else if blas.Lsame((pivTNG), func() *byte {y := byte(' '); return &y }()) {
		(*iPVTNG) = 0
	} else if blas.Lsame((pivTNG), func() *byte {y := byte('L'); return &y }()) {
		(*iPVTNG) = 1
		(*nPVts) = (*(m))
	} else if blas.Lsame((pivTNG), func() *byte {y := byte('R'); return &y }()) {
		(*iPVTNG) = 2
		(*nPVts) = (*(n))
	} else if blas.Lsame((pivTNG), func() *byte {y := byte('B'); return &y }()) {
		(*iPVTNG) = 3
		(*nPVts) = (Min((*(n)), (*(m))))
	} else if blas.Lsame((pivTNG), func() *byte {y := byte('F'); return &y }()) {
		(*iPVTNG) = 3
		(*nPVts) = (Min((*(n)), (*(m))))
	} else {
		(*iPVTNG) = -1
	}
	//
	//     Decode GRADE
	//
	if blas.Lsame((GRADE), func() *byte {y := byte('N'); return &y }()) {
		(*iGRADE) = 0
	} else if blas.Lsame((GRADE), func() *byte {y := byte('L'); return &y }()) {
		(*iGRADE) = 1
	} else if blas.Lsame((GRADE), func() *byte {y := byte('R'); return &y }()) {
		(*iGRADE) = 2
	} else if blas.Lsame((GRADE), func() *byte {y := byte('B'); return &y }()) {
		(*iGRADE) = 3
	} else if blas.Lsame((GRADE), func() *byte {y := byte('E'); return &y }()) {
		(*iGRADE) = 4
	} else if (*blas.Lsame((GRADE), func() *byte {y := byte('H'); return &y }())) || (*blas.Lsame((GRADE), func() *byte {y := byte('S'); return &y }())) {
		(*iGRADE) = 5
	} else {
		(*iGRADE) = -1
	}
	//
	//     Decode PACK
	//
	if blas.Lsame((pack), func() *byte {y := byte('N'); return &y }()) {
		(*ipack) = 0
	} else if blas.Lsame((pack), func() *byte {y := byte('U'); return &y }()) {
		(*ipack) = 1
	} else if blas.Lsame((pack), func() *byte {y := byte('L'); return &y }()) {
		(*ipack) = 2
	} else if blas.Lsame((pack), func() *byte {y := byte('C'); return &y }()) {
		(*ipack) = 3
	} else if blas.Lsame((pack), func() *byte {y := byte('R'); return &y }()) {
		(*ipack) = 4
	} else if blas.Lsame((pack), func() *byte {y := byte('B'); return &y }()) {
		(*ipack) = 5
	} else if blas.Lsame((pack), func() *byte {y := byte('Q'); return &y }()) {
		(*ipack) = 6
	} else if blas.Lsame((pack), func() *byte {y := byte('Z'); return &y }()) {
		(*ipack) = 7
	} else {
		(*ipack) = -1
	}
	//
	//     Set certain internal parameters
	//
	(*mnmin) = (Min((*(m)), (*(n))))
	(*klL) = (Min((*(kl)), (*(m))-1))
	(*kuu) = (Min((*(ku)), (*(n))-1))
	//
	//     If inv(dl) is used, check to see if DL has a zero entry.
	//
	(*Dzero) = false
	if (*iGRADE) == 4 && (*(modeL)) == 0 {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			if (*(dl))[(*i)-1] == (*zero) {
				(*Dzero) = true
			}
			//Label10:
		}
	}
	//
	//     Check values in ipivOT
	//
	(*BADPVT) = false
	if (*iPVTNG) > 0 {
		for (*j) = 1; (*j) <= (*nPVts); (*j)++ {
			if (*(ipivOT))[(*j)-1] <= 0 || (*(ipivOT))[(*j)-1] > (*nPVts) {
				(*BADPVT) = true
			}
			//Label20:
		}
	}
	//
	//     Set info if an error
	//
	if (*(m)) < 0 {
		(*(info)) = -1
	} else if (*(m)) != (*(n)) && (*isym) == 0 {
		(*(info)) = -1
	} else if (*(n)) < 0 {
		(*(info)) = -2
	} else if (*idist) == -1 {
		(*(info)) = -3
	} else if (*isym) == -1 {
		(*(info)) = -5
	} else if (*(mode)) < -6 || (*(mode)) > 6 {
		(*(info)) = -7
	} else if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && (*(cond)) < (*one) {
		(*(info)) = -8
	} else if ((*(mode)) != -6 && (*(mode)) != 0 && (*(mode)) != 6) && (*irsign) == -1 {
		(*(info)) = -10
	} else if (*iGRADE) == -1 || ((*iGRADE) == 4 && (*(m)) != (*(n))) || (((*iGRADE) >= 1 && (*iGRADE) <= 4) && (*isym) == 0) {
		(*(info)) = -11
	} else if (*iGRADE) == 4 && (*Dzero) {
		(*(info)) = -12
	} else if ((*iGRADE) == 1 || (*iGRADE) == 3 || (*iGRADE) == 4 || (*iGRADE) == 5) && ((*(modeL)) < -6 || (*(modeL)) > 6) {
		(*(info)) = -13
	} else if ((*iGRADE) == 1 || (*iGRADE) == 3 || (*iGRADE) == 4 || (*iGRADE) == 5) && ((*(modeL)) != -6 && (*(modeL)) != 0 && (*(modeL)) != 6) && (*(condL)) < (*one) {
		(*(info)) = -14
	} else if ((*iGRADE) == 2 || (*iGRADE) == 3) && ((*(modeR)) < -6 || (*(modeR)) > 6) {
		(*(info)) = -16
	} else if ((*iGRADE) == 2 || (*iGRADE) == 3) && ((*(modeR)) != -6 && (*(modeR)) != 0 && (*(modeR)) != 6) && (*(condR)) < (*one) {
		(*(info)) = -17
	} else if (*iPVTNG) == -1 || ((*iPVTNG) == 3 && (*(m)) != (*(n))) || (((*iPVTNG) == 1 || (*iPVTNG) == 2) && (*isym) == 0) {
		(*(info)) = -18
	} else if (*iPVTNG) != 0 && (*BADPVT) {
		(*(info)) = -19
	} else if (*(kl)) < 0 {
		(*(info)) = -20
	} else if (*(ku)) < 0 || ((*isym) == 0 && (*(kl)) != (*(ku))) {
		(*(info)) = -21
	} else if (*(sparse)) < (*zero) || (*(sparse)) > (*one) {
		(*(info)) = -22
	} else if (*ipack) == -1 || (((*ipack) == 1 || (*ipack) == 2 || (*ipack) == 5 || (*ipack) == 6) && (*isym) == 1) || ((*ipack) == 3 && (*isym) == 1 && ((*(kl)) != 0 || (*(m)) != (*(n)))) || ((*ipack) == 4 && (*isym) == 1 && ((*(ku)) != 0 || (*(m)) != (*(n)))) {
		(*(info)) = -24
	} else if (((*ipack) == 0 || (*ipack) == 1 || (*ipack) == 2) && (*(lda)) < MAX(1, (*(m)))) || (((*ipack) == 3 || (*ipack) == 4) && (*(lda)) < 1) || (((*ipack) == 5 || (*ipack) == 6) && (*(lda)) < (*kuu)+1) || ((*ipack) == 7 && (*(lda)) < (*klL)+(*kuu)+1) {
		(*(info)) = -26
	}
	//
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("DLAtmR"); return &y }(), -(*(info)))
		return
	}
	//
	//     Decide if we can pivot consistently
	//
	(*FULBND) = false
	if (*kuu) == (*(n))-1 && (*klL) == (*(m))-1 {
		(*FULBND) = true
	}
	//
	//     Initialize random number generator
	//
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*(iseed))[(*i)-1] = (MOD(ABS(((*(iseed))[(*i)-1])), int(4096)))
		//Label30:
	}
	//
	(*(iseed))[3] = 2*((*(iseed))[3]/2) + 1
	//
	//     2)      Set up d, DL, and DR, if indicated.
	//
	//             Compute D according to cond and mode
	//
	Dlatm1((mode), (cond), irsign, idist, (iseed), (d), mnmin, (info))
	if (*(info)) != 0 {
		(*(info)) = 1
		return
	}
	if (*(mode)) != 0 && (*(mode)) != -6 && (*(mode)) != 6 {
		//
		//        Scale by dmax
		//
		(*temp) = (ABS(((*(d))[0])))
		for (*i) = 2; (*i) <= (*mnmin); (*i)++ {
			(*temp) = (MAX((*temp), ABS(((*(d))[(*i)-1]))))
			//Label40:
		}
		if (*temp) == (*zero) && (*(dmax)) != (*zero) {
			(*(info)) = 2
			return
		}
		if (*temp) != (*zero) {
			(*alpha) = (*(dmax)) / (*temp)
		} else {
			(*alpha) = (*one)
		}
		for (*i) = 1; (*i) <= (*mnmin); (*i)++ {
			(*(d))[(*i)-1] = (*alpha) * (*(d))[(*i)-1]
			//Label50:
		}
		//
	}
	//
	//     Compute DL if grading set
	//
	if (*iGRADE) == 1 || (*iGRADE) == 3 || (*iGRADE) == 4 || (*iGRADE) == 5 {
		Dlatm1((modeL), (condL), func() *int {y := 0; return &y }(), idist, (iseed), (dl), (m), (info))
		if (*(info)) != 0 {
			(*(info)) = 3
			return
		}
	}
	//
	//     Compute DR if grading set
	//
	if (*iGRADE) == 2 || (*iGRADE) == 3 {
		Dlatm1((modeR), (condR), func() *int {y := 0; return &y }(), idist, (iseed), (dr), (n), (info))
		if (*(info)) != 0 {
			(*(info)) = 4
			return
		}
	}
	//
	//     3)     Generate iwork if pivoting
	//
	if (*iPVTNG) > 0 {
		for (*i) = 1; (*i) <= (*nPVts); (*i)++ {
			(*(iwork))[(*i)-1] = (*i)
			//Label60:
		}
		if *FULBND {
			for (*i) = 1; (*i) <= (*nPVts); (*i)++ {
				(*k) = (*(ipivOT))[(*i)-1]
				(*j) = (*(iwork))[(*i)-1]
				(*(iwork))[(*i)-1] = (*(iwork))[(*k)-1]
				(*(iwork))[(*k)-1] = (*j)
				//Label70:
			}
		} else {
			for (*i) = (*nPVts); (*i) <= 1; (*i) += -1 {
				(*k) = (*(ipivOT))[(*i)-1]
				(*j) = (*(iwork))[(*i)-1]
				(*(iwork))[(*i)-1] = (*(iwork))[(*k)-1]
				(*(iwork))[(*k)-1] = (*j)
				//Label80:
			}
		}
	}
	//
	//     4)      Generate matrices for each kind of PACKing
	//             Always sweep matrix columnwise (if symmetric, upper
	//             half only) so that matrix generated does not depend
	//             on PACK
	//
	if *FULBND {
		//
		//        Use dlatm3 so matrices generated with differing pivOting only
		//        differ only in the order of their rows and/or columns.
		//
		if (*ipack) == 0 {
			if (*isym) == 0 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = 1; (*i) <= (*j); (*i)++ {
						(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*(a))[(*isub)-1][(*jsub)-1] = (*temp)
						(*(a))[(*jsub)-1][(*isub)-1] = (*temp)
						//Label90:
					}
					//Label100:
				}
			} else if (*isym) == 1 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = 1; (*i) <= (*(m)); (*i)++ {
						(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*(a))[(*isub)-1][(*jsub)-1] = (*temp)
						//Label110:
					}
					//Label120:
				}
			}
			//
		} else if (*ipack) == 1 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					(*mnSUB) = (Min((*isub), (*jsub)))
					(*mXSUB) = (MAX((*isub), (*jsub)))
					(*(a))[(*mnSUB)-1][(*mXSUB)-1] = (*temp)
					if (*mnSUB) != (*mXSUB) {
						(*(a))[(*mXSUB)-1][(*mnSUB)-1] = (*zero)
					}
					//Label130:
				}
				//Label140:
			}
			//
		} else if (*ipack) == 2 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					(*mnSUB) = (Min((*isub), (*jsub)))
					(*mXSUB) = (MAX((*isub), (*jsub)))
					(*(a))[(*mXSUB)-1][(*mnSUB)-1] = (*temp)
					if (*mnSUB) != (*mXSUB) {
						(*(a))[(*mnSUB)-1][(*mXSUB)-1] = (*zero)
					}
					//Label150:
				}
				//Label160:
			}
			//
		} else if (*ipack) == 3 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					//
					//                 Compute K = location of (isub,jsub) entry in packed
					//                 array
					//
					(*mnSUB) = (Min((*isub), (*jsub)))
					(*mXSUB) = (MAX((*isub), (*jsub)))
					(*k) = (*mXSUB)*((*mXSUB)-1)/2 + (*mnSUB)
					//
					//                 Convert K to (Iisub,Jjsub) location
					//
					(*jjsub) = ((*k)-1)/(*(lda)) + 1
					(*iisub) = (*k) - (*(lda))*((*jjsub)-1)
					//
					(*(a))[(*iisub)-1][(*jjsub)-1] = (*temp)
					//Label170:
				}
				//Label180:
			}
			//
		} else if (*ipack) == 4 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					//
					//                 Compute K = location of (I,j) entry in packed array
					//
					(*mnSUB) = (Min((*isub), (*jsub)))
					(*mXSUB) = (MAX((*isub), (*jsub)))
					if (*mnSUB) == 1 {
						(*k) = (*mXSUB)
					} else {
						(*k) = (*(n))*((*(n))+1)/2 - ((*(n))-(*mnSUB)+1)*((*(n))-(*mnSUB)+2)/2 + (*mXSUB) - (*mnSUB) + 1
					}
					//
					//                 Convert K to (Iisub,Jjsub) location
					//
					(*jjsub) = ((*k)-1)/(*(lda)) + 1
					(*iisub) = (*k) - (*(lda))*((*jjsub)-1)
					//
					(*(a))[(*iisub)-1][(*jjsub)-1] = (*temp)
					//Label190:
				}
				//Label200:
			}
			//
		} else if (*ipack) == 5 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
					if (*i) < 1 {
						(*(a))[(*j)-(*i)+0][(*i)+(*(n))-1] = (*zero)
					} else {
						(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*mnSUB) = (Min((*isub), (*jsub)))
						(*mXSUB) = (MAX((*isub), (*jsub)))
						(*(a))[(*mXSUB)-(*mnSUB)+0][(*mnSUB)-1] = (*temp)
					}
					//Label210:
				}
				//Label220:
			}
			//
		} else if (*ipack) == 6 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
					(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					(*mnSUB) = (Min((*isub), (*jsub)))
					(*mXSUB) = (MAX((*isub), (*jsub)))
					(*(a))[(*mnSUB)-(*mXSUB)+(*kuu)+0][(*mXSUB)-1] = (*temp)
					//Label230:
				}
				//Label240:
			}
			//
		} else if (*ipack) == 7 {
			//
			if (*isym) == 0 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
						(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*mnSUB) = (Min((*isub), (*jsub)))
						(*mXSUB) = (MAX((*isub), (*jsub)))
						(*(a))[(*mnSUB)-(*mXSUB)+(*kuu)+0][(*mXSUB)-1] = (*temp)
						if (*i) < 1 {
							(*(a))[(*j)-(*i)+1+(*kuu)-1][(*i)+(*(n))-1] = (*zero)
						}
						if (*i) >= 1 && (*mnSUB) != (*mXSUB) {
							(*(a))[(*mXSUB)-(*mnSUB)+1+(*kuu)-1][(*mnSUB)-1] = (*temp)
						}
						//Label250:
					}
					//Label260:
				}
			} else if (*isym) == 1 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = (*j) - (*kuu); (*i) <= (*j)+(*klL); (*i)++ {
						(*temp) = (*dlatm3((m), (n), i, j, isub, jsub, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*(a))[(*isub)-(*jsub)+(*kuu)+0][(*jsub)-1] = (*temp)
						//Label270:
					}
					//Label280:
				}
			}
			//
		}
		//
	} else {
		//
		//        Use Dlatm2
		//
		if (*ipack) == 0 {
			if (*isym) == 0 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = 1; (*i) <= (*j); (*i)++ {
						(*(a))[(*i)-1][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						(*(a))[(*j)-1][(*i)-1] = (*(a))[(*i)-1][(*j)-1]
						//Label290:
					}
					//Label300:
				}
			} else if (*isym) == 1 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = 1; (*i) <= (*(m)); (*i)++ {
						(*(a))[(*i)-1][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						//Label310:
					}
					//Label320:
				}
			}
			//
		} else if (*ipack) == 1 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*(a))[(*i)-1][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					if (*i) != (*j) {
						(*(a))[(*j)-1][(*i)-1] = (*zero)
					}
					//Label330:
				}
				//Label340:
			}
			//
		} else if (*ipack) == 2 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*(a))[(*j)-1][(*i)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					if (*i) != (*j) {
						(*(a))[(*i)-1][(*j)-1] = (*zero)
					}
					//Label350:
				}
				//Label360:
			}
			//
		} else if (*ipack) == 3 {
			//
			(*isub) = 0
			(*jsub) = 1
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = 1; (*i) <= (*j); (*i)++ {
					(*isub) = (*isub) + 1
					if (*isub) > (*(lda)) {
						(*isub) = 1
						(*jsub) = (*jsub) + 1
					}
					(*(a))[(*isub)-1][(*jsub)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					//Label370:
				}
				//Label380:
			}
			//
		} else if (*ipack) == 4 {
			//
			if (*isym) == 0 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = 1; (*i) <= (*j); (*i)++ {
						//
						//                    Compute K = location of (I,j) entry in packed array
						//
						if (*i) == 1 {
							(*k) = (*j)
						} else {
							(*k) = (*(n))*((*(n))+1)/2 - ((*(n))-(*i)+1)*((*(n))-(*i)+2)/2 + (*j) - (*i) + 1
						}
						//
						//                    Convert K to (isub,jsub) location
						//
						(*jsub) = ((*k)-1)/(*(lda)) + 1
						(*isub) = (*k) - (*(lda))*((*jsub)-1)
						//
						(*(a))[(*isub)-1][(*jsub)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						//Label390:
					}
					//Label400:
				}
			} else {
				(*isub) = 0
				(*jsub) = 1
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = (*j); (*i) <= (*(m)); (*i)++ {
						(*isub) = (*isub) + 1
						if (*isub) > (*(lda)) {
							(*isub) = 1
							(*jsub) = (*jsub) + 1
						}
						(*(a))[(*isub)-1][(*jsub)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						//Label410:
					}
					//Label420:
				}
			}
			//
		} else if (*ipack) == 5 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
					if (*i) < 1 {
						(*(a))[(*j)-(*i)+0][(*i)+(*(n))-1] = (*zero)
					} else {
						(*(a))[(*j)-(*i)+0][(*i)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					}
					//Label430:
				}
				//Label440:
			}
			//
		} else if (*ipack) == 6 {
			//
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
					(*(a))[(*i)-(*j)+(*kuu)+0][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
					//Label450:
				}
				//Label460:
			}
			//
		} else if (*ipack) == 7 {
			//
			if (*isym) == 0 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = (*j) - (*kuu); (*i) <= (*j); (*i)++ {
						(*(a))[(*i)-(*j)+(*kuu)+0][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						if (*i) < 1 {
							(*(a))[(*j)-(*i)+1+(*kuu)-1][(*i)+(*(n))-1] = (*zero)
						}
						if (*i) >= 1 && (*i) != (*j) {
							(*(a))[(*j)-(*i)+1+(*kuu)-1][(*i)-1] = (*(a))[(*i)-(*j)+(*kuu)+0][(*j)-1]
						}
						//Label470:
					}
					//Label480:
				}
			} else if (*isym) == 1 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					for (*i) = (*j) - (*kuu); (*i) <= (*j)+(*klL); (*i)++ {
						(*(a))[(*i)-(*j)+(*kuu)+0][(*j)-1] = (*Dlatm2((m), (n), i, j, (kl), (ku), idist, (iseed), (d), igrade, (dl), (dr), ipvtng, (iwork), (sparse)))
						//Label490:
					}
					//Label500:
				}
			}
			//
		}
		//
	}
	//
	//     5)      Scaling the norm
	//
	if (*ipack) == 0 {
		(*onorm) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (m), (n), (a), (lda), tempa))
	} else if (*ipack) == 1 {
		(*onorm) = (*Dlansy(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('U'); return &y }(), (n), (a), (lda), tempa))
	} else if (*ipack) == 2 {
		(*onorm) = (*Dlansy(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('L'); return &y }(), (n), (a), (lda), tempa))
	} else if (*ipack) == 3 {
		(*onorm) = (*Dlansp(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('U'); return &y }(), (n), (a), tempa))
	} else if (*ipack) == 4 {
		(*onorm) = (*Dlansp(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('L'); return &y }(), (n), (a), tempa))
	} else if (*ipack) == 5 {
		(*onorm) = (*Dlansb(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('L'); return &y }(), (n), klL, (a), (lda), tempa))
	} else if (*ipack) == 6 {
		(*onorm) = (*Dlansb(func() *byte {y := byte('M'); return &y }(), func() *byte {y := byte('U'); return &y }(), (n), kuu, (a), (lda), tempa))
	} else if (*ipack) == 7 {
		(*onorm) = (*DLANGB(func() *byte {y := byte('M'); return &y }(), (n), klL, kuu, (a), (lda), tempa))
	}
	//
	if (*(anorm)) >= (*zero) {
		//
		if (*(anorm)) > (*zero) && (*onorm) == (*zero) {
			//
			//           Desired scaling impossible
			//
			(*(info)) = 5
			return
			//
		} else if ((*(anorm)) > (*one) && (*onorm) < (*one)) || ((*(anorm)) < (*one) && (*onorm) > (*one)) {
			//
			//           Scale carefully to avoid over / underflow
			//
			if (*ipack) <= 2 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					Dscal((m), (*one)/(*onorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					Dscal((m), (anorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					//Label510:
				}
				//
			} else if (*ipack) == 3 || (*ipack) == 4 {
				//
				Dscal((*(n))*((*(n))+1)/2, (*one)/(*onorm), (a), func() *int {y := 1; return &y }())
				Dscal((*(n))*((*(n))+1)/2, (anorm), (a), func() *int {y := 1; return &y }())
				//
			} else if (*ipack) >= 5 {
				//
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					Dscal((*klL)+(*kuu)+1, (*one)/(*onorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					Dscal((*klL)+(*kuu)+1, (anorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					//Label520:
				}
				//
			}
			//
		} else {
			//
			//           Scale straightforwardly
			//
			if (*ipack) <= 2 {
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					Dscal((m), (*(anorm))/(*onorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					//Label530:
				}
				//
			} else if (*ipack) == 3 || (*ipack) == 4 {
				//
				Dscal((*(n))*((*(n))+1)/2, (*(anorm))/(*onorm), (a), func() *int {y := 1; return &y }())
				//
			} else if (*ipack) >= 5 {
				//
				for (*j) = 1; (*j) <= (*(n)); (*j)++ {
					Dscal((*klL)+(*kuu)+1, (*(anorm))/(*onorm), &((*(a))[0][(*j)-1]), func() *int {y := 1; return &y }())
					//Label540:
				}
			}
			//
		}
		//
	}
	//
	//     End of DLAtmR
	//
}
