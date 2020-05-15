package goblas

import 

// Dchktz tests Dtzrzf.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchktz( dotype, nm, mval, nn, nval, thresh, tsterr, a,
//                          copya, s, tau, work, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nm, nn, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            mval(*), nval(*)
//       DOUBLE PRECISION   a(*), copya(*), S(*),
//      $                   tau(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchktz tests Dtzrzf.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] dotype
// \verbatim
//          dotype is LOGICAL array, dimension (ntypes)
//          The matrix types to be used for testing.  Matrices of type j
//          (for 1 <= j <= ntypes) are used for testing if dotype(j) =
//          .TRUE.; if dotype(j) = .FALSE., then type j is not used.
// \endverbatim
//
// \param[in] nm
// \verbatim
//          nm is intEGER
//          The number of values of M contained in the vector mval.
// \endverbatim
//
// \param[in] mval
// \verbatim
//          mval is intEGER array, dimension (nm)
//          The values of the matrix row dimension M.
// \endverbatim
//
// \param[in] nn
// \verbatim
//          nn is intEGER
//          The number of values of N contained in the vector nval.
// \endverbatim
//
// \param[in] nval
// \verbatim
//          nval is intEGER array, dimension (nn)
//          The values of the matrix column dimension N.
// \endverbatim
//
// \param[in] thresh
// \verbatim
//          thresh is DOUBLE PRECISION
//          The threshold value for the test ratios.  A result is
//          included in the output file if result >= thresh.  To have
//          every test ratio printed, use thresh = 0.
// \endverbatim
//
// \param[in] tsterr
// \verbatim
//          tsterr is LOGICAL
//          Flag that indicates whether error exits are to be tested.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (mmax*nmax)
//          where mmax is the maximum value of M in mval and nmax is the
//          maximum value of N in nval.
// \endverbatim
//
// \param[out] copya
// \verbatim
//          copya is DOUBLE PRECISION array, dimension (mmax*nmax)
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension
//                      (min(mmax,nmax))
// \endverbatim
//
// \param[out] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (mmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (mmax*nmax + 4*nmax + mmax)
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number for output.
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
func Dchktz(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, thresh *float64, tsterr *bool, a *[]float64, copya *[]float64, s *[]float64, tau *[]float64, work *[]float64, nout *int) {
	ntypes := new(int)
	ntests := new(int)
	one := new(float64)
	zero := new(float64)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	im := new(int)
	imode := new(int)
	in := new(int)
	info := new(int)
	k := new(int)
	lda := new(int)
	lwork := new(int)
	m := new(int)
	mnmin := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nrun := new(int)
	eps := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iseedy := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	result := func() *[]float64 {
		arr := make([]float64, 3)
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	infot := new(int)
	iounit := new(int)
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.iounit = new(float64)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new(int)
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
	(*ntypes) = 3
	(*ntests) = 3
	(*one) = 1.0
	(*zero) = 0.0
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
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
	infot = common.infoc.infot
	iounit = common.infoc.iounit
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Data statements ..
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3] = 1988, 1989, 1990, 1991
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y :=[]byte("TZ"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-1] = (*iseedy)[(*i)-1]
		//Label10:
	}
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     Test the error exits
	//
	if *(tsterr) {
		Derrtz(path, (nout))
	}
	(*infot) = 0
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		//
		//        Do for each value of M in mval.
		//
		(*m) = (*(mval))[(*im)-1]
		(*lda) = (MAX(1, (*m)))
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			//
			//           Do for each value of N in nval for which M .LE. N.
			//
			(*n) = (*(nval))[(*in)-1]
			(*mnmin) = (Min((*m), (*n)))
			(*lwork) = (*MAX(func() *int {y := 1; return &y }(), (*n)*(*n)+4*(*m)+(*n), (*m)*(*n)+2*(*mnmin)+4*(*n)))
			//
			if (*m) <= (*n) {
				for (*imode) = 1; (*imode) <= (*ntypes); (*imode)++ {
					if !(*(dotype))[(*imode)-1] {
						goto Label50
					}
					//
					//                 Do for each type of singular value distribution.
					//                    0:  zero matrix
					//                    1:  one small singular value
					//                    2:  exponential distribution
					//
					(*mode) = (*imode) - 1
					//
					//                 Test DTZRQF
					//
					//                 Generate test matrix of size m by n using
					//                 singular value distribution indicated by `mode'.
					//
					if (*mode) == 0 {
						Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, zero, zero, (a), lda)
						for (*i) = 1; (*i) <= (*mnmin); (*i)++ {
							(*(s))[(*i)-1] = (*zero)
							//Label30:
						}
					} else {
						Dlatms(m, n, func() *[]byte {y :=[]byte("Uniform"); return &y }(), iseed, func() *[]byte {y :=[]byte("Nonsymmetric"); return &y }(), (s), imode, (*one)/(*eps), one, m, n, func() *[]byte {y :=[]byte("No packing"); return &y }(), (a), lda, (work), info)
						Dgeqr2(m, n, (a), lda, (work), &((*(work))[(*mnmin)+0]), info)
						Dlaset(func() *[]byte {y :=[]byte("Lower"); return &y }(), (*m)-1, n, zero, zero, &((*(a))[1]), lda)
						Dlaord(func() *[]byte {y :=[]byte("Decreasing"); return &y }(), mnmin, (s), func() *int {y := 1; return &y }())
					}
					//
					//                 Save A and its singular values
					//
					Dlacpy(func() *[]byte {y :=[]byte("All"); return &y }(), m, n, (a), lda, (copya), lda)
					//
					//                 Call Dtzrzf to reduce the upper trapezoidal matrix to
					//                 upper triangular form.
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dtzrzf"); return &y }()
					Dtzrzf(m, n, (a), lda, (tau), (work), lwork, info)
					//
					//                 Compute norm(svd(a) - svd(r))
					//
					(*result)[0] = (*Dqrt12(m, m, (a), lda, (s), (work), lwork))
					//
					//                 Compute norm( A - R*Q)
					//
					(*result)[1] = (*Drzt01(m, n, (copya), (a), lda, (tau), (work), lwork))
					//
					//                 Compute norm(Q'*Q - I).
					//
					(*result)[2] = (*Drzt02(m, n, (a), lda, (tau), (work), lwork))
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = 1; (*k) <= (*ntests); (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {y :=[]byte(" M =%5d, N =%5d, type %2d, test %2d, ratio =%12.5f\n"); return &y }(), (*m), (*n), (*imode), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label40:
					}
					(*nrun) = (*nrun) + 3
				Label50:
				}
			}
			//Label60:
		}
		//Label70:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	//
	//     End if Dchktz
	//
}
