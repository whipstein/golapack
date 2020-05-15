package goblas

import 

// Dchkrq tests Dgerqf, Dorgrq and DORMRQ.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkrq( dotype, nm, mval, nn, nval, nnb, nbval, nxval,
//                          nrhs, thresh, tsterr, nmax, a, af, aq, ar, ac,
//                          b, x, xact, tau, work, rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nm, nmax, nn, nnb, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), mval(*), nbval(*), nval(*),
//      $                   nxval(*)
//       DOUBLE PRECISION   a(*), Ac(*), af(*), aq(*), Ar(*),
//      $                   B(*), rwork(*), tau(*), work(*),
//      $                   X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkrq tests Dgerqf, Dorgrq and DORMRQ.
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
// \param[in] nnb
// \verbatim
//          nnb is intEGER
//          The number of values of nb and nx contained in the
//          vectors nbval and nxval.  The blocking parameters are used
//          in pairs (nb,nx).
// \endverbatim
//
// \param[in] nbval
// \verbatim
//          nbval is intEGER array, dimension (nnb)
//          The values of the blocksize nb.
// \endverbatim
//
// \param[in] nxval
// \verbatim
//          nxval is intEGER array, dimension (nnb)
//          The values of the crossover point nx.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand side vectors to be generated for
//          each linear system.
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
// \param[in] nmax
// \verbatim
//          nmax is intEGER
//          The maximum value permitted for M or n, used in dimensioning
//          the work arrays.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] aq
// \verbatim
//          aq is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] AR
// \verbatim
//          AR is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] AC
// \verbatim
//          AC is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (nmax*nrhs)
// \endverbatim
//
// \param[out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (nmax*nrhs)
// \endverbatim
//
// \param[out] xact
// \verbatim
//          xact is DOUBLE PRECISION array, dimension (nmax*nrhs)
// \endverbatim
//
// \param[out] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (nmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (nmax)
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is intEGER array, dimension (nmax)
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
func Dchkrq(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nxval *[]int, nrhs *int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, af *[]float64, aq *[]float64, ar *[]float64, ac *[]float64, b *[]float64, x *[]float64, xact *[]float64, tau *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	ntests := new(int)
	ntypes := new(int)
	zero := new(float64)
	dist := new(byte)
	_type := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	ik := new(int)
	im := new(int)
	imat := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	k := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	lwork := new(int)
	m := new(int)
	minmn := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nk := new(int)
	nrun := new(int)
	nt := new(int)
	nx := new(int)
	anorm := new(float64)
	cndnum := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iseedy := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	kval := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	result := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	infot := new(int)
	nunit := new(int)
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.nunit = new(float64)
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
	(*ntests) = 7
	(*ntypes) = 8
	(*zero) = 0.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
	infot = common.infoc.infot
	nunit = common.infoc.nunit
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
	(*path)[1] = *func() *[]byte {y :=[]byte("RQ"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-1] = (*iseedy)[(*i)-1]
		//Label10:
	}
	//
	//     Test the error exits
	//
	if *(tsterr) {
		Derrrq(path, (nout))
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	//
	(*lda) = (*(nmax))
	(*lwork) = (*(nmax)) * MAX((*(nmax)), (*(nrhs)))
	//
	//     Do for each value of M in mval.
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		(*m) = (*(mval))[(*im)-1]
		//
		//        Do for each value of N in nval.
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			(*n) = (*(nval))[(*in)-1]
			(*minmn) = (Min((*m), (*n)))
			for (*imat) = 1; (*imat) <= (*ntypes); (*imat)++ {
				//
				//              Do the tests only if dotype( imat) is true.
				//
				if !(*(dotype))[(*imat)-1] {
					goto Label50
				}
				//
				//              Set up parameters with Dlatb4 and generate a test matrix
				//              with Dlatms.
				//
				Dlatb4(path, imat, m, n, _type, kl, ku, anorm, mode, cndnum, dist)
				//
				(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y }()
				Dlatms(m, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, func() *[]byte {y :=[]byte("No packing"); return &y }(), (a), lda, (work), info)
				//
				//              Check error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					goto Label50
				}
				//
				//              Set some values for K: the first value must be minmn,
				//              corresponding to the call of drqt01; other values are
				//              used in the calls of drqt02, and must not exceed minmn.
				//
				(*kval)[0] = (*minmn)
				(*kval)[1] = 0
				(*kval)[2] = 1
				(*kval)[3] = (*minmn) / 2
				if (*minmn) == 0 {
					(*nk) = 1
				} else if (*minmn) == 1 {
					(*nk) = 2
				} else if (*minmn) <= 3 {
					(*nk) = 3
				} else {
					(*nk) = 4
				}
				//
				//              Do for each value of K in kval
				//
				for (*ik) = 1; (*ik) <= (*nk); (*ik)++ {
					(*k) = (*kval)[(*ik)-1]
					//
					//                 Do for each pair of values (nb,nx) in nbval and nxval.
					//
					for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
						(*nb) = (*(nbval))[(*inb)-1]
						Xlaenv(func() *int {y := 1; return &y }(), nb)
						(*nx) = (*(nxval))[(*inb)-1]
						Xlaenv(func() *int {y := 3; return &y }(), nx)
						for (*i) = 1; (*i) <= (*ntests); (*i)++ {
							(*result)[(*i)-1] = (*zero)
						}
						(*nt) = 2
						if (*ik) == 1 {
							//
							//                       Test Dgerqf
							//
							drqt01(m, n, (a), (af), (aq), (ar), lda, (tau), (work), lwork, (rwork), &((*result)[0]))
						} else if (*m) <= (*n) {
							//
							//                       Test Dorgrq, using factorization
							//                       returned by drqt01
							//
							drqt02(m, n, k, (a), (af), (aq), (ar), lda, (tau), (work), lwork, (rwork), &((*result)[0]))
						}
						if (*m) >= (*k) {
							//
							//                       Test Dormrq, using factorization returned
							//                       by drqt01
							//
							drqt03(m, n, k, (af), (ac), (ar), (aq), lda, (tau), (work), lwork, (rwork), &((*result)[2]))
							(*nt) = (*nt) + 4
							//
							//                       If M>=N and K=N, call Dgerqs to solve a system
							//                       with nrhs right hand sides and compute the
							//                       residual.
							//
							if (*k) == (*m) && (*inb) == 1 {
								//
								//                          Generate a solution and set the right
								//                          hand side.
								//
								(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y }()
								Dlarhs(path, func() *[]byte {y :=[]byte("New"); return &y }(), func() *[]byte {y :=[]byte("Full"); return &y }(), func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), (nrhs), (a), lda, (xact), lda, (b), lda, iseed, info)
								//
								Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, (nrhs), (b), lda, &((*(x))[(*n)-(*m)+0]), lda)
								(*srnamt) = *func() *[]byte {y :=[]byte("Dgerqs"); return &y }()
								Dgerqs(m, n, (nrhs), (af), lda, (tau), (x), lda, (work), lwork, info)
								//
								//                          Check error code from Dgerqs.
								//
								if (*info) != 0 {
									Alaerh(path, func() *[]byte {y :=[]byte("Dgerqs"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, (nrhs), -1, nb, imat, nfail, nerrs, (nout))
								}
								//
								Dget02(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, (nrhs), (a), lda, (x), lda, (b), lda, (rwork), &((*result)[6]))
								(*nt) = (*nt) + 1
							}
						}
						//
						//                    Print information about the tests that did not
						//                    pass the threshold.
						//
						for (*i) = 1; (*i) <= (*nt); (*i)++ {
							if (*result)[(*i)-1] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Alahd((nout), path)
								}
								WRITE((*(nout)), *func() *[]byte {
									y :=[]byte(" M=%5d, N=%5d, K=%5d, nb=%4d, nx=%5d, type %2d, test(%2d)=%12.5f\n")
									return &y
								}(), (*m), (*n), (*k), (*nb), (*nx), (*imat), (*i), (*result)[(*i)-1])
								(*nfail) = (*nfail) + 1
							}
							//Label20:
						}
						(*nrun) = (*nrun) + (*nt)
						//Label30:
					}
					//Label40:
				}
			Label50:
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
	return
	//
	//     End of Dchkrq
	//
}
