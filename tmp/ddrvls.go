package goblas

import 

// Ddrvls tests the least squares driver routines Dgels, Dgetsls, Dgelss, Dgelsy,
// and DgelsD.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvls( dotype, NM, mval, nn, nval, nns, nsval, nnb,
//                          nbval, nxVAL, thresh, tsterr, a, copya, b,
//                          copyb, c, S, copys, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       inTEGER            NM, nn, nnb, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       inTEGER            mval(*), nbval(*), nsval(*),
//      $                   nval(*), nxVAL(*)
//       DOUBLE PRECISION   a(*), B(*), c(*), copya(*), copyb(*),
//      $                   copys(*), S(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvls tests the least squares driver routines Dgels, Dgetsls, Dgelss, Dgelsy,
// and DgelsD.
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
//          The matrix of type j is generated as follows:
//          j=1: A = U*D*V where U and V are random orthogonal matrices
//               and D has random entries (> 0.1) taken from a uniform
//               distribution (0,1). A is full rank.
//          j=2: The same of 1, but A is scaled up.
//          j=3: The same of 1, but A is scaled down.
//          j=4: A = U*D*V where U and V are random orthogonal matrices
//               and D has 3*min(m,N)/4 random entries (> 0.1) taken
//               from a uniform distribution (0,1) and the remaining
//               entries set to 0. A is rank-deficient.
//          j=5: The same of 4, but A is scaled up.
//          j=6: The same of 5, but A is scaled down.
// \endverbatim
//
// \param[in] NM
// \verbatim
//          NM is inTEGER
//          The number of values of M contained in the vector mval.
// \endverbatim
//
// \param[in] mval
// \verbatim
//          mval is inTEGER array, dimension (nm)
//          The values of the matrix row dimension M.
// \endverbatim
//
// \param[in] nn
// \verbatim
//          nn is inTEGER
//          The number of values of N contained in the vector nval.
// \endverbatim
//
// \param[in] nval
// \verbatim
//          nval is inTEGER array, dimension (nn)
//          The values of the matrix column dimension N.
// \endverbatim
//
// \param[in] nns
// \verbatim
//          nns is inTEGER
//          The number of values of nrhs contained in the vector nsval.
// \endverbatim
//
// \param[in] nsval
// \verbatim
//          nsval is inTEGER array, dimension (nns)
//          The values of the number of right hand sides nrhs.
// \endverbatim
//
// \param[in] nnb
// \verbatim
//          nnb is inTEGER
//          The number of values of nb and nx contained in the
//          vectors nbval and nxVAL.  The blocking parameters are used
//          in pairs (nb,nx).
// \endverbatim
//
// \param[in] nbval
// \verbatim
//          nbval is inTEGER array, dimension (nnb)
//          The values of the blocksize nb.
// \endverbatim
//
// \param[in] nxVAL
// \verbatim
//          nxVAL is inTEGER array, dimension (nnb)
//          The values of the crossover point nx.
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
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (mmax*nsmax)
//          where mmax is the maximum value of M in mval and NSMAX is the
//          maximum value of nrhs in nsval.
// \endverbatim
//
// \param[out] copyb
// \verbatim
//          copyb is DOUBLE PRECISION array, dimension (mmax*nsmax)
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION array, dimension (mmax*nsmax)
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension
//                      (min(mmax,nmax))
// \endverbatim
//
// \param[out] copys
// \verbatim
//          copys is DOUBLE PRECISION array, dimension
//                      (min(mmax,nmax))
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is inTEGER
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
// \date June 2017
//
// \ingroup double_lin
//
//  =====================================================================
func Ddrvls(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, nns *int, nsval *[]int, nnb *int, nbval *[]int, nxval *[]int, thresh *float64, tsterr *bool, a *[]float64, copya *[]float64, b *[]float64, copyb *[]float64, c *[]float64, s *[]float64, copys *[]float64, nout *int) {
	ntests := new(int)
	smlsiz := new(int)
	one := new(float64)
	two := new(float64)
	zero := new(float64)
	trans := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	crank := new(int)
	i := new(int)
	im := new(int)
	imb := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	ins := new(int)
	irank := new(int)
	iscale := new(int)
	itran := new(int)
	itype := new(int)
	j := new(int)
	k := new(int)
	lda := new(int)
	ldb := new(int)
	ldwork := new(int)
	lwlsy := new(int)
	lwork := new(int)
	m := new(int)
	mnmin := new(int)
	n := new(int)
	nb := new(int)
	ncols := new(int)
	nerrs := new(int)
	nfail := new(int)
	nrhs := new(int)
	nrows := new(int)
	nrun := new(int)
	rank := new(int)
	mb := new(int)
	mmax := new(int)
	nmax := new(int)
	nsmax := new(int)
	liwork := new(int)
	lworkDgels := new(int)
	lworkDgetsls := new(int)
	lworkDgelss := new(int)
	lworkDgelsy := new(int)
	lworkDgelsd := new(int)
	eps := new(float64)
	norma := new(float64)
	normb := new(float64)
	rcond := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iseedy := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iwq := func() *[]int {
		arr := make([]int, 1)
		return &arr
	}()
	result := func() *[]float64 {
		arr := make([]float64, 16)
		return &arr
	}()
	wq := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	allocatable := new(float64)
	allocatable := new(int)
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
	//  -- lapACK test routine (version 3.7.1) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2017
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*ntests) = 16
	(*smlsiz) = 25
	(*one) = 1.0
	(*two) = 2.0
	(*zero) = 0.0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Allocatable Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Scalars in Common ..
	//     ..
	//     .. Common blocks ..
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
	(*path)[1] = *func() *[]byte {y :=[]byte("LS"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-1] = (*iseedy)[(*i)-1]
		//Label10:
	}
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	//
	//     Threshold for rank estimation
	//
	(*rcond) = SQRt((*eps)) - (SQRt((*eps))-(*eps))/2
	//
	//     Test the error exits
	//
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	Xlaenv(func() *int {y := 9; return &y }(), smlsiz)
	if *(tsterr) {
		Derrls(path, (nout))
	}
	//
	//     Print the header if NM = 0 or nn = 0 and thresh = 0.
	//
	if ((*(nm)) == 0 || (*(nn)) == 0) && (*(thresh)) == (*zero) {
		Alahd((nout), path)
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	Xlaenv(func() *int {y := 9; return &y }(), smlsiz)
	//
	//     Compute maximal workspace needed for all routines
	//
	(*nmax) = 0
	(*mmax) = 0
	(*nsmax) = 0
	for (*i) = 1; (*i) <= (*(nm)); (*i)++ {
		if (*(mval))[(*i)-1] > (*mmax) {
			(*mmax) = (*(mval))[(*i)-1]
		}
	}
	for (*i) = 1; (*i) <= (*(nn)); (*i)++ {
		if (*(nval))[(*i)-1] > (*nmax) {
			(*nmax) = (*(nval))[(*i)-1]
		}
	}
	for (*i) = 1; (*i) <= (*(nns)); (*i)++ {
		if (*(nsval))[(*i)-1] > (*nsmax) {
			(*nsmax) = (*(nsval))[(*i)-1]
		}
	}
	(*m) = (*mmax)
	(*n) = (*nmax)
	(*nrhs) = (*nsmax)
	(*mnmin) = (MAX(Min((*m), (*n)), 1))
	//
	//     Compute workspace needed for routines
	//     Dqrt14, Dqrt17 (two side cases), Dqrt15 and Dqrt12
	//
	(*lwork) = (*MAX(func() *int {y := 1; return &y }(), ((*m)+(*n))*(*nrhs), ((*n)+(*nrhs))*((*m)+2), ((*m)+(*nrhs))*((*n)+2), MAX((*m)+(*mnmin), (*nrhs)*(*mnmin), 2*(*n)+(*m)), MAX((*m)*(*n)+4*(*mnmin)+MAX((*m), (*n)), (*m)*(*n)+2*(*mnmin)+4*(*n))))
	(*liwork) = 1
	//
	//     Iterate through all test cases and compute necessary workspace
	//     sizes for ?GELS, ?GETSLS, ?GELSY, ?GELSS and ?GELSD routines.
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		(*m) = (*(mval))[(*im)-1]
		(*lda) = (MAX(1, (*m)))
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			(*n) = (*(nval))[(*in)-1]
			(*mnmin) = (MAX(Min((*m), (*n)), 1))
			(*ldb) = (*MAX(func() *int {y := 1; return &y }(), m, N))
			for (*inS) = 1; (*inS) <= (*(nns)); (*inS)++ {
				(*nrhs) = (*(nsval))[(*inS)-1]
				for (*irank) = 1; (*irank) <= 2; (*irank)++ {
					for (*iscale) = 1; (*iscale) <= 3; (*iscale)++ {
						(*itype) = ((*irank)-1)*3 + (*iscale)
						if (*(dotype))[(*itype)-1] {
							if (*irank) == 1 {
								for (*itran) = 1; (*itran) <= 2; (*itran)++ {
									if (*itran) == 1 {
										(*trans) = 'N'
									} else {
										(*trans) = 'T'
									}
									//
									//                             Compute workspace needed for Dgels
									Dgels(trans, m, n, nrhs, (a), lda, (b), ldb, wq, -1, info)
									(*lwork_Dgels) = (*int(&((*wq)[0])))
									//                             Compute workspace needed for Dgetsls
									Dgetsls(trans, m, n, nrhs, (a), lda, (b), ldb, wq, -1, info)
									(*lwork_Dgetsls) = (*int(&((*wq)[0])))
								}
							}
							//                       Compute workspace needed for Dgelsy
							Dgelsy(m, n, nrhs, (a), lda, (b), ldb, iwq, rcond, crank, wq, -1, info)
							(*lworkDgelsy) = (*int(&((*wq)[0])))
							//                       Compute workspace needed for Dgelss
							Dgelss(m, n, nrhs, (a), lda, (b), ldb, (s), rcond, crank, wq, -1, info)
							(*lworkDgelss) = (*int(&((*wq)[0])))
							//                       Compute workspace needed for DgelsD
							Dgelsd(m, n, nrhs, (a), lda, (b), ldb, (s), rcond, crank, wq, -1, iwq, info)
							(*lworkDgelsd) = (*int(&((*wq)[0])))
							//                       Compute liwork workspace needed for Dgelsy and DgelsD
							(*liwork) = (*MAX(liwork, n, &((*iwq)[0])))
							//                       Compute lwork workspace needed for all functions
							(*lwork) = (*MAX(lwork, lwork_Dgels, lwork_Dgetsls, lworkDgelsy, lworkDgelss, lworkDgelsd))
						}
					}
				}
			}
		}
	}
	//
	(*lwlsy) = (*lwork)
	//
	ALLOCATE(work(lwork))
	ALLOCATE(iwork(liwork))
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		(*m) = (*(mval))[(*im)-1]
		(*lda) = (MAX(1, (*m)))
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			(*n) = (*(nval))[(*in)-1]
			(*mnmin) = (MAX(Min((*m), (*n)), 1))
			(*ldb) = (*MAX(func() *int {y := 1; return &y }(), m, N))
			(*mb) = ((*mnmin) + 1)
			//
			for (*inS) = 1; (*inS) <= (*(nns)); (*inS)++ {
				(*nrhs) = (*(nsval))[(*inS)-1]
				//
				for (*irank) = 1; (*irank) <= 2; (*irank)++ {
					for (*iscale) = 1; (*iscale) <= 3; (*iscale)++ {
						(*itype) = ((*irank)-1)*3 + (*iscale)
						if !(*(dotype))[(*itype)-1] {
							goto Label110
						}
						//
						if (*irank) == 1 {
							//
							//                       Test Dgels
							//
							//                       Generate a matrix of scaling type iscale
							//
							Dqrt13(iscale, m, n, (copya), lda, norma, iseed)
							for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
								(*nb) = (*(nbval))[(*inb)-1]
								Xlaenv(func() *int {y := 1; return &y }(), nb)
								Xlaenv(func() *int {y := 3; return &y }(), &((*(nxVAL))[(*inb)-1]))
								//
								for (*itran) = 1; (*itran) <= 2; (*itran)++ {
									if (*itran) == 1 {
										(*trans) = 'N'
										(*nrows) = (*m)
										(*ncols) = (*n)
									} else {
										(*trans) = 'T'
										(*nrows) = (*n)
										(*ncols) = (*m)
									}
									(*ldwork) = (MAX(1, (*ncols)))
									//
									//                             Set up a consistent rhs
									//
									if (*ncols) > 0 {
										Dlarnv(func() *int {y := 2; return &y }(), iseed, (*ncols)*(*nrhs), work)
										Dscal((*ncols)*(*nrhs), (*one)/DBLE((*ncols)), work, func() *int {y := 1; return &y }())
									}
									Dgemm(trans, func() *[]byte {y :=[]byte("No transpose"); return &y }(), nrows, nrhs, ncols, one, (copya), lda, work, ldwork, zero, (b), ldb)
									Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (b), ldb, (copyb), ldb)
									//
									//                             Solve LS or overdetermined system
									//
									if (*m) > 0 && (*n) > 0 {
										Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (copya), lda, (a), lda)
										Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (copyb), ldb, (b), ldb)
									}
									(*srnamt) = *func() *[]byte {y :=[]byte("Dgels "); return &y }()
									Dgels(trans, m, n, nrhs, (a), lda, (b), ldb, work, lwork, info)
									if (*info) != 0 {
										Alaerh(path, func() *[]byte {y :=[]byte("Dgels "); return &y }(), info, func() *int {y := 0; return &y }(), trans, m, n, nrhs, -1, nb, itype, nfail, nerrs, (nout))
									}
									//
									//                             Check correctness of results
									//
									(*ldwork) = (MAX(1, (*nrows)))
									if (*nrows) > 0 && (*nrhs) > 0 {
										Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (copyb), ldb, (c), ldb)
									}
									Dqrt16(trans, m, n, nrhs, (copya), lda, (b), ldb, (c), ldb, work, &((*result)[0]))
									//
									if ((*itran) == 1 && (*m) >= (*n)) || ((*itran) == 2 && (*m) < (*n)) {
										//
										//                                Solving LS system
										//
										(*result)[1] = (*Dqrt17(trans, func() *int {y := 1; return &y }(), m, n, nrhs, (copya), lda, (b), ldb, (copyb), ldb, (c), work, lwork))
									} else {
										//
										//                                Solving overdetermined system
										//
										(*result)[1] = (*Dqrt14(trans, m, n, nrhs, (copya), lda, (b), ldb, work, lwork))
									}
									//
									//                             Print information about the tests that
									//                             did not pass the threshold.
									//
									for (*k) = 1; (*k) <= 2; (*k)++ {
										if (*result)[(*k)-1] >= (*(thresh)) {
											if (*nfail) == 0 && (*nerrs) == 0 {
												Alahd((nout), path)
											}
											WRITE((*(nout)), *func() *[]byte {
												y :=[]byte(" trans='%c', M=%5d, N=%5d, nrhs=%4d, nb=%4d, type%2d, test(%2d)=%12.5f\n")
												return &y
											}(), (*trans), (*m), (*n), (*nrhs), (*nb), (*itype), (*k), (*result)[(*k)-1])
											(*nfail) = (*nfail) + 1
										}
										//Label20:
									}
									(*nrun) = (*nrun) + 2
									//Label30:
								}
								//Label40:
							}
							//
							//
							//                       Test Dgetsls
							//
							//                       Generate a matrix of scaling type iscale
							//
							Dqrt13(iscale, m, n, (copya), lda, norma, iseed)
							for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
								(*mb) = (*(nbval))[(*inb)-1]
								Xlaenv(func() *int {y := 1; return &y }(), mb)
								for (*imb) = 1; (*imb) <= (*(nnb)); (*imb)++ {
									(*nb) = (*(nbval))[(*imb)-1]
									Xlaenv(func() *int {y := 2; return &y }(), nb)
									//
									for (*itran) = 1; (*itran) <= 2; (*itran)++ {
										if (*itran) == 1 {
											(*trans) = 'N'
											(*nrows) = (*m)
											(*ncols) = (*n)
										} else {
											(*trans) = 'T'
											(*nrows) = (*n)
											(*ncols) = (*m)
										}
										(*ldwork) = (MAX(1, (*ncols)))
										//
										//                             Set up a consistent rhs
										//
										if (*ncols) > 0 {
											Dlarnv(func() *int {y := 2; return &y }(), iseed, (*ncols)*(*nrhs), work)
											Dscal((*ncols)*(*nrhs), (*one)/DBLE((*ncols)), work, func() *int {y := 1; return &y }())
										}
										Dgemm(trans, func() *[]byte {y :=[]byte("No transpose"); return &y }(), nrows, nrhs, ncols, one, (copya), lda, work, ldwork, zero, (b), ldb)
										Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (b), ldb, (copyb), ldb)
										//
										//                             Solve LS or overdetermined system
										//
										if (*m) > 0 && (*n) > 0 {
											Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (copya), lda, (a), lda)
											Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (copyb), ldb, (b), ldb)
										}
										(*srnamt) = *func() *[]byte {y :=[]byte("Dgetsls "); return &y }()
										Dgetsls(trans, m, n, nrhs, (a), lda, (b), ldb, work, lwork, info)
										if (*info) != 0 {
											Alaerh(path, func() *[]byte {y :=[]byte("Dgetsls "); return &y }(), info, func() *int {y := 0; return &y }(), trans, m, n, nrhs, -1, nb, itype, nfail, nerrs, (nout))
										}
										//
										//                             Check correctness of results
										//
										(*ldwork) = (MAX(1, (*nrows)))
										if (*nrows) > 0 && (*nrhs) > 0 {
											Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), nrows, nrhs, (copyb), ldb, (c), ldb)
										}
										Dqrt16(trans, m, n, nrhs, (copya), lda, (b), ldb, (c), ldb, work, &((*result)[14]))
										//
										if ((*itran) == 1 && (*m) >= (*n)) || ((*itran) == 2 && (*m) < (*n)) {
											//
											//                                Solving LS system
											//
											(*result)[15] = (*Dqrt17(trans, func() *int {y := 1; return &y }(), m, n, nrhs, (copya), lda, (b), ldb, (copyb), ldb, (c), work, lwork))
										} else {
											//
											//                                Solving overdetermined system
											//
											(*result)[15] = (*Dqrt14(trans, m, n, nrhs, (copya), lda, (b), ldb, work, lwork))
										}
										//
										//                             Print information about the tests that
										//                             did not pass the threshold.
										//
										for (*k) = 15; (*k) <= 16; (*k)++ {
											if (*result)[(*k)-1] >= (*(thresh)) {
												if (*nfail) == 0 && (*nerrs) == 0 {
													Alahd((nout), path)
												}
												WRITE((*(nout)), *func() *[]byte {
													y :=[]byte(" trans='%c M=%5d, N=%5d, nrhs=%4d, mb=%4d, nb=%4d, type%2d, test(%2d)=%12.5f\n")
													return &y
												}(), (*trans), (*m), (*n), (*nrhs), (*mb), (*nb), (*itype), (*k), (*result)[(*k)-1])
												(*nfail) = (*nfail) + 1
											}
											//Label50:
										}
										(*nrun) = (*nrun) + 2
										//Label60:
									}
									//Label62:
								}
								//Label65:
							}
						}
						//
						//                    Generate a matrix of scaling type iscale and rank
						//                    type irank.
						//
						Dqrt15(iscale, irank, m, n, nrhs, (copya), lda, (copyb), ldb, (copys), rank, norma, normb, iseed, work, lwork)
						//
						//                    workspace used: MAX(M+Min(m,N),nrhs*Min(m,N),2*n+M)
						//
						(*ldwork) = (MAX(1, (*m)))
						//
						//                    Loop for testing different block sizes.
						//
						for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
							(*nb) = (*(nbval))[(*inb)-1]
							Xlaenv(func() *int {y := 1; return &y }(), nb)
							Xlaenv(func() *int {y := 3; return &y }(), &((*(nxVAL))[(*inb)-1]))
							//
							//                       Test Dgelsy
							//
							//                       Dgelsy:  Compute the minimum-norm solution X
							//                       to min( norm( A * X - B))
							//                       using the rank-revealing orthogonal
							//                       factorization.
							//
							//                       Initialize vector iwork.
							//
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								iwork(j) = 0
								//Label70:
							}
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (copya), lda, (a), lda)
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, (b), ldb)
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("Dgelsy"); return &y }()
							Dgelsy(m, n, nrhs, (a), lda, (b), ldb, iwork, rcond, crank, work, lwlsy, info)
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("Dgelsy"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, nrhs, -1, nb, itype, nfail, nerrs, (nout))
							}
							//
							//                       Test 3:  Compute relative error in svd
							//                                workspace: M*n + 4*Min(m,N) + MAX(m,N)
							//
							(*result)[2] = (*Dqrt12(crank, crank, (a), lda, (copys), work, lwork))
							//
							//                       Test 4:  Compute error in solution
							//                                workspace:  M*nrhs + M
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, work, ldwork)
							Dqrt16(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, ldwork, work((*m)*(*nrhs)+1), &((*result)[3]))
							//
							//                       Test 5:  Check norm of r'*A
							//                                workspace: nrhs*(M+N)
							//
							(*result)[4] = (*zero)
							if (*m) > (*crank) {
								(*result)[4] = (*Dqrt17(func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *int {y := 1; return &y }(), m, n, nrhs, (copya), lda, (b), ldb, (copyb), ldb, (c), work, lwork))
							}
							//
							//                       Test 6:  Check if x is in the rowspace of A
							//                                workspace: (M+nrhs)*(N+2)
							//
							(*result)[5] = (*zero)
							//
							if (*n) > (*crank) {
								(*result)[5] = (*Dqrt14(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, lwork))
							}
							//
							//                       Test Dgelss
							//
							//                       Dgelss:  Compute the minimum-norm solution X
							//                       to min( norm( A * X - B))
							//                       using the SVD.
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (copya), lda, (a), lda)
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, (b), ldb)
							(*srnamt) = *func() *[]byte {y :=[]byte("Dgelss"); return &y }()
							Dgelss(m, n, nrhs, (a), lda, (b), ldb, (s), rcond, crank, work, lwork, info)
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("Dgelss"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, nrhs, -1, nb, itype, nfail, nerrs, (nout))
							}
							//
							//                       workspace used: 3*min(m,n) +
							//                                       max(2*min(m,n),nrhs,max(m,n))
							//
							//                       Test 7:  Compute relative error in svd
							//
							if (*rank) > 0 {
								Daxpy(mnmin, -(*one), (copys), func() *int {y := 1; return &y }(), (s), func() *int {y := 1; return &y }())
								(*result)[6] = Dasum(mnmin, (s), func() *int {y := 1; return &y }()) / Dasum(mnmin, (copys), func() *int {y := 1; return &y }()) / ((*eps) * DBLE((*mnmin)))
							} else {
								(*result)[6] = (*zero)
							}
							//
							//                       Test 8:  Compute error in solution
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, work, ldwork)
							Dqrt16(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, ldwork, work((*m)*(*nrhs)+1), &((*result)[7]))
							//
							//                       Test 9:  Check norm of r'*A
							//
							(*result)[8] = (*zero)
							if (*m) > (*crank) {
								(*result)[8] = (*Dqrt17(func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *int {y := 1; return &y }(), m, n, nrhs, (copya), lda, (b), ldb, (copyb), ldb, (c), work, lwork))
							}
							//
							//                       Test 10:  Check if x is in the rowspace of A
							//
							(*result)[9] = (*zero)
							if (*n) > (*crank) {
								(*result)[9] = (*Dqrt14(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, lwork))
							}
							//
							//                       Test DgelsD
							//
							//                       DgelsD:  Compute the minimum-norm solution X
							//                       to min( norm( A * X - B)) using a
							//                       divide and conquer SVD.
							//
							//                       Initialize vector iwork.
							//
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								iwork(j) = 0
								//Label80:
							}
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (copya), lda, (a), lda)
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, (b), ldb)
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("DgelsD"); return &y }()
							Dgelsd(m, n, nrhs, (a), lda, (b), ldb, (s), rcond, crank, work, lwork, iwork, info)
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("DgelsD"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, nrhs, -1, nb, itype, nfail, nerrs, (nout))
							}
							//
							//                       Test 11:  Compute relative error in svd
							//
							if (*rank) > 0 {
								Daxpy(mnmin, -(*one), (copys), func() *int {y := 1; return &y }(), (s), func() *int {y := 1; return &y }())
								(*result)[10] = Dasum(mnmin, (s), func() *int {y := 1; return &y }()) / Dasum(mnmin, (copys), func() *int {y := 1; return &y }()) / ((*eps) * DBLE((*mnmin)))
							} else {
								(*result)[10] = (*zero)
							}
							//
							//                       Test 12:  Compute error in solution
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, nrhs, (copyb), ldb, work, ldwork)
							Dqrt16(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, ldwork, work((*m)*(*nrhs)+1), &((*result)[11]))
							//
							//                       Test 13:  Check norm of r'*A
							//
							(*result)[12] = (*zero)
							if (*m) > (*crank) {
								(*result)[12] = (*Dqrt17(func() *[]byte {y :=[]byte("No transpose"); return &y }(), func() *int {y := 1; return &y }(), m, n, nrhs, (copya), lda, (b), ldb, (copyb), ldb, (c), work, lwork))
							}
							//
							//                       Test 14:  Check if x is in the rowspace of A
							//
							(*result)[13] = (*zero)
							if (*n) > (*crank) {
								(*result)[13] = (*Dqrt14(func() *[]byte {y :=[]byte("No transpose"); return &y }(), m, n, nrhs, (copya), lda, (b), ldb, work, lwork))
							}
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = 3; (*k) <= 14; (*k)++ {
								if (*result)[(*k)-1] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Alahd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte {y :=[]byte(" M=%5d, N=%5d, nrhs=%4d, nb=%4d, type%2d, test(%2d)=%12.5f\n"); return &y }(), (*m), (*n), (*nrhs), (*nb), (*itype), (*k), (*result)[(*k)-1])
									(*nfail) = (*nfail) + 1
								}
								//Label90:
							}
							(*nrun) = (*nrun) + 12
							//
							//Label100:
						}
					Label110:
					}
					//Label120:
				}
				//Label130:
			}
			//Label140:
		}
		//Label150:
	}
	//
	//     Print a summary of the results.
	//
	Alasvm(path, (nout), nfail, nrun, nerrs)
	//
	//
	DEALLOCATE(work)
	DEALLOCATE(iwork)
	return
	//
	//     End of Ddrvls
	//
}
