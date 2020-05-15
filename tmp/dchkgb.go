package goblas

import 

// Dchkgb tests Dgbtrf, -TRS, -RFS, and -CON
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkgb( dotype, nm, mval, nn, nval, nnb, nbval, nns,
//                          nsval, thresh, tsterr, a, la, afac, lafac, b,
//                          x, xact, work, rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            la, lafac, nm, nn, nnb, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), mval(*), nbval(*), nsval(*),
//      $                   nval(*)
//       DOUBLE PRECISION   a(*), afac(*), B(*), rwork(*),
//      $                   work(*), X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkgb tests Dgbtrf, -TRS, -RFS, and -CON
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
//          The number of values of nb contained in the vector nbval.
// \endverbatim
//
// \param[in] nbval
// \verbatim
//          nbval is intEGER array, dimension (nnb)
//          The values of the blocksize nb.
// \endverbatim
//
// \param[in] nns
// \verbatim
//          nns is intEGER
//          The number of values of nrhs contained in the vector nsval.
// \endverbatim
//
// \param[in] nsval
// \verbatim
//          nsval is intEGER array, dimension (nns)
//          The values of the number of right hand sides nrhs.
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
//          A is DOUBLE PRECISION array, dimension (LA)
// \endverbatim
//
// \param[in] LA
// \verbatim
//          LA is intEGER
//          The length of the array A.  LA >= (klmax+kuMAX+1)*nMAX
//          where klmax is the largest entry in the local array klval,
//                kuMAX is the largest entry in the local array kuval and
//                nMAX is the largest entry in the input array nval.
// \endverbatim
//
// \param[out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (lafac)
// \endverbatim
//
// \param[in] lafac
// \verbatim
//          lafac is intEGER
//          The length of the array afac. lafac >= (2*klMAX+kuMAX+1)*nMAX
//          where klmax is the largest entry in the local array klval,
//                kuMAX is the largest entry in the local array kuval and
//                nMAX is the largest entry in the input array nval.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (nMAX*nsmax)
//          where NSMAX is the largest entry in nsval.
// \endverbatim
//
// \param[out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (nMAX*nsmax)
// \endverbatim
//
// \param[out] xact
// \verbatim
//          xact is DOUBLE PRECISION array, dimension (nMAX*nsmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (nMAX*max(3,NSMAX,nMAX))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension
//                      (max(nMAX,2*nsmax))
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is intEGER array, dimension (2*nMAX)
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
func Dchkgb(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, a *[]float64, la *int, afac *[]float64, lafac *int, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	nbw := new(int)
	ntran := new(int)
	trfcon := new(bool)
	zerot := new(bool)
	dist := new(byte)
	norm := new(byte)
	trans := new(byte)
	_type := new(byte)
	xtype := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	i1 := new(int)
	i2 := new(int)
	ikl := new(int)
	iku := new(int)
	im := new(int)
	imat := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	ioff := new(int)
	irhs := new(int)
	itran := new(int)
	izero := new(int)
	j := new(int)
	k := new(int)
	kl := new(int)
	koff := new(int)
	ku := new(int)
	lda := new(int)
	ldafac := new(int)
	ldb := new(int)
	m := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nkl := new(int)
	nku := new(int)
	nrhs := new(int)
	nrun := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	anormi := new(float64)
	anormo := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	rcondi := new(float64)
	rcondo := new(float64)
	transs := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iseedy := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	klval := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	kuval := func() *[]int {
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
	common.infoc.lerr = new(bool)
	common.infoc.ok = new(bool)
	common.infoc.nunit = new(int)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new([]byte)
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
	(*zero) = 0.0e+0
	(*ntypes) = 8
	(*ntests) = 7
	(*nbw) = 4
	(*ntran) = 3
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
	nunit = common.infoc.nunit
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Data statements ..
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3], (*transs)[0], (*transs)[1], (*transs)[2] = 1988, 1989, 1990, 1991, 'N', 'T', 'C'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y :=[]byte("GB"); return &y }()
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
		Derrge(path, (nout))
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	//
	//     Initialize the first value for the lower and upper bandwidths.
	//
	(*klval)[0] = 0
	(*kuval)[0] = 0
	//
	//     Do for each value of M in mval
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		(*m) = (*(mval))[(*im)-1]
		//
		//        Set values to use for the lower bandwidth.
		//
		(*klval)[1] = (*m) + ((*m)+1)/4
		//
		//        klval( 2) = MAX( M-1, 0)
		//
		(*klval)[2] = (3*(*m) - 1) / 4
		(*klval)[3] = ((*m) + 1) / 4
		//
		//        Do for each value of N in nval
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			(*n) = (*(nval))[(*in)-1]
			(*xtype) = 'N'
			//
			//           Set values to use for the upper bandwidth.
			//
			(*kuval)[1] = (*n) + ((*n)+1)/4
			//
			//           kuval( 2) = MAX( N-1, 0)
			//
			(*kuval)[2] = (3*(*n) - 1) / 4
			(*kuval)[3] = ((*n) + 1) / 4
			//
			//           Set limits on the number of loop iterations.
			//
			(*nkl) = (Min((*m)+1, int(4)))
			if (*n) == 0 {
				(*nkl) = 2
			}
			(*nku) = (Min((*n)+1, int(4)))
			if (*m) == 0 {
				(*nku) = 2
			}
			(*nimat) = (*ntypes)
			if (*m) <= 0 || (*n) <= 0 {
				(*nimat) = 1
			}
			//
			for (*ikl) = 1; (*ikl) <= (*nkl); (*ikl)++ {
				//
				//              Do for kl = 0, (5*m+1)/4, (3M-1)/4, and (M+1)/4. This
				//              order makes it easier to skip redundant values for small
				//              values of M.
				//
				(*kl) = (*klval)[(*ikl)-1]
				for (*iku) = 1; (*iku) <= (*nku); (*iku)++ {
					//
					//                 Do for ku = 0, (5*n+1)/4, (3N-1)/4, and (N+1)/4. This
					//                 order makes it easier to skip redundant values for
					//                 small values of N.
					//
					(*ku) = (*kuval)[(*iku)-1]
					//
					//                 Check that A and afac are big enough to generate this
					//                 matrix.
					//
					(*lda) = (*kl) + (*ku) + 1
					(*ldafac) = 2*(*kl) + (*ku) + 1
					if ((*lda)*(*n)) > (*la) || ((*ldafac)*(*n)) > (*(lafac)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Alahd((nout), path)
						}
						if (*n)*((*kl)+(*ku)+1) > (*la) {
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" *** In Dchkgb, LA=%5d is too small for M=%5d, N=%5d, kl=%4d, ku=%4d\n ==> increase LA to at least %5d\n")
								return &y
							}(), (*la), (*m), (*n), (*kl), (*ku), (*n)*((*kl)+(*ku)+1))
							(*nerrs) = (*nerrs) + 1
						}
						if (*n)*(2*(*kl)+(*ku)+1) > (*(lafac)) {
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" *** In Dchkgb, lafac=%5d is too small for M=%5d, N=%5d, kl=%4d, ku=%4d\n ==> increase lafac to at least %5d\n")
								return &y
							}(), (*(lafac)), (*m), (*n), (*kl), (*ku), (*n)*(2*(*kl)+(*ku)+1))
							(*nerrs) = (*nerrs) + 1
						}
						goto Label130
					}
					//
					for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
						//
						//                    Do the tests only if dotype( imat) is true.
						//
						if !(*(dotype))[(*imat)-1] {
							goto Label120
						}
						//
						//                    Skip types 2, 3, or 4 if the matrix size is too
						//                    small.
						//
						(*zerot) = (*imat) >= 2 && (*imat) <= 4
						if (*zerot) && (*n) < (*imat)-1 {
							goto Label120
						}
						//
						if !(*zerot) || !(*(dotype))[0] {
							//
							//                       Set up parameters with Dlatb4 and generate a
							//                       test matrix with Dlatms.
							//
							Dlatb4(path, imat, m, n, _type, kl, ku, anorm, mode, cndnum, dist)
							//
							(*koff) = (MAX(1, (*ku)+2-(*n)))
							for (*i) = 1; (*i) <= (*koff)-1; (*i)++ {
								(*a)[(*i)-1] = (*zero)
								//Label20:
							}
							(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y }()
							Dlatms(m, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, func() *byte {y := byte('Z'); return &y }(), &((*a)[(*koff)-1]), lda, (work), info)
							//
							//                       Check the error code from Dlatms.
							//
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), m, n, kl, ku, -1, imat, nfail, nerrs, (nout))
								goto Label120
							}
						} else if (*izero) > 0 {
							//
							//                       Use the same matrix for types 3 and 4 as for
							//                       type 2 by copying back the zeroed out column.
							//
							Dcopy((*i2)-(*i1)+1, (b), func() *int {y := 1; return &y }(), &((*a)[(*ioff)+(*i1)-1]), func() *int {y := 1; return &y }())
						}
						//
						//                    For types 2, 3, and 4, zero one or more columns of
						//                    the matrix to test that info is returned correctly.
						//
						(*izero) = 0
						if *zerot {
							if (*imat) == 2 {
								(*izero) = 1
							} else if (*imat) == 3 {
								(*izero) = (Min((*m), (*n)))
							} else {
								(*izero) = Min((*m), (*n))/2 + 1
							}
							(*ioff) = ((*izero) - 1) * (*lda)
							if (*imat) < 4 {
								//
								//                          Store the column to be zeroed out in B.
								//
								(*i1) = (MAX(1, (*ku)+2-(*izero)))
								(*i2) = (Min((*kl)+(*ku)+1, (*ku)+1+((*m)-(*izero))))
								Dcopy((*i2)-(*i1)+1, &((*a)[(*ioff)+(*i1)-1]), func() *int {y := 1; return &y }(), (b), func() *int {y := 1; return &y }())
								//
								for (*i) = (*i1); (*i) <= (*i2); (*i)++ {
									(*a)[(*ioff)+(*i)-1] = (*zero)
									//Label30:
								}
							} else {
								for (*j) = (*izero); (*j) <= (*n); (*j)++ {
									for (*i) = MAX(1, (*ku)+2-(*j)); (*i) <= (Min((*kl)+(*ku)+1, (*ku)+1+((*m)-(*j)))); (*i)++ {
										(*a)[(*ioff)+(*i)-1] = (*zero)
										//Label40:
									}
									(*ioff) = (*ioff) + (*lda)
									//Label50:
								}
							}
						}
						//
						//                    These lines, if used in place of the calls in the
						//                    loop over inb, cause the code to bomb on a Sun
						//                    SPARCstation.
						//
						//                     anormo = dlangb( 'O', n, kl, ku, a, lda, rwork)
						//                     anormi = dlangb( 'I', n, kl, ku, a, lda, rwork)
						//
						//                    Do for each blocksize in nbval
						//
						for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
							(*nb) = (*(nbval))[(*inb)-1]
							Xlaenv(func() *int {y := 1; return &y }(), nb)
							//
							//                       Compute the LU factorization of the band matrix.
							//
							if (*m) > 0 && (*n) > 0 {
								Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (*kl)+(*ku)+1, n, (a), lda, &((*(afac))[(*kl)+0]), ldafac)
							}
							(*srnamt) = *func() *[]byte {y :=[]byte("Dgbtrf"); return &y }()
							Dgbtrf(m, n, kl, ku, (afac), ldafac, (iwork), info)
							//
							//                       Check error code from Dgbtrf.
							//
							if (*info) != (*izero) {
								Alaerh(path, func() *[]byte {y :=[]byte("Dgbtrf"); return &y }(), info, izero, func() *byte {y := byte(' '); return &y }(), m, n, kl, ku, nb, imat, nfail, nerrs, (nout))
							}
							(*trfcon) = false
							//
							//+    TEST 1
							//                       Re_construct matrix from factors and compute
							//                       residual.
							//
							Dgbt01(m, n, kl, ku, (a), lda, (afac), ldafac, (iwork), (work), &((*result)[0]))
							//
							//                       Print information about the tests so far that
							//                       did not pass the threshold.
							//
							if (*result)[0] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Alahd((nout), path)
								}
								WRITE((*(nout)), *func() *[]byte {
									y :=[]byte(" M =%5d, N =%5d, kl=%5d, ku=%5d, nb =%4d, type %1d, test(%1d)=%12.5f\n")
									return &y
								}(), (*m), (*n), (*kl), (*ku), (*nb), (*imat), 1, (*result)[0])
								(*nfail) = (*nfail) + 1
							}
							(*nrun) = (*nrun) + 1
							//
							//                       Skip the remaining tests if this is not the
							//                       first block size or if M .ne. N.
							//
							if (*inb) > 1 || (*m) != (*n) {
								goto Label110
							}
							//
							(*anormo) = (*dlangb(func() *byte {y := byte('O'); return &y }(), n, kl, ku, (a), lda, (rwork)))
							(*anormi) = (*dlangb(func() *byte {y := byte('I'); return &y }(), n, kl, ku, (a), lda, (rwork)))
							//
							if (*info) == 0 {
								//
								//                          Form the inverse of A so we can get a good
								//                          estimate of cndnum = norm(a) * norm(inv(a)).
								//
								(*ldb) = (MAX(1, (*n)))
								Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n, n, zero, one, (work), ldb)
								(*srnamt) = *func() *[]byte {y :=[]byte("Dgbtrs"); return &y }()
								Dgbtrs(func() *[]byte {y :=[]byte("No transpose"); return &y }(), n, kl, ku, n, (afac), ldafac, (iwork), (work), ldb, info)
								//
								//                          Compute the 1-norm condition number of A.
								//
								(*ainvnm) = (*Dlange(func() *byte {y := byte('O'); return &y }(), n, n, (work), ldb, (rwork)))
								if (*anormo) <= (*zero) || (*ainvnm) <= (*zero) {
									(*rcondo) = (*one)
								} else {
									(*rcondo) = ((*one) / (*anormo)) / (*ainvnm)
								}
								//
								//                          Compute the infinity-norm condition number of
								//                          A.
								//
								(*ainvnm) = (*Dlange(func() *byte {y := byte('I'); return &y }(), n, n, (work), ldb, (rwork)))
								if (*anormi) <= (*zero) || (*ainvnm) <= (*zero) {
									(*rcondi) = (*one)
								} else {
									(*rcondi) = ((*one) / (*anormi)) / (*ainvnm)
								}
							} else {
								//
								//                          Do only the condition estimate if info.NE.0.
								//
								(*trfcon) = true
								(*rcondo) = (*zero)
								(*rcondi) = (*zero)
							}
							//
							//                       Skip the solve tests if the matrix is singular.
							//
							if *trfcon {
								goto Label90
							}
							//
							for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
								(*nrhs) = (*(nsval))[(*irhs)-1]
								(*xtype) = 'N'
								//
								for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
									(*trans) = (*transs)[(*itran)-1]
									if (*itran) == 1 {
										(*rcondc) = (*rcondo)
										(*norm) = 'O'
									} else {
										(*rcondc) = (*rcondi)
										(*norm) = 'I'
									}
									//
									//+    TEST 2:
									//                             Solve and compute residual for A * X = B.
									//
									(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y }()
									Dlarhs(path, xtype, func() *byte {y := byte(' '); return &y }(), trans, n, n, kl, ku, nrhs, (a), lda, (xact), ldb, (b), ldb, iseed, info)
									(*xtype) = 'C'
									Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), ldb, (x), ldb)
									//
									(*srnamt) = *func() *[]byte {y :=[]byte("Dgbtrs"); return &y }()
									Dgbtrs(trans, n, kl, ku, nrhs, (afac), ldafac, (iwork), (x), ldb, info)
									//
									//                             Check error code from Dgbtrs.
									//
									if (*info) != 0 {
										Alaerh(path, func() *[]byte {y :=[]byte("Dgbtrs"); return &y }(), info, func() *int {y := 0; return &y }(), trans, n, n, kl, ku, -1, imat, nfail, nerrs, (nout))
									}
									//
									Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), ldb, (work), ldb)
									Dgbt02(trans, m, n, kl, ku, nrhs, (a), lda, (x), ldb, (work), ldb, &((*result)[1]))
									//
									//+    TEST 3:
									//                             Check solution from generated exact
									//                             solution.
									//
									Dget04(n, nrhs, (x), ldb, (xact), ldb, rcondc, &((*result)[2]))
									//
									//+    TESts 4, 5, 6:
									//                             Use iterative refinement to improve the
									//                             solution.
									//
									(*srnamt) = *func() *[]byte {y :=[]byte("Dgbrfs"); return &y }()
									Dgbrfs(trans, n, kl, ku, nrhs, (a), lda, (afac), ldafac, (iwork), (b), ldb, (x), ldb, (rwork), &((*(rwork))[(*nrhs)+0]), (work), &((*(iwork))[(*n)+0]), info)
									//
									//                             Check error code from Dgbrfs.
									//
									if (*info) != 0 {
										Alaerh(path, func() *[]byte {y :=[]byte("Dgbrfs"); return &y }(), info, func() *int {y := 0; return &y }(), trans, n, n, kl, ku, nrhs, imat, nfail, nerrs, (nout))
									}
									//
									Dget04(n, nrhs, (x), ldb, (xact), ldb, rcondc, &((*result)[3]))
									Dgbt05(trans, n, kl, ku, nrhs, (a), lda, (b), ldb, (x), ldb, (xact), ldb, (rwork), &((*(rwork))[(*nrhs)+0]), &((*result)[4]))
									for (*k) = 2; (*k) <= 6; (*k)++ {
										if (*result)[(*k)-1] >= (*(thresh)) {
											if (*nfail) == 0 && (*nerrs) == 0 {
												Alahd((nout), path)
											}
											WRITE((*(nout)), *func() *[]byte {
												y :=[]byte(" trans='%c', N=%5d, kl=%5d, ku=%5d, nrhs=%3d, type %1d, test(%1d)=%12.5f\n")
												return &y
											}(), (*trans), (*n), (*kl), (*ku), (*nrhs), (*imat), (*k), (*result)[(*k)-1])
											(*nfail) = (*nfail) + 1
										}
										//Label60:
									}
									(*nrun) = (*nrun) + 5
									//Label70:
								}
								//Label80:
							}
							//
							//+    TEST 7:
							//                          Get an estimate of rcond = 1/cndnum.
							//
						Label90:
							;
							for (*itran) = 1; (*itran) <= 2; (*itran)++ {
								if (*itran) == 1 {
									(*anorm) = (*anormo)
									(*rcondc) = (*rcondo)
									(*norm) = 'O'
								} else {
									(*anorm) = (*anormi)
									(*rcondc) = (*rcondi)
									(*norm) = 'I'
								}
								(*srnamt) = *func() *[]byte {y :=[]byte("Dgbcon"); return &y }()
								Dgbcon(norm, n, kl, ku, (afac), ldafac, (iwork), anorm, rcond, (work), &((*(iwork))[(*n)+0]), info)
								//
								//                             Check error code from Dgbcon.
								//
								if (*info) != 0 {
									Alaerh(path, func() *[]byte {y :=[]byte("Dgbcon"); return &y }(), info, func() *int {y := 0; return &y }(), norm, n, n, kl, ku, -1, imat, nfail, nerrs, (nout))
								}
								//
								(*result)[6] = (*Dget06(rcond, rcondc))
								//
								//                          Print information about the tests that did
								//                          not pass the threshold.
								//
								if (*result)[6] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Alahd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte {
										y :=[]byte(" norm ='%c', N=%5d, kl=%5d, ku=%5d,           type %1d, test(%1d)=%12.5f\n")
										return &y
									}(), (*norm), (*n), (*kl), (*ku), (*imat), 7, (*result)[6])
									(*nfail) = (*nfail) + 1
								}
								(*nrun) = (*nrun) + 1
								//Label100:
							}
							//
						Label110:
						}
					Label120:
					}
				Label130:
				}
				//Label140:
			}
			//Label150:
		}
		//Label160:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	//
	return
	//
	//     End of Dchkgb
	//
}
