package goblas

import 

// Dchkge tests Dgetrf, -tri, -TRS, -RFS, and -CON.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkge( dotype, nm, mval, nn, nval, nnb, nbval, nns,
//                          nsval, thresh, tsterr, nmax, a, afac, ainv, b,
//                          x, xact, work, rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nm, nmax, nn, nnb, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), mval(*), nbval(*), nsval(*),
//      $                   nval(*)
//       DOUBLE PRECISION   a(*), afac(*), ainv(*), B(*),
//      $                   rwork(*), work(*), X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkge tests Dgetrf, -tri, -TRS, -RFS, and -CON.
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
//          nbval is intEGER array, dimension (nbval)
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
// \param[out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (nmax*nsmax)
//          where NSMAX is the largest entry in nsval.
// \endverbatim
//
// \param[out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (nmax*nsmax)
// \endverbatim
//
// \param[out] xact
// \verbatim
//          xact is DOUBLE PRECISION array, dimension (nmax*nsmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (nmax*max(3,NSMAX))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension
//                      (max(2*nmax,2*nsmax+Nwork))
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is intEGER array, dimension (2*nmax)
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
func Dchkge(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
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
	im := new(int)
	imat := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	ioff := new(int)
	irhs := new(int)
	itran := new(int)
	izero := new(int)
	k := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	lwork := new(int)
	m := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nrhs := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	anormi := new(float64)
	anormo := new(float64)
	cndnum := new(float64)
	dummy := new(float64)
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
	result := func() *[]float64 {
		arr := make([]float64, 8)
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
	(*ntypes) = 11
	(*ntests) = 8
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
	(*path)[1] = *func() *[]byte {y :=[]byte("GE"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-(1)] = (*iseedy)[(*i)-(1)]
		//Label10:
	}
	//
	//     Test the error exits
	//
	Xlaenv(func() *int {y := 1; return &y }(), func() *int {y := 1; return &y }())
	if *(tsterr) {
		Derrge(path, (nout))
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	//
	//     Do for each value of M in mval
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		(*m) = (*(mval))[(*im)-(1)]
		(*lda) = (MAX(1, (*m)))
		//
		//        Do for each value of N in nval
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			(*n) = (*(nval))[(*in)-(1)]
			(*xtype) = 'N'
			(*nimat) = (*ntypes)
			if (*m) <= 0 || (*n) <= 0 {
				(*nimat) = 1
			}
			//
			for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
				//
				//              Do the tests only if dotype( imat) is true.
				//
				if !(*(dotype))[(*imat)-(1)] {
					goto Label100
				}
				//
				//              Skip types 5, 6, or 7 if the matrix size is too small.
				//
				(*zerot) = (*imat) >= 5 && (*imat) <= 7
				if (*zerot) && (*n) < (*imat)-4 {
					goto Label100
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
					goto Label100
				}
				//
				//              For types 5-7, zero one or more columns of the matrix to
				//              test that info is returned correctly.
				//
				if *zerot {
					if (*imat) == 5 {
						(*izero) = 1
					} else if (*imat) == 6 {
						(*izero) = (Min((*m), (*n)))
					} else {
						(*izero) = Min((*m), (*n))/2 + 1
					}
					(*ioff) = ((*izero) - 1) * (*lda)
					if (*imat) < 7 {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*(a))[(*ioff)+(*i)-(1)] = (*zero)
							//Label20:
						}
					} else {
						Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), m, (*n)-(*izero)+1, zero, zero, &((*(a))[(*ioff)+0]), lda)
					}
				} else {
					(*izero) = 0
				}
				//
				//              These lines, if used in place of the calls in the DO 60
				//              loop, cause the code to bomb on a Sun SPARCstation.
				//
				//               anormo = Dlange( 'O', m, n, a, lda, rwork)
				//               anormi = Dlange( 'I', m, n, a, lda, rwork)
				//
				//              Do for each blocksize in nbval
				//
				for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
					(*nb) = (*(nbval))[(*inb)-(1)]
					Xlaenv(func() *int {y := 1; return &y }(), nb)
					//
					//                 Compute the LU factorization of the matrix.
					//
					Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (a), lda, (afac), lda)
					(*srnamt) = *func() *[]byte {y :=[]byte("Dgetrf"); return &y }()
					Dgetrf(m, n, (afac), lda, (iwork), info)
					//
					//                 Check error code from Dgetrf.
					//
					if (*info) != (*izero) {
						Alaerh(path, func() *[]byte {y :=[]byte("Dgetrf"); return &y }(), info, izero, func() *byte {y := byte(' '); return &y }(), m, n, -1, -1, nb, imat, nfail, nerrs, (nout))
					}
					(*trfcon) = false
					//
					//+    TEST 1
					//                 Re_construct matrix from factors and compute residual.
					//
					Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, (afac), lda, (ainv), lda)
					Dget01(m, n, (a), lda, (ainv), lda, (iwork), (rwork), &((*result)[0]))
					(*nt) = 1
					//
					//+    TEST 2
					//                 Form the inverse if the factorization was successful
					//                 and compute the residual.
					//
					if (*m) == (*n) && (*info) == 0 {
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, n, (afac), lda, (ainv), lda)
						(*srnamt) = *func() *[]byte {y :=[]byte("Dgetri"); return &y }()
						(*nrhs) = (*(nsval))[0]
						(*lwork) = (*(nmax)) * MAX(3, (*nrhs))
						Dgetri(n, (ainv), lda, (iwork), (work), lwork, info)
						//
						//                    Check error code from Dgetri.
						//
						if (*info) != 0 {
							Alaerh(path, func() *[]byte {y :=[]byte("Dgetri"); return &y }(), info, func() *int {y := 0; return &y }(), func() *byte {y := byte(' '); return &y }(), n, n, -1, -1, nb, imat, nfail, nerrs, (nout))
						}
						//
						//                    Compute the residual for the matrix times its
						//                    inverse.  Also compute the 1-norm condition number
						//                    of A.
						//
						Dget03(n, (a), lda, (ainv), lda, (work), lda, (rwork), rcondo, &((*result)[1]))
						(*anormo) = (*Dlange(func() *byte {y := byte('O'); return &y }(), m, n, (a), lda, (rwork)))
						//
						//                    Compute the infinity-norm condition number of A.
						//
						(*anormi) = (*Dlange(func() *byte {y := byte('I'); return &y }(), m, n, (a), lda, (rwork)))
						(*ainvnm) = (*Dlange(func() *byte {y := byte('I'); return &y }(), n, n, (ainv), lda, (rwork)))
						if (*anormi) <= (*zero) || (*ainvnm) <= (*zero) {
							(*rcondi) = (*one)
						} else {
							(*rcondi) = ((*one) / (*anormi)) / (*ainvnm)
						}
						(*nt) = 2
					} else {
						//
						//                    Do only the condition estimate if info > 0.
						//
						(*trfcon) = true
						(*anormo) = (*Dlange(func() *byte {y := byte('O'); return &y }(), m, n, (a), lda, (rwork)))
						(*anormi) = (*Dlange(func() *byte {y := byte('I'); return &y }(), m, n, (a), lda, (rwork)))
						(*rcondo) = (*zero)
						(*rcondi) = (*zero)
					}
					//
					//                 Print information about the tests so far that did not
					//                 pass the threshold.
					//
					for (*k) = 1; (*k) <= (*nt); (*k)++ {
						if (*result)[(*k)-(1)] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {y :=[]byte(" M = %5d, N =%5d, nb =%4d, type %2d, test(%2d) =%12.5f\n"); return &y }(), (*m), (*n), (*nb), (*imat), (*k), (*result)[(*k)-(1)])
							(*nfail) = (*nfail) + 1
						}
						//Label30:
					}
					(*nrun) = (*nrun) + (*nt)
					//
					//                 Skip the remaining tests if this is not the first
					//                 block size or if M .ne. N.  Skip the solve tests if
					//                 the matrix is singular.
					//
					if (*inb) > 1 || (*m) != (*n) {
						goto Label90
					}
					if *trfcon {
						goto Label70
					}
					//
					for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
						(*nrhs) = (*(nsval))[(*irhs)-(1)]
						(*xtype) = 'N'
						//
						for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
							(*trans) = (*transs)[(*itran)-(1)]
							if (*itran) == 1 {
								(*rcondc) = (*rcondo)
							} else {
								(*rcondc) = (*rcondi)
							}
							//
							//+    TEST 3
							//                       Solve and compute residual for A * X = B.
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y }()
							Dlarhs(path, xtype, func() *byte {y := byte(' '); return &y }(), trans, n, n, kl, ku, nrhs, (a), lda, (xact), lda, (b), lda, iseed, info)
							(*xtype) = 'C'
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), lda, (x), lda)
							(*srnamt) = *func() *[]byte {y :=[]byte("Dgetrs"); return &y }()
							Dgetrs(trans, n, nrhs, (afac), lda, (iwork), (x), lda, info)
							//
							//                       Check error code from Dgetrs.
							//
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("Dgetrs"); return &y }(), info, func() *int {y := 0; return &y }(), trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, (nout))
							}
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), lda, (work), lda)
							Dget02(trans, n, n, nrhs, (a), lda, (x), lda, (work), lda, (rwork), &((*result)[2]))
							//
							//+    TEST 4
							//                       Check solution from generated exact solution.
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[3]))
							//
							//+    TESts 5, 6, and 7
							//                       Use iterative refinement to improve the
							//                       solution.
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("Dgerfs"); return &y }()
							Dgerfs(trans, n, nrhs, (a), lda, (afac), lda, (iwork), (b), lda, (x), lda, (rwork), &((*(rwork))[(*nrhs)+0]), (work), &((*(iwork))[(*n)+0]), info)
							//
							//                       Check error code from Dgerfs.
							//
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("Dgerfs"); return &y }(), info, func() *int {y := 0; return &y }(), trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, (nout))
							}
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[4]))
							Dget07(trans, n, nrhs, (a), lda, (b), lda, (x), lda, (xact), lda, (rwork), &(true), &((*(rwork))[(*nrhs)+0]), &((*result)[5]))
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = 3; (*k) <= 7; (*k)++ {
								if (*result)[(*k)-(1)] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Alahd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte {y :=[]byte(" trans='%c', N =%5d, nrhs=%3d, type %2d, test(%2d) =%12.5f\n"); return &y }(), (*trans), (*n), (*nrhs), (*imat), (*k), (*result)[(*k)-(1)])
									(*nfail) = (*nfail) + 1
								}
								//Label40:
							}
							(*nrun) = (*nrun) + 5
							//Label50:
						}
						//Label60:
					}
					//
					//+    TEST 8
					//                    Get an estimate of rcond = 1/cndnum.
					//
				Label70:
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
						(*srnamt) = *func() *[]byte {y :=[]byte("Dgecon"); return &y }()
						Dgecon(norm, n, (afac), lda, anorm, rcond, (work), &((*(iwork))[(*n)+0]), info)
						//
						//                       Check error code from Dgecon.
						//
						if (*info) != 0 {
							Alaerh(path, func() *[]byte {y :=[]byte("Dgecon"); return &y }(), info, func() *int {y := 0; return &y }(), norm, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
						}
						//
						//                       This line is needed on a Sun SPARCstation.
						//
						(*dummy) = (*rcond)
						//
						(*result)[7] = (*Dget06(rcond, rcondc))
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						if (*result)[7] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {y :=[]byte(" norm ='%c', N =%5d,           type %2d, test(%2d) =%12.5f\n"); return &y }(), (*norm), (*n), (*imat), 8, (*result)[7])
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 1
						//Label80:
					}
				Label90:
				}
			Label100:
			}
			//Label110:
		}
		//Label120:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchkge
	//
}
