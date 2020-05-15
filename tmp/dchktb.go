package goblas

import 

// Dchktb tests Dtbtrs, -RFS, and -CON, and Dlatbs.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchktb( dotype, nn, nval, nns, nsval, thresh, tsterr,
//                          nmax, AB, ainv, b, x, xact, work, rwork, iwork,
//                          nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nmax, nn, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), nsval(*), nval(*)
//       DOUBLE PRECISION   AB(*), ainv(*), B(*), rwork(*),
//      $                   work(*), X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchktb tests Dtbtrs, -RFS, and -CON, and Dlatbs.
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
//          The leading dimension of the work arrays.
//          nmax >= the maximum value of N in nval.
// \endverbatim
//
// \param[out] AB
// \verbatim
//          AB is DOUBLE PRECISION array, dimension (nmax*nmax)
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
//                      (max(nmax,2*nsmax))
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
func Dchktb(dotype *[]bool, nn *int, nval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, nmax *int, ab *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	ntype1 := new(int)
	ntypes := new(int)
	ntests := new(int)
	ntran := new(int)
	one := new(float64)
	zero := new(float64)
	diag := new(byte)
	norm := new(byte)
	trans := new(byte)
	uplo := new(byte)
	xtype := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	idiag := new(int)
	ik := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	irhs := new(int)
	itran := new(int)
	iuplo := new(int)
	j := new(int)
	k := new(int)
	kd := new(int)
	lda := new(int)
	ldab := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nimat2 := new(int)
	nk := new(int)
	nrhs := new(int)
	nrun := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	rcondi := new(float64)
	rcondo := new(float64)
	scale := new(float64)
	transs := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	uplos := func() *[]byte {
		arr := make([]byte, 2)
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
	(*ntype1) = 9
	(*ntypes) = 17
	(*ntests) = 8
	(*ntran) = 3
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
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
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Data statements ..
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3] = 1988, 1989, 1990, 1991
	(*uplos)[0], (*uplos)[1], (*transs)[0], (*transs)[1], (*transs)[2] = 'U', 'L', 'N', 'T', 'C'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte{y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte{y := []byte("TB"); return &y }()
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
		Derrtr(path, (nout))
	}
	(*infot) = 0
	//
	for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
		//
		//        Do for each value of N in nval
		//
		(*n) = (*(nval))[(*in)-1]
		(*lda) = (MAX(1, (*n)))
		(*xtype) = 'N'
		(*nimat) = (*ntype1)
		(*nimat2) = (*ntypes)
		if (*n) <= 0 {
			(*nimat) = 1
			(*nimat2) = (*ntype1) + 1
		}
		//
		(*nk) = (Min((*n)+1, int(4)))
		for (*ik) = 1; (*ik) <= (*nk); (*ik)++ {
			//
			//           Do for kd = 0, n, (3N-1)/4, and (N+1)/4. This order makes
			//           it easier to skip redundant values for small values of N.
			//
			if (*ik) == 1 {
				(*kd) = 0
			} else if (*ik) == 2 {
				(*kd) = (MAX((*n), 0))
			} else if (*ik) == 3 {
				(*kd) = (3*(*n) - 1) / 4
			} else if (*ik) == 4 {
				(*kd) = ((*n) + 1) / 4
			}
			(*ldab) = (*kd) + 1
			//
			for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
				//
				//              Do the tests only if dotype( imat) is true.
				//
				if !(*(dotype))[(*imat)-1] {
					goto Label90
				}
				//
				for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
					//
					//                 Do first for uplo = 'U', then for uplo = 'L'
					//
					(*uplo) = (*uplos)[(*iuplo)-1]
					//
					//                 Call Dlattb to generate a triangular test matrix.
					//
					(*srnamt) = *func() *[]byte{y := []byte("Dlattb"); return &y }()
					Dlattb(imat, uplo, func() *[]byte{y := []byte("No transpose"); return &y }(), diag, iseed, n, kd, (ab), ldab, (x), (work), info)
					//
					//                 Set idiag = 1 for non-unit matrices, 2 for unit.
					//
					if blas.Lsame(diag, func() *byte{y := byte('N'); return &y }()) {
						(*idiag) = 1
					} else {
						(*idiag) = 2
					}
					//
					//                 Form the inverse of A so we can get a good estimate
					//                 of rcondc = 1/(norm(a) * norm(inv(a))).
					//
					Dlaset(func() *[]byte{y := []byte("Full"); return &y }(), n, n, zero, one, (ainv), lda)
					if blas.Lsame(uplo, func() *byte{y := byte('U'); return &y }()) {
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							Dtbsv(uplo, func() *[]byte{y := []byte("No transpose"); return &y }(), diag, j, kd, (ab), ldab, &((*(ainv))[((*j)-1)*(*lda)+0]), func() *int{y := 1; return &y }())
						//Label20:
						}
					} else {
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							Dtbsv(uplo, func() *[]byte{y := []byte("No transpose"); return &y }(), diag, (*n)-(*j)+1, kd, &((*(ab))[((*j)-1)*(*ldab)+0]), ldab, &((*(ainv))[((*j)-1)*(*lda)+(*j)-1]), func() *int{y := 1; return &y }())
						//Label30:
						}
					}
					//
					//                 Compute the 1-norm condition number of A.
					//
					(*anorm) = (*Dlantb(func() *byte{y := byte('1'); return &y }(), uplo, diag, n, kd, (ab), ldab, (rwork)))
					(*ainvnm) = (*Dlantr(func() *byte{y := byte('1'); return &y }(), uplo, diag, n, n, (ainv), lda, (rwork)))
					if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
						(*rcondo) = (*one)
					} else {
						(*rcondo) = ((*one) / (*anorm)) / (*ainvnm)
					}
					//
					//                 Compute the infinity-norm condition number of A.
					//
					(*anorm) = (*Dlantb(func() *byte{y := byte('I'); return &y }(), uplo, diag, n, kd, (ab), ldab, (rwork)))
					(*ainvnm) = (*Dlantr(func() *byte{y := byte('I'); return &y }(), uplo, diag, n, n, (ainv), lda, (rwork)))
					if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
						(*rcondi) = (*one)
					} else {
						(*rcondi) = ((*one) / (*anorm)) / (*ainvnm)
					}
					//
					for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
						(*nrhs) = (*(nsval))[(*irhs)-1]
						(*xtype) = 'N'
						//
						for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
							//
							//                    Do for op(a) = a, A**T, or A**H.
							//
							(*trans) = (*transs)[(*itran)-1]
							if (*itran) == 1 {
								(*norm) = 'O'
								(*rcondc) = (*rcondo)
							} else {
								(*norm) = 'I'
								(*rcondc) = (*rcondi)
							}
							//
							//+    TEST 1
							//                    Solve and compute residual for op(a)*x = b.
							//
							(*srnamt) = *func() *[]byte{y := []byte("Dlarhs"); return &y }()
							Dlarhs(path, xtype, uplo, trans, n, n, kd, idiag, nrhs, (ab), ldab, (xact), lda, (b), lda, iseed, info)
							(*xtype) = 'C'
							Dlacpy(func() *[]byte{y := []byte("Full"); return &y }(), n, nrhs, (b), lda, (x), lda)
							//
							(*srnamt) = *func() *[]byte{y := []byte("Dtbtrs"); return &y }()
							Dtbtrs(uplo, trans, diag, n, kd, nrhs, (ab), ldab, (x), lda, info)
							//
							//                    Check error code from Dtbtrs.
							//
							if (*info) != 0 {
								Alaerh ((*path), "Dtbtrs", (*info), 0, (*uplo) append ( append ([] byte { }, append ( append ([] byte { },), (*trans))), (*diag)), (*N), (*N), (*kd), (*kd), (*nrhs), (*imat), (*nfail), (*nerrs), (*nout))
							}
							//
							Dtbt02(uplo, trans, diag, n, kd, nrhs, (ab), ldab, (x), lda, (b), lda, (work), &((*result)[0]))
							//
							//+    TEST 2
							//                    Check solution from generated exact solution.
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[1]))
							//
							//+    TESts 3, 4, and 5
							//                    Use iterative refinement to improve the solution
							//                    and compute error bounds.
							//
							(*srnamt) = *func() *[]byte{y := []byte("Dtbrfs"); return &y }()
							Dtbrfs(uplo, trans, diag, n, kd, nrhs, (ab), ldab, (b), lda, (x), lda, (rwork), &((*(rwork))[(*nrhs)+0]), (work), (iwork), info)
							//
							//                    Check error code from Dtbrfs.
							//
							if (*info) != 0 {
								Alaerh ((*path), "Dtbrfs", (*info), 0, (*uplo) append ( append ([] byte { }, append ( append ([] byte { },), (*trans))), (*diag)), (*N), (*N), (*kd), (*kd), (*nrhs), (*imat), (*nfail), (*nerrs), (*nout))
							}
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[2]))
							Dtbt05(uplo, trans, diag, n, kd, nrhs, (ab), ldab, (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*nrhs)+0]), &((*result)[3]))
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = 1; (*k) <= 5; (*k)++ {
								if (*result)[(*k)-1] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Alahd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte{y := []byte(" uplo='%c', trans='%c',        diag='%c', N=%5d, kd=%5d, nrhs=%5d, type %2d, test(%2d)=%12.5f\n"); return &y }(), (*uplo), (*trans), (*diag), (*n), (*kd), (*nrhs), (*imat), (*k), (*result)[(*k)-1])
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
					//+    TEST 6
					//                    Get an estimate of rcond = 1/cndnum.
					//
					for (*itran) = 1; (*itran) <= 2; (*itran)++ {
						if (*itran) == 1 {
							(*norm) = 'O'
							(*rcondc) = (*rcondo)
						} else {
							(*norm) = 'I'
							(*rcondc) = (*rcondi)
						}
						(*srnamt) = *func() *[]byte{y := []byte("Dtbcon"); return &y }()
						Dtbcon(NORM, uplo, diag, n, kd, (ab), ldab, rcond, (work), (iwork), info)
						//
						//                    Check error code from Dtbcon.
						//
						if (*info) != 0 {
							Alaerh ((*path), "Dtbcon", (*info), 0, (*NORM) append ( append ([] byte { }, append ( append ([] byte { },), (*uplo))), (*diag)), (*N), (*N), (*kd), (*kd), - 1, (*imat), (*nfail), (*nerrs), (*nout))
						}
						//
						Dtbt06(rcond, rcondc, uplo, diag, n, kd, (ab), ldab, (rwork), &((*result)[5]))
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						if (*result)[5] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c',%5d,%5d,  ...), type %2d, test(%2d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dtbcon"); return &y }(), (*norm), (*uplo), (*diag), (*n), (*kd), (*imat), 6, (*result)[5])
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 1
					//Label70:
					}
				//Label80:
				}
			Label90:
			}
			//
			//           Use pathological test matrices to test Dlatbs.
			//
			for (*imat) = (*ntype1) + 1; (*imat) <= (*nimat2); (*imat)++ {
				//
				//              Do the tests only if dotype( imat) is true.
				//
				if !(*(dotype))[(*imat)-1] {
					goto Label120
				}
				//
				for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
					//
					//                 Do first for uplo = 'U', then for uplo = 'L'
					//
					(*uplo) = (*uplos)[(*iuplo)-1]
					for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
						//
						//                    Do for op(a) = a, A**T, and A**H.
						//
						(*trans) = (*transs)[(*itran)-1]
						//
						//                    Call Dlattb to generate a triangular test matrix.
						//
						(*srnamt) = *func() *[]byte{y := []byte("Dlattb"); return &y }()
						Dlattb(imat, uplo, trans, diag, iseed, n, kd, (ab), ldab, (x), (work), info)
						//
						//+    TEST 7
						//                    Solve the system op(a)*x = b
						//
						(*srnamt) = *func() *[]byte{y := []byte("Dlatbs"); return &y }()
						Dcopy(n, (x), func() *int{y := 1; return &y }(), (b), func() *int{y := 1; return &y }())
						Dlatbs(uplo, trans, diag, func() *byte{y := byte('N'); return &y }(), n, kd, (ab), ldab, (b), scale, (rwork), info)
						//
						//                    Check error code from Dlatbs.
						//
						if (*info) != 0 {
							Alaerh ((*path), "Dlatbs", (*info), 0, append ( append ([] byte { }, (*uplo)), (*trans)) append ( append ([] byte { }, append ( append ([] byte { },), (*diag))), "N"), (*N), (*N), (*kd), (*kd), - 1, (*imat), (*nfail), (*nerrs), (*nout))
						}
						//
						Dtbt03(uplo, trans, diag, n, kd, func() *int{y := 1; return &y }(), (ab), ldab, scale, (rwork), one, (b), lda, (x), lda, (work), &((*result)[6]))
						//
						//+    TEST 8
						//                    Solve op(a)*x = b again with NOrmin = 'Y'.
						//
						Dcopy(n, (x), func() *int{y := 1; return &y }(), (b), func() *int{y := 1; return &y }())
						Dlatbs(uplo, trans, diag, func() *byte{y := byte('Y'); return &y }(), n, kd, (ab), ldab, (b), scale, (rwork), info)
						//
						//                    Check error code from Dlatbs.
						//
						if (*info) != 0 {
							Alaerh ((*path), "Dlatbs", (*info), 0, append ( append ([] byte { }, (*uplo)), (*trans)) append ( append ([] byte { }, append ( append ([] byte { },), (*diag))), "Y"), (*N), (*N), (*kd), (*kd), - 1, (*imat), (*nfail), (*nerrs), (*nout))
						}
						//
						Dtbt03(uplo, trans, diag, n, kd, func() *int{y := 1; return &y }(), (ab), ldab, scale, (rwork), one, (b), lda, (x), lda, (work), &((*result)[7]))
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						if (*result)[6] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c', '%c',%5d,%5d, ... ),  type %2d, test(%1d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dlatbs"); return &y }(), (*uplo), (*trans), (*diag), 'N', (*n), (*kd), (*imat), 7, (*result)[6])
							(*nfail) = (*nfail) + 1
						}
						if (*result)[7] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c', '%c',%5d,%5d, ... ),  type %2d, test(%1d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dlatbs"); return &y }(), (*uplo), (*trans), (*diag), 'Y', (*n), (*kd), (*imat), 8, (*result)[7])
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 2
					//Label100:
					}
				//Label110:
				}
			Label120:
			}
		//Label130:
		}
	//Label140:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchktb
	//
}
