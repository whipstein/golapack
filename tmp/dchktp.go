package goblas

import 

// Dchktp tests Dtptri, -TRS, -RFS, and -CON, and Dlatps
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchktp( dotype, nn, nval, nns, nsval, thresh, tsterr,
//                          nmax, AP, ainvp, b, x, xact, work, rwork,
//                          iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nmax, nn, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), nsval(*), nval(*)
//       DOUBLE PRECISION   ainvp(*), AP(*), B(*), rwork(*),
//      $                   work(*), X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchktp tests Dtptri, -TRS, -RFS, and -CON, and Dlatps
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
//          The leading dimension of the work arrays.  nmax >= the
//          maximumm value of N in nval.
// \endverbatim
//
// \param[out] AP
// \verbatim
//          AP is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
// \endverbatim
//
// \param[out] ainvp
// \verbatim
//          ainvp is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
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
// \param[out] iwork
// \verbatim
//          iwork is intEGER array, dimension (nmax)
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension
//                      (max(nmax,2*nsmax))
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
func Dchktp(dotype *[]bool, nn *int, nval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, nmax *int, ap *[]float64, ainvp *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
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
	imat := new(int)
	in := new(int)
	info := new(int)
	irhs := new(int)
	itran := new(int)
	iuplo := new(int)
	k := new(int)
	lap := new(int)
	lda := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
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
		arr := make([]float64, 9)
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
	(*ntype1) = 10
	(*ntypes) = 18
	(*ntests) = 9
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
	(*path)[1] = *func() *[]byte{y := []byte("TP"); return &y }()
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
		(*lap) = (*lda) * ((*lda) + 1) / 2
		(*xtype) = 'N'
		//
		for (*imat) = 1; (*imat) <= (*ntype1); (*imat)++ {
			//
			//           Do the tests only if dotype( imat) is true.
			//
			if !(*(dotype))[(*imat)-1] {
				goto Label70
			}
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				//
				//              Do first for uplo = 'U', then for uplo = 'L'
				//
				(*uplo) = (*uplos)[(*iuplo)-1]
				//
				//              Call Dlattp to generate a triangular test matrix.
				//
				(*srnamt) = *func() *[]byte{y := []byte("Dlattp"); return &y }()
				Dlattp(imat, uplo, func() *[]byte{y := []byte("No transpose"); return &y }(), diag, iseed, n, (ap), (x), (work), info)
				//
				//              Set idiag = 1 for non-unit matrices, 2 for unit.
				//
				if blas.Lsame(diag, func() *byte{y := byte('N'); return &y }()) {
					(*idiag) = 1
				} else {
					(*idiag) = 2
				}
				//
				//+    TEST 1
				//              Form the inverse of A.
				//
				if (*n) > 0 {
					Dcopy(lap, (ap), func() *int{y := 1; return &y }(), (ainvp), func() *int{y := 1; return &y }())
				}
				(*srnamt) = *func() *[]byte{y := []byte("Dtptri"); return &y }()
				Dtptri(uplo, diag, n, (ainvp), info)
				//
				//              Check error code from Dtptri.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte{y := []byte("Dtptri"); return &y }(), info, func() *int{y := 0; return &y }(), append(append([]byte{}, (*uplo)), (*diag)), n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
				}
				//
				//              Compute the infinity-norm condition number of A.
				//
				(*anorm) = (*Dlantp(func() *byte{y := byte('I'); return &y }(), uplo, diag, n, (ap), (rwork)))
				(*ainvnm) = (*Dlantp(func() *byte{y := byte('I'); return &y }(), uplo, diag, n, (ainvp), (rwork)))
				if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
					(*rcondi) = (*one)
				} else {
					(*rcondi) = ((*one) / (*anorm)) / (*ainvnm)
				}
				//
				//              Compute the residual for the triangular matrix times its
				//              inverse.  Also compute the 1-norm condition number of A.
				//
				Dtpt01(uplo, diag, n, (ap), (ainvp), rcondo, (rwork), &((*result)[0]))
				//
				//              Print the test ratio if it is .GE. thresh.
				//
				if (*result)[0] >= (*(thresh)) {
					if (*nfail) == 0 && (*nerrs) == 0 {
						Alahd((nout), path)
					}
					WRITE((*(nout)), *func() *[]byte{y := []byte(" uplo='%c', diag='%c', N=%5d, type %2d, test(%2d)= %12.5f\n"); return &y }(), (*uplo), (*diag), (*n), (*imat), 1, (*result)[0])
					(*nfail) = (*nfail) + 1
				}
				(*nrun) = (*nrun) + 1
				//
				for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
					(*nrhs) = (*(nsval))[(*irhs)-1]
					(*xtype) = 'N'
					//
					for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
						//
						//                 Do for op(a) = a, A**T, or A**H.
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
						//+    TEST 2
						//                 Solve and compute residual for op(a)*x = b.
						//
						(*srnamt) = *func() *[]byte{y := []byte("Dlarhs"); return &y }()
						Dlarhs(path, xtype, uplo, trans, n, n, func() *int{y := 0; return &y }(), idiag, nrhs, (ap), lap, (xact), lda, (b), lda, iseed, info)
						(*xtype) = 'C'
						Dlacpy(func() *[]byte{y := []byte("Full"); return &y }(), n, nrhs, (b), lda, (x), lda)
						//
						(*srnamt) = *func() *[]byte{y := []byte("Dtptrs"); return &y }()
						Dtptrs(uplo, trans, diag, n, nrhs, (ap), (x), lda, info)
						//
						//                 Check error code from Dtptrs.
						//
						if (*info) != 0 {
							Alaerh ((*path), "Dtptrs", (*info), 0, (*uplo) append ( append ([] byte { }, append ( append ([] byte { },), (*trans))), (*diag)), (*N), (*N), - 1, - 1, - 1, (*imat), (*nfail), (*nerrs), (*nout))
						}
						//
						Dtpt02(uplo, trans, diag, n, nrhs, (ap), (x), lda, (b), lda, (work), &((*result)[1]))
						//
						//+    TEST 3
						//                 Check solution from generated exact solution.
						//
						Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[2]))
						//
						//+    TESts 4, 5, and 6
						//                 Use iterative refinement to improve the solution and
						//                 compute error bounds.
						//
						(*srnamt) = *func() *[]byte{y := []byte("Dtprfs"); return &y }()
						Dtprfs(uplo, trans, diag, n, nrhs, (ap), (b), lda, (x), lda, (rwork), &((*(rwork))[(*nrhs)+0]), (work), (iwork), info)
						//
						//                 Check error code from Dtprfs.
						//
						if (*info) != 0 {
							Alaerh ((*path), "Dtprfs", (*info), 0, (*uplo) append ( append ([] byte { }, append ( append ([] byte { },), (*trans))), (*diag)), (*N), (*N), - 1, - 1, (*nrhs), (*imat), (*nfail), (*nerrs), (*nout))
						}
						//
						Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[3]))
						Dtpt05(uplo, trans, diag, n, nrhs, (ap), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*nrhs)+0]), &((*result)[4]))
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						for (*k) = 2; (*k) <= 6; (*k)++ {
							if (*result)[(*k)-1] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Alahd((nout), path)
								}
								WRITE((*(nout)), *func() *[]byte{y := []byte(" uplo='%c', trans='%c', diag='%c', N=%5d', nrhs=%5d, type %2d, test(%2d)= %12.5f\n"); return &y }(), (*uplo), (*trans), (*diag), (*n), (*nrhs), (*imat), (*k), (*result)[(*k)-1])
								(*nfail) = (*nfail) + 1
							}
						//Label20:
						}
						(*nrun) = (*nrun) + 5
					//Label30:
					}
				//Label40:
				}
				//
				//+    TEST 7
				//                 Get an estimate of rcond = 1/cndnum.
				//
				for (*itran) = 1; (*itran) <= 2; (*itran)++ {
					if (*itran) == 1 {
						(*norm) = 'O'
						(*rcondc) = (*rcondo)
					} else {
						(*norm) = 'I'
						(*rcondc) = (*rcondi)
					}
					//
					(*srnamt) = *func() *[]byte{y := []byte("Dtpcon"); return &y }()
					Dtpcon(NORM, uplo, diag, n, (ap), rcond, (work), (iwork), info)
					//
					//                 Check error code from Dtpcon.
					//
					if (*info) != 0 {
						Alaerh ((*path), "Dtpcon", (*info), 0, (*NORM) append ( append ([] byte { }, append ( append ([] byte { },), (*uplo))), (*diag)), (*N), (*N), - 1, - 1, - 1, (*imat), (*nfail), (*nerrs), (*nout))
					}
					//
					Dtpt06(rcond, rcondc, uplo, diag, n, (ap), (rwork), &((*result)[6]))
					//
					//                 Print the test ratio if it is .GE. thresh.
					//
					if (*result)[6] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Alahd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c',%5d, ...), type %2d, test(%2d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dtpcon"); return &y }(), (*norm), (*uplo), (*diag), (*n), (*imat), 7, (*result)[6])
						(*nfail) = (*nfail) + 1
					}
					(*nrun) = (*nrun) + 1
				//Label50:
				}
			//Label60:
			}
		Label70:
		}
		//
		//        Use pathological test matrices to test Dlatps.
		//
		for (*imat) = (*ntype1) + 1; (*imat) <= (*ntypes); (*imat)++ {
			//
			//           Do the tests only if dotype( imat) is true.
			//
			if !(*(dotype))[(*imat)-1] {
				goto Label100
			}
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				//
				//              Do first for uplo = 'U', then for uplo = 'L'
				//
				(*uplo) = (*uplos)[(*iuplo)-1]
				for (*itran) = 1; (*itran) <= (*ntran); (*itran)++ {
					//
					//                 Do for op(a) = a, A**T, or A**H.
					//
					(*trans) = (*transs)[(*itran)-1]
					//
					//                 Call Dlattp to generate a triangular test matrix.
					//
					(*srnamt) = *func() *[]byte{y := []byte("Dlattp"); return &y }()
					Dlattp(imat, uplo, trans, diag, iseed, n, (ap), (x), (work), info)
					//
					//+    TEST 8
					//                 Solve the system op(a)*x = b.
					//
					(*srnamt) = *func() *[]byte{y := []byte("Dlatps"); return &y }()
					Dcopy(n, (x), func() *int{y := 1; return &y }(), (b), func() *int{y := 1; return &y }())
					Dlatps(uplo, trans, diag, func() *byte{y := byte('N'); return &y }(), n, (ap), (b), scale, (rwork), info)
					//
					//                 Check error code from Dlatps.
					//
					if (*info) != 0 {
						Alaerh ((*path), "Dlatps", (*info), 0, append ( append ([] byte { }, (*uplo)), (*trans)) append ( append ([] byte { }, append ( append ([] byte { },), (*diag))), "N"), (*N), (*N), - 1, - 1, - 1, (*imat), (*nfail), (*nerrs), (*nout))
					}
					//
					Dtpt03(uplo, trans, diag, n, func() *int{y := 1; return &y }(), (ap), scale, (rwork), one, (b), lda, (x), lda, (work), &((*result)[7]))
					//
					//+    TEST 9
					//                 Solve op(a)*x = b again with NOrmin = 'Y'.
					//
					Dcopy(n, (x), func() *int{y := 1; return &y }(), &((*(b))[(*n)+0]), func() *int{y := 1; return &y }())
					Dlatps(uplo, trans, diag, func() *byte{y := byte('Y'); return &y }(), n, (ap), &((*(b))[(*n)+0]), scale, (rwork), info)
					//
					//                 Check error code from Dlatps.
					//
					if (*info) != 0 {
						Alaerh ((*path), "Dlatps", (*info), 0, append ( append ([] byte { }, (*uplo)), (*trans)) append ( append ([] byte { }, append ( append ([] byte { },), (*diag))), "Y"), (*N), (*N), - 1, - 1, - 1, (*imat), (*nfail), (*nerrs), (*nout))
					}
					//
					Dtpt03(uplo, trans, diag, n, func() *int{y := 1; return &y }(), (ap), scale, (rwork), one, &((*(b))[(*n)+0]), lda, (x), lda, (work), &((*result)[8]))
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					if (*result)[7] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Alahd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c', '%c',%5d, ...), type %2d, test(%2d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dlatps"); return &y }(), (*uplo), (*trans), (*diag), 'N', (*n), (*imat), 8, (*result)[7])
						(*nfail) = (*nfail) + 1
					}
					if (*result)[8] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Alahd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte{y := []byte(" %s( '%c', '%c', '%c', '%c',%5d, ...), type %2d, test(%2d)=%12.5f\n"); return &y }(), *func() *[]byte{y := []byte("Dlatps"); return &y }(), (*uplo), (*trans), (*diag), 'Y', (*n), (*imat), 9, (*result)[8])
						(*nfail) = (*nfail) + 1
					}
					(*nrun) = (*nrun) + 2
				//Label80:
				}
			//Label90:
			}
		Label100:
		}
	//Label110:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchktp
	//
}
