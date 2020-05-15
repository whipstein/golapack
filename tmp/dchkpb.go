package goblas

import 

// Dchkpb tests DPBTRF, -TRS, -RFS, and -CON.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkpb( dotype, nn, nval, nnb, nbval, nns, nsval,
//                          thresh, tsterr, nmax, a, afac, ainv, b, x,
//                          xact, work, rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nmax, nn, nnb, nns, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), nbval(*), nsval(*), nval(*)
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
// Dchkpb tests DPBTRF, -TRS, -RFS, and -CON.
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
//          The values of the matrix dimension N.
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
//          The maximum value permitted for n, used in dimensioning the
//          work arrays.
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
func Dchkpb(dotype *[]bool, nn *int, nval *[]int, nnb *int, nbval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	nbw := new(int)
	zerot := new(bool)
	dist := new(byte)
	packit := new(byte)
	_type := new(byte)
	uplo := new(byte)
	xtype := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	i1 := new(int)
	i2 := new(int)
	ikd := new(int)
	imat := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	ioff := new(int)
	irhs := new(int)
	iuplo := new(int)
	iw := new(int)
	izero := new(int)
	k := new(int)
	kd := new(int)
	kl := new(int)
	koff := new(int)
	ku := new(int)
	lda := new(int)
	ldab := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nkd := new(int)
	nrhs := new(int)
	nrun := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	iseedy := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	kdval := func() *[]int {
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
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
	(*ntypes) = 8
	(*ntests) = 7
	(*nbw) = 4
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
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3] = 1988, 1989, 1990, 1991
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y :=[]byte("PB"); return &y }()
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
		Derrpo(path, (nout))
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	(*kdval)[0] = 0
	//
	//     Do for each value of N in nval
	//
	for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
		(*n) = (*(nval))[(*in)-1]
		(*lda) = (MAX((*n), 1))
		(*xtype) = 'N'
		//
		//        Set limits on the number of loop iterations.
		//
		(*nkd) = (MAX(1, Min((*n), int(4))))
		(*nimat) = (*ntypes)
		if (*n) == 0 {
			(*nimat) = 1
		}
		//
		(*kdval)[1] = (*n) + ((*n)+1)/4
		(*kdval)[2] = (3*(*n) - 1) / 4
		(*kdval)[3] = ((*n) + 1) / 4
		//
		for (*ikd) = 1; (*ikd) <= (*nkd); (*ikd)++ {
			//
			//           Do for kd = 0, (5*n+1)/4, (3N-1)/4, and (N+1)/4. This order
			//           makes it easier to skip redundant values for small values
			//           of N.
			//
			(*kd) = (*kdval)[(*ikd)-1]
			(*ldab) = (*kd) + 1
			//
			//           Do first for uplo = 'U', then for uplo = 'L'
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				(*koff) = 1
				if (*iuplo) == 1 {
					(*uplo) = 'U'
					(*koff) = (MAX(1, (*kd)+2-(*n)))
					(*packit) = 'Q'
				} else {
					(*uplo) = 'L'
					(*packit) = 'B'
				}
				//
				for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
					//
					//                 Do the tests only if dotype( imat) is true.
					//
					if !(*(dotype))[(*imat)-1] {
						goto Label60
					}
					//
					//                 Skip types 2, 3, or 4 if the matrix size is too small.
					//
					(*zerot) = (*imat) >= 2 && (*imat) <= 4
					if (*zerot) && (*n) < (*imat)-1 {
						goto Label60
					}
					//
					if !(*zerot) || !(*(dotype))[0] {
						//
						//                    Set up parameters with Dlatb4 and generate a test
						//                    matrix with Dlatms.
						//
						Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
						//
						(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y }()
						Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kd, kd, packit, &((*(a))[(*koff)-1]), ldab, (work), info)
						//
						//                    Check error code from Dlatms.
						//
						if (*info) != 0 {
							Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, kd, kd, -1, imat, nfail, nerrs, (nout))
							goto Label60
						}
					} else if (*izero) > 0 {
						//
						//                    Use the same matrix for types 3 and 4 as for type
						//                    2 by copying back the zeroed out column,
						//
						(*iw) = 2*(*lda) + 1
						if (*iuplo) == 1 {
							(*ioff) = ((*izero)-1)*(*ldab) + (*kd) + 1
							Dcopy((*izero)-(*i1), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*ioff)-(*izero)+(*i1)-1]), func() *int {y := 1; return &y }())
							(*iw) = (*iw) + (*izero) - (*i1)
							Dcopy((*i2)-(*izero)+1, &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*ioff)-1]), MAX((*ldab)-1, 1))
						} else {
							(*ioff) = ((*i1)-1)*(*ldab) + 1
							Dcopy((*izero)-(*i1), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*ioff)+(*izero)-(*i1)-1]), MAX((*ldab)-1, 1))
							(*ioff) = ((*izero)-1)*(*ldab) + 1
							(*iw) = (*iw) + (*izero) - (*i1)
							Dcopy((*i2)-(*izero)+1, &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }(), &((*(a))[(*ioff)-1]), func() *int {y := 1; return &y }())
						}
					}
					//
					//                 For types 2-4, zero one row and column of the matrix
					//                 to test that info is returned correctly.
					//
					(*izero) = 0
					if *zerot {
						if (*imat) == 2 {
							(*izero) = 1
						} else if (*imat) == 3 {
							(*izero) = (*n)
						} else {
							(*izero) = (*n)/2 + 1
						}
						//
						//                    Save the zeroed out row and column in work(*,3)
						//
						(*iw) = 2 * (*lda)
						for (*i) = 1; (*i) <= (Min(2*(*kd)+1, (*n))); (*i)++ {
							(*(work))[(*iw)+(*i)-1] = (*zero)
							//Label20:
						}
						(*iw) = (*iw) + 1
						(*i1) = (MAX((*izero)-(*kd), 1))
						(*i2) = (Min((*izero)+(*kd), (*n)))
						//
						if (*iuplo) == 1 {
							(*ioff) = ((*izero)-1)*(*ldab) + (*kd) + 1
							Dswap((*izero)-(*i1), &((*(a))[(*ioff)-(*izero)+(*i1)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }())
							(*iw) = (*iw) + (*izero) - (*i1)
							Dswap((*i2)-(*izero)+1, &((*(a))[(*ioff)-1]), MAX((*ldab)-1, 1), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }())
						} else {
							(*ioff) = ((*i1)-1)*(*ldab) + 1
							Dswap((*izero)-(*i1), &((*(a))[(*ioff)+(*izero)-(*i1)-1]), MAX((*ldab)-1, 1), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }())
							(*ioff) = ((*izero)-1)*(*ldab) + 1
							(*iw) = (*iw) + (*izero) - (*i1)
							Dswap((*i2)-(*izero)+1, &((*(a))[(*ioff)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*iw)-1]), func() *int {y := 1; return &y }())
						}
					}
					//
					//                 Do for each value of nb in nbval
					//
					for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
						(*nb) = (*(nbval))[(*inb)-1]
						Xlaenv(func() *int {y := 1; return &y }(), nb)
						//
						//                    Compute the L*l' or U'*U factorization of the band
						//                    matrix.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (*kd)+1, n, (a), ldab, (afac), ldab)
						(*srnamt) = *func() *[]byte {y :=[]byte("DPBTRF"); return &y }()
						Dpbtrf(uplo, n, kd, (afac), ldab, info)
						//
						//                    Check error code from DPBTRF.
						//
						if (*info) != (*izero) {
							Alaerh(path, func() *[]byte {y :=[]byte("DPBTRF"); return &y }(), info, izero, uplo, n, n, kd, kd, nb, imat, nfail, nerrs, (nout))
							goto Label50
						}
						//
						//                    Skip the tests if info is not 0.
						//
						if (*info) != 0 {
							goto Label50
						}
						//
						//+    TEST 1
						//                    Re_construct matrix from factors and compute
						//                    residual.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), (*kd)+1, n, (afac), ldab, (ainv), ldab)
						Dpbt01(uplo, n, kd, (a), ldab, (ainv), ldab, (rwork), &((*result)[0]))
						//
						//                    Print the test ratio if it is .GE. thresh.
						//
						if (*result)[0] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" uplo='%c', N=%5d, kd=%5d, nb=%4d, type %2d, test %2d, ratio= %12.5f\n")
								return &y
							}(), (*uplo), (*n), (*kd), (*nb), (*imat), 1, (*result)[0])
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 1
						//
						//                    Only do other tests if this is the first blocksize.
						//
						if (*inb) > 1 {
							goto Label50
						}
						//
						//                    Form the inverse of A so we can get a good estimate
						//                    of rcondc = 1/(norm(a) * norm(inv(a))).
						//
						Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), n, n, zero, one, (ainv), lda)
						(*srnamt) = *func() *[]byte {y :=[]byte("DPBTRS"); return &y }()
						Dpbtrs(uplo, n, kd, n, (afac), ldab, (ainv), lda, info)
						//
						//                    Compute rcondc = 1/(norm(a) * norm(inv(a))).
						//
						(*anorm) = (*Dlansb(func() *byte {y := byte('1'); return &y }(), uplo, n, kd, (a), ldab, (rwork)))
						(*ainvnm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), n, n, (ainv), lda, (rwork)))
						if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
							(*rcondc) = (*one)
						} else {
							(*rcondc) = ((*one) / (*anorm)) / (*ainvnm)
						}
						//
						for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
							(*nrhs) = (*(nsval))[(*irhs)-1]
							//
							//+    TEST 2
							//                    Solve and compute residual for A * X = B.
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y }()
							Dlarhs(path, xtype, uplo, func() *byte {y := byte(' '); return &y }(), n, n, kd, kd, nrhs, (a), ldab, (xact), lda, (b), lda, iseed, info)
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), lda, (x), lda)
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("DPBTRS"); return &y }()
							Dpbtrs(uplo, n, kd, nrhs, (afac), ldab, (x), lda, info)
							//
							//                    Check error code from DPBTRS.
							//
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("DPBTRS"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, (nout))
							}
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, nrhs, (b), lda, (work), lda)
							Dpbt02(uplo, n, kd, nrhs, (a), ldab, (x), lda, (work), lda, (rwork), &((*result)[1]))
							//
							//+    TEST 3
							//                    Check solution from generated exact solution.
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[2]))
							//
							//+    TESts 4, 5, and 6
							//                    Use iterative refinement to improve the solution.
							//
							(*srnamt) = *func() *[]byte {y :=[]byte("DPBRFS"); return &y }()
							Dpbrfs(uplo, n, kd, nrhs, (a), ldab, (afac), ldab, (b), lda, (x), lda, (rwork), &((*(rwork))[(*nrhs)+0]), (work), (iwork), info)
							//
							//                    Check error code from DPBRFS.
							//
							if (*info) != 0 {
								Alaerh(path, func() *[]byte {y :=[]byte("DPBRFS"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, (nout))
							}
							//
							Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[3]))
							Dpbt05(uplo, n, kd, nrhs, (a), ldab, (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*nrhs)+0]), &((*result)[4]))
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = 2; (*k) <= 6; (*k)++ {
								if (*result)[(*k)-1] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Alahd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte {
										y :=[]byte(" uplo='%c', N=%5d, kd=%5d, nrhs=%3d, type %2d, test(%2d) = %12.5f\n")
										return &y
									}(), (*uplo), (*n), (*kd), (*nrhs), (*imat), (*k), (*result)[(*k)-1])
									(*nfail) = (*nfail) + 1
								}
								//Label30:
							}
							(*nrun) = (*nrun) + 5
							//Label40:
						}
						//
						//+    TEST 7
						//                    Get an estimate of rcond = 1/cndnum.
						//
						(*srnamt) = *func() *[]byte {y :=[]byte("DPBCON"); return &y }()
						Dpbcon(uplo, n, kd, (afac), ldab, anorm, rcond, (work), (iwork), info)
						//
						//                    Check error code from DPBCON.
						//
						if (*info) != 0 {
							Alaerh(path, func() *[]byte {y :=[]byte("DPBCON"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, kd, kd, -1, imat, nfail, nerrs, (nout))
						}
						//
						(*result)[6] = (*Dget06(rcond, rcondc))
						//
						//                    Print the test ratio if it is .GE. thresh.
						//
						if (*result)[6] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" uplo='%c', N=%5d, kd=%5d,           type %2d, test(%2d) = %12.5f\n")
								return &y
							}(), (*uplo), (*n), (*kd), (*imat), 7, (*result)[6])
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 1
					Label50:
					}
				Label60:
				}
				//Label70:
			}
			//Label80:
		}
		//Label90:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchkpb
	//
}
