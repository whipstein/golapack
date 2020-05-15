package goblas

import 

// Ddrvsp tests the driver routines Dspsv and -SVX.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvsp( dotype, nn, nval, nrhs, thresh, tsterr, nmax,
//                          a, afac, ainv, b, x, xact, work, rwork, iwork,
//                          nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       inTEGER            nmax, nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       inTEGER            iwork(*), nval(*)
//       DOUBLE PRECISION   a(*), afac(*), ainv(*), B(*),
//      $                   rwork(*), work(*), X(*), XACt(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvsp tests the driver routines Dspsv and -SVX.
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
//          nn is inTEGER
//          The number of values of N contained in the vector nval.
// \endverbatim
//
// \param[in] nval
// \verbatim
//          nval is inTEGER array, dimension (nn)
//          The values of the matrix dimension N.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
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
//          nmax is inTEGER
//          The maximum value permitted for n, used in dimensioning the
//          work arrays.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
// \endverbatim
//
// \param[out] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
// \endverbatim
//
// \param[out] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
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
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (nmax*max(2,nrhs))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (nmax+2*nrhs)
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is inTEGER array, dimension (2*nmax)
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
// \date December 2016
//
// \ingroup double_lin
//
//  =====================================================================
func Ddrvsp(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	nfact := new(int)
	zerot := new(bool)
	dist := new(byte)
	fact := new(byte)
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
	ifact := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	ioff := new(int)
	iuplo := new(int)
	izero := new(int)
	j := new(int)
	k := new(int)
	k1 := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	lwork := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	npp := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	facts := func() *[]byte {
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
		arr := make([]float64, 6)
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
	(*ntypes) = 10
	(*ntests) = 6
	(*nfact) = 2
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Scalars in Common ..
	//     ..
	//     .. Common blocks ..
	infot = common.infoc.infot
	nunit = common.infoc.nunit
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Data statements ..
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3] = 1988, 1989, 1990, 1991
	(*facts)[0], (*facts)[1] = 'F', 'N'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y}()
	(*path)[1] = *func() *[]byte {y :=[]byte("SP"); return &y}()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-1] = (*iseedy)[(*i)-1]
		//Label10:
	}
	(*lwork) = (MAX(2*(*(nmax)), (*(nmax))*(*(nrhs))))
	//
	//     Test the error exits
	//
	if *(tsterr) {
		Derrvx(path, (nout))
	}
	(*infot) = 0
	//
	//     Do for each value of N in nval
	//
	for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
		(*n) = (*(nval))[(*in)-1]
		(*lda) = (MAX((*n), 1))
		(*npp) = (*n) * ((*n) + 1) / 2
		(*xtype) = 'N'
		(*nimat) = (*ntypes)
		if (*n) <= 0 {
			(*nimat) = 1
		}
		//
		for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
			//
			//           Do the tests only if dotype( imat) is true.
			//
			if !(*(dotype))[(*imat)-1] {
				goto Label170
			}
			//
			//           Skip types 3, 4, 5, or 6 if the matrix size is too small.
			//
			(*zerot) = (*imat) >= 3 && (*imat) <= 6
			if (*zerot) && (*n) < (*imat)-2 {
				goto Label170
			}
			//
			//           Do first for uplo = 'U', then for uplo = 'L'
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				if (*iuplo) == 1 {
					(*uplo) = 'U'
					(*packit) = 'C'
				} else {
					(*uplo) = 'L'
					(*packit) = 'R'
				}
				//
				//              Set up parameters with Dlatb4 and generate a test matrix
				//              with Dlatms.
				//
				Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
				//
				(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y}()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, packit, (a), lda, (work), info)
				//
				//              Check error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y}(), info, func() *int {y := 0; return &y}(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					goto Label160
				}
				//
				//              For types 3-6, zero one or more rows and columns of the
				//              matrix to test that info is returned correctly.
				//
				if *zerot {
					if (*imat) == 3 {
						(*izero) = 1
					} else if (*imat) == 4 {
						(*izero) = (*n)
					} else {
						(*izero) = (*n)/2 + 1
					}
					//
					if (*imat) < 6 {
						//
						//                    Set row and column izero to zero.
						//
						if (*iuplo) == 1 {
							(*ioff) = ((*izero) - 1) * (*izero) / 2
							for (*i) = 1; (*i) <= (*izero)-1; (*i)++ {
								(*(a))[(*ioff)+(*i)-1] = (*zero)
								//Label20:
							}
							(*ioff) = (*ioff) + (*izero)
							for (*i) = (*izero); (*i) <= (*n); (*i)++ {
								(*(a))[(*ioff)-1] = (*zero)
								(*ioff) = (*ioff) + (*i)
								//Label30:
							}
						} else {
							(*ioff) = (*izero)
							for (*i) = 1; (*i) <= (*izero)-1; (*i)++ {
								(*(a))[(*ioff)-1] = (*zero)
								(*ioff) = (*ioff) + (*n) - (*i)
								//Label40:
							}
							(*ioff) = (*ioff) - (*izero)
							for (*i) = (*izero); (*i) <= (*n); (*i)++ {
								(*(a))[(*ioff)+(*i)-1] = (*zero)
								//Label50:
							}
						}
					} else {
						(*ioff) = 0
						if (*iuplo) == 1 {
							//
							//                       Set the first izero rows and columns to zero.
							//
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								(*i2) = (Min((*j), (*izero)))
								for (*i) = 1; (*i) <= (*i2); (*i)++ {
									(*(a))[(*ioff)+(*i)-1] = (*zero)
									//Label60:
								}
								(*ioff) = (*ioff) + (*j)
								//Label70:
							}
						} else {
							//
							//                       Set the last izero rows and columns to zero.
							//
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								(*i1) = (MAX((*j), (*izero)))
								for (*i) = (*i1); (*i) <= (*n); (*i)++ {
									(*(a))[(*ioff)+(*i)-1] = (*zero)
									//Label80:
								}
								(*ioff) = (*ioff) + (*n) - (*j)
								//Label90:
							}
						}
					}
				} else {
					(*izero) = 0
				}
				//
				for (*ifact) = 1; (*ifact) <= (*nfact); (*ifact)++ {
					//
					//                 Do first for fact = 'F', then for other values.
					//
					(*fact) = (*facts)[(*ifact)-1]
					//
					//                 Compute the condition number for comparison with
					//                 the value returned by Dspsvx.
					//
					if *zerot {
						if (*ifact) == 1 {
							goto Label150
						}
						(*rcondc) = (*zero)
						//
					} else if (*ifact) == 1 {
						//
						//                    Compute the 1-norm of A.
						//
						(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y}(), uplo, n, (a), (rwork)))
						//
						//                    Factor the matrix A.
						//
						Dcopy(npp, (a), func() *int {y := 1; return &y}(), (afac), func() *int {y := 1; return &y}())
						Dsptrf(uplo, n, (afac), (iwork), info)
						//
						//                    Compute inv(a) and take its norm.
						//
						Dcopy(npp, (afac), func() *int {y := 1; return &y}(), (ainv), func() *int {y := 1; return &y}())
						Dsptri(uplo, n, (ainv), (iwork), (work), info)
						(*ainvnm) = (*Dlansp(func() *byte {y := byte('1'); return &y}(), uplo, n, (ainv), (rwork)))
						//
						//                    Compute the 1-norm condition number of A.
						//
						if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
							(*rcondc) = (*one)
						} else {
							(*rcondc) = ((*one) / (*anorm)) / (*ainvnm)
						}
					}
					//
					//                 Form an exact solution and set the right hand side.
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y}()
					Dlarhs(path, xtype, uplo, func() *byte {y := byte(' '); return &y}(), n, n, kl, ku, (nrhs), (a), lda, (xact), lda, (b), lda, iseed, info)
					(*xtype) = 'C'
					//
					//                 --- Test Dspsv  ---
					//
					if (*ifact) == 2 {
						Dcopy(npp, (a), func() *int {y := 1; return &y}(), (afac), func() *int {y := 1; return &y}())
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (x), lda)
						//
						//                    Factor the matrix and solve the system using Dspsv.
						//
						(*srnamt) = *func() *[]byte {y :=[]byte("Dspsv "); return &y}()
						Dspsv(uplo, n, (nrhs), (afac), (iwork), (x), lda, info)
						//
						//                    Adjust the expected value of info to account for
						//                    pivoting.
						//
						(*k) = (*izero)
						if (*k) > 0 {
						Label100:
							;
							if (*(iwork))[(*k)-1] < 0 {
								if (*(iwork))[(*k)-1] != -(*k) {
									(*k) = -(*(iwork))[(*k)-1]
									goto Label100
								}
							} else if (*(iwork))[(*k)-1] != (*k) {
								(*k) = (*(iwork))[(*k)-1]
								goto Label100
							}
						}
						//
						//                    Check error code from Dspsv .
						//
						if (*info) != (*k) {
							Alaerh(path, func() *[]byte {y :=[]byte("Dspsv "); return &y}(), info, k, uplo, n, n, -1, -1, (nrhs), imat, nfail, nerrs, (nout))
							goto Label120
						} else if (*info) != 0 {
							goto Label120
						}
						//
						//                    Re_construct matrix from factors and compute
						//                    residual.
						//
						Dspt01(uplo, n, (a), (afac), (iwork), (ainv), lda, (rwork), &((*result)[0]))
						//
						//                    Compute residual of the computed solution.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
						Dppt02(uplo, n, (nrhs), (a), (x), lda, (work), lda, (rwork), &((*result)[1]))
						//
						//                    Check solution from generated exact solution.
						//
						Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
						(*nt) = 3
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						for (*k) = 1; (*k) <= (*nt); (*k)++ {
							if (*result)[(*k)-1] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Aladhd((nout), path)
								}
								WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %s, uplo='%c', N =%5d, type %2d, test %2d, ratio =%12.5f\n"); return &y}(), *func() *[]byte {y :=[]byte("Dspsv "); return &y}(), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
								(*nfail) = (*nfail) + 1
							}
							//Label110:
						}
						(*nrun) = (*nrun) + (*nt)
					Label120:
					}
					//
					//                 --- Test Dspsvx ---
					//
					if (*ifact) == 2 && (*npp) > 0 {
						Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), npp, func() *int {y := 1; return &y}(), zero, zero, (afac), npp)
					}
					Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), zero, zero, (x), lda)
					//
					//                 Solve the system and compute the condition number and
					//                 error bounds using Dspsvx.
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dspsvx"); return &y}()
					Dspsvx(fact, uplo, n, (nrhs), (a), (afac), (iwork), (b), lda, (x), lda, rcond, (rwork), &((*(rwork))[(*(nrhs))+0]), (work), &((*(iwork))[(*n)+0]), info)
					//
					//                 Adjust the expected value of info to account for
					//                 pivoting.
					//
					(*k) = (*izero)
					if (*k) > 0 {
					Label130:
						;
						if (*(iwork))[(*k)-1] < 0 {
							if (*(iwork))[(*k)-1] != -(*k) {
								(*k) = -(*(iwork))[(*k)-1]
								goto Label130
							}
						} else if (*(iwork))[(*k)-1] != (*k) {
							(*k) = (*(iwork))[(*k)-1]
							goto Label130
						}
					}
					//
					//                 Check the error code from Dspsvx.
					//
					if (*info) != (*k) {
						Alaerh(path, func() *[]byte {y :=[]byte("Dspsvx"); return &y}(), info, k, append(append([]byte{}, (*fact)), (*uplo)), n, n, -1, -1, (nrhs), imat, nfail, nerrs, (nout))
						goto Label150
					}
					//
					if (*info) == 0 {
						if (*ifact) >= 2 {
							//
							//                       Re_construct matrix from factors and compute
							//                       residual.
							//
							Dspt01(uplo, n, (a), (afac), (iwork), (ainv), lda, &((*(rwork))[2*(*(nrhs))+0]), &((*result)[0]))
							(*k1) = 1
						} else {
							(*k1) = 2
						}
						//
						//                    Compute residual of the computed solution.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
						Dppt02(uplo, n, (nrhs), (a), (x), lda, (work), lda, &((*(rwork))[2*(*(nrhs))+0]), &((*result)[1]))
						//
						//                    Check solution from generated exact solution.
						//
						Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
						//
						//                    Check the error bounds from iterative refinement.
						//
						Dppt05(uplo, n, (nrhs), (a), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*(nrhs))+0]), &((*result)[3]))
					} else {
						(*k1) = 6
					}
					//
					//                 Compare rcond from Dspsvx with the computed value
					//                 in rcondc.
					//
					(*result)[5] = (*Dget06(rcond, rcondc))
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = (*k1); (*k) <= 6; (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Aladhd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" %s, fact='%c', uplo='%c', N =%5d, type %2d, test %2d, ratio =%12.5f\n")
								return &y
							}(), *func() *[]byte {y :=[]byte("Dspsvx"); return &y}(), (*fact), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label140:
					}
					(*nrun) = (*nrun) + 7 - (*k1)
					//
				Label150:
				}
				//
			Label160:
			}
		Label170:
		}
		//Label180:
	}
	//
	//     Print a summary of the results.
	//
	Alasvm(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Ddrvsp
	//
}
