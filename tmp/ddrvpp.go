package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Ddrvpp tests the driver routines Dppsv and -SVX.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvpp( dotype, nn, nval, nrhs, thresh, tsterr, nmax,
//                          a, afac, asav, b, bsav, x, xact, S, work,
//                          rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       inTEGER            nmax, nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       inTEGER            iwork(*), nval(*)
//       DOUBLE PRECISION   a(*), afac(*), asav(*), B(*),
//      $                   bsav(*), rwork(*), S(*), work(*),
//      $                   X(*), XACt(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvpp tests the driver routines Dppsv and -SVX.
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
// \param[out] asav
// \verbatim
//          asav is DOUBLE PRECISION array, dimension
//                      (nmax*(nmax+1)/2)
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (nmax*nrhs)
// \endverbatim
//
// \param[out] bsav
// \verbatim
//          bsav is DOUBLE PRECISION array, dimension (nmax*nrhs)
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
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension (nmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (nmax*max(3,nrhs))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (nmax+2*nrhs)
// \endverbatim
//
// \param[out] iwork
// \verbatim
//          iwork is inTEGER array, dimension (nmax)
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
func Ddrvpp(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, asav *[]float64, b *[]float64, bsav *[]float64, x *[]float64, xact *[]float64, s *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	equil := new(bool)
	nofact := new(bool)
	prefac := new(bool)
	zerot := new(bool)
	dist := new(byte)
	equed := new(byte)
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
	iequed := new(int)
	ifact := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	ioff := new(int)
	iuplo := new(int)
	izero := new(int)
	k := new(int)
	k1 := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfact := new(int)
	nfail := new(int)
	nimat := new(int)
	npp := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	amax := new(float64)
	anorm := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	amax := new(float64)
	scond := new(float64)
	equeds := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	facts := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	packs := func() *[]byte {
		arr := make([]byte, 2)
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
	(*ntypes) = 9
	(*ntests) = 6
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
	(*uplos)[0], (*uplos)[1], (*facts)[0], (*facts)[1], (*facts)[2], (*packs)[0], (*packs)[1], (*equeds)[0], (*equeds)[1] = 'U', 'L', 'F', 'N', 'E', 'C', 'R', 'N', 'Y'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y := []byte("PP"); return &y }()
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
				goto Label130
			}
			//
			//           Skip types 3, 4, or 5 if the matrix size is too small.
			//
			(*zerot) = (*imat) >= 3 && (*imat) <= 5
			if (*zerot) && (*n) < (*imat)-2 {
				goto Label130
			}
			//
			//           Do first for uplo = 'U', then for uplo = 'L'
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				(*uplo) = (*uplos)[(*iuplo)-1]
				(*packit) = (*packs)[(*iuplo)-1]
				//
				//              Set up parameters with Dlatb4 and generate a test matrix
				//              with Dlatms.
				//
				Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
				(*rcondc) = (*one) / (*cndnum)
				//
				(*srnamt) = *func() *[]byte {y := []byte("Dlatms"); return &y }()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, packit, (a), lda, (work), info)
				//
				//              Check error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y := []byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					goto Label120
				}
				//
				//              For types 3-5, zero one row and column of the matrix to
				//              test that info is returned correctly.
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
					//                 Set row and column izero of A to 0.
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
					(*izero) = 0
				}
				//
				//              Save a copy of the matrix A in asav.
				//
				Dcopy(npp, (a), func() *int {y := 1; return &y }(), (asav), func() *int {y := 1; return &y }())
				//
				for (*iequed) = 1; (*iequed) <= 2; (*iequed)++ {
					(*equed) = (*equeds)[(*iequed)-1]
					if (*iequed) == 1 {
						(*nfact) = 3
					} else {
						(*nfact) = 1
					}
					//
					for (*ifact) = 1; (*ifact) <= (*nfact); (*ifact)++ {
						(*fact) = (*facts)[(*ifact)-1]
						(*prefac) = (*blas.Lsame(fact, func() *byte {y := byte('F'); return &y }()))
						(*nofact) = (*blas.Lsame(fact, func() *byte {y := byte('N'); return &y }()))
						(*equil) = (*blas.Lsame(fact, func() *byte {y := byte('E'); return &y }()))
						//
						if *zerot {
							if *prefac {
								goto Label100
							}
							(*rcondc) = (*zero)
							//
						} else if !blas.Lsame(fact, func() *byte {y := byte('N'); return &y }()) {
							//
							//                       Compute the condition number for comparison with
							//                       the value returned by Dppsvx (fact = 'N' reuses
							//                       the condition number from the previous iteration
							//                       with fact = 'F').
							//
							Dcopy(npp, (asav), func() *int {y := 1; return &y }(), (afac), func() *int {y := 1; return &y }())
							if (*equil) || (*iequed) > 1 {
								//
								//                          Compute row and column scale factors to
								//                          equilibrate the matrix A.
								//
								Dppequ(uplo, n, (afac), (s), scond, amax, info)
								if (*info) == 0 && (*n) > 0 {
									if (*iequed) > 1 {
										(*scond) = (*zero)
									}
									//
									//                             Equilibrate the matrix.
									//
									Dlaqsp(uplo, n, (afac), (s), scond, amax, equed)
								}
							}
							//
							//                       Save the condition number of the
							//                       non-equilibrated system for use in Dget04.
							//
							if *equil {
								(*amax) = (*rcondc)
							}
							//
							//                       Compute the 1-norm of A.
							//
							(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), uplo, n, (afac), (rwork)))
							//
							//                       Factor the matrix A.
							//
							Dpptrf(uplo, n, (afac), info)
							//
							//                       Form the inverse of A.
							//
							Dcopy(npp, (afac), func() *int {y := 1; return &y }(), (a), func() *int {y := 1; return &y }())
							Dpptri(uplo, n, (a), info)
							//
							//                       Compute the 1-norm condition number of A.
							//
							(*ainvnm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), uplo, n, (a), (rwork)))
							if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
								(*rcondc) = (*one)
							} else {
								(*rcondc) = ((*one) / (*anorm)) / (*ainvnm)
							}
						}
						//
						//                    Restore the matrix A.
						//
						Dcopy(npp, (asav), func() *int {y := 1; return &y }(), (a), func() *int {y := 1; return &y }())
						//
						//                    Form an exact solution and set the right hand side.
						//
						(*srnamt) = *func() *[]byte {y := []byte("Dlarhs"); return &y }()
						Dlarhs(path, xtype, uplo, func() *byte {y := byte(' '); return &y }(), n, n, kl, ku, (nrhs), (a), lda, (xact), lda, (b), lda, iseed, info)
						(*xtype) = 'C'
						Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (bsav), lda)
						//
						if *nofact {
							//
							//                       --- Test Dppsv  ---
							//
							//                       Compute the L*L' or U'*U factorization of the
							//                       matrix and solve the system.
							//
							Dcopy(npp, (a), func() *int {y := 1; return &y }(), (afac), func() *int {y := 1; return &y }())
							Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (x), lda)
							//
							(*srnamt) = *func() *[]byte {y := []byte("Dppsv "); return &y }()
							Dppsv(uplo, n, (nrhs), (afac), (x), lda, info)
							//
							//                       Check error code from Dppsv .
							//
							if (*info) != (*izero) {
								Alaerh(path, func() *[]byte {y := []byte("Dppsv "); return &y }(), info, izero, uplo, n, n, -1, -1, (nrhs), imat, nfail, nerrs, (nout))
								goto Label70
							} else if (*info) != 0 {
								goto Label70
							}
							//
							//                       Re_construct matrix from factors and compute
							//                       residual.
							//
							Dppt01(uplo, n, (a), (afac), (rwork), &((*result)[0]))
							//
							//                       Compute residual of the computed solution.
							//
							Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (work), lda)
							Dppt02(uplo, n, (nrhs), (a), (x), lda, (work), lda, (rwork), &((*result)[1]))
							//
							//                       Check solution from generated exact solution.
							//
							Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
							(*nt) = 3
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = 1; (*k) <= (*nt); (*k)++ {
								if (*result)[(*k)-1] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Aladhd((nout), path)
									}
									WRITE((*(nout)), *func() *[]byte {y := []byte(" %s, uplo='%c', N =%5d, type %1d, test(%1d)=%12.5f\n"); return &y }(), *func() *[]byte {y := []byte("Dppsv "); return &y }(), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
									(*nfail) = (*nfail) + 1
								}
								//Label60:
							}
							(*nrun) = (*nrun) + (*nt)
						Label70:
						}
						//
						//                    --- Test Dppsvx ---
						//
						if !(*prefac) && (*npp) > 0 {
							Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), npp, func() *int {y := 1; return &y }(), zero, zero, (afac), npp)
						}
						Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), zero, zero, (x), lda)
						if (*iequed) > 1 && (*n) > 0 {
							//
							//                       Equilibrate the matrix if fact='F' and
							//                       equed='Y'.
							//
							Dlaqsp(uplo, n, (a), (s), scond, amax, equed)
						}
						//
						//                    Solve the system and compute the condition number
						//                    and error bounds using Dppsvx.
						//
						(*srnamt) = *func() *[]byte {y := []byte("Dppsvx"); return &y }()
						Dppsvx(fact, uplo, n, (nrhs), (a), (afac), equed, (s), (b), lda, (x), lda, rcond, (rwork), &((*(rwork))[(*(nrhs))+0]), (work), (iwork), info)
						//
						//                    Check the error code from Dppsvx.
						//
						if (*info) != (*izero) {
							Alaerh(path, func() *[]byte {y := []byte("Dppsvx"); return &y }(), info, izero, append(append([]byte{}, (*fact)), (*uplo)), n, n, -1, -1, (nrhs), imat, nfail, nerrs, (nout))
							goto Label90
						}
						//
						if (*info) == 0 {
							if !(*prefac) {
								//
								//                          Re_construct matrix from factors and compute
								//                          residual.
								//
								Dppt01(uplo, n, (a), (afac), &((*(rwork))[2*(*(nrhs))+0]), &((*result)[0]))
								(*k1) = 1
							} else {
								(*k1) = 2
							}
							//
							//                       Compute residual of the computed solution.
							//
							Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (bsav), lda, (work), lda)
							Dppt02(uplo, n, (nrhs), (asav), (x), lda, (work), lda, &((*(rwork))[2*(*(nrhs))+0]), &((*result)[1]))
							//
							//                       Check solution from generated exact solution.
							//
							if (*nofact) || ((*prefac) && blas.Lsame(equed, func() *byte {y := byte('N'); return &y }())) {
								Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
							} else {
								Dget04(n, (nrhs), (x), lda, (xact), lda, amax, &((*result)[2]))
							}
							//
							//                       Check the error bounds from iterative
							//                       refinement.
							//
							Dppt05(uplo, n, (nrhs), (asav), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*(nrhs))+0]), &((*result)[3]))
						} else {
							(*k1) = 6
						}
						//
						//                    Compare rcond from Dppsvx with the computed value
						//                    in rcondc.
						//
						(*result)[5] = (*Dget06(rcond, rcondc))
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						for (*k) = (*k1); (*k) <= 6; (*k)++ {
							if (*result)[(*k)-1] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Aladhd((nout), path)
								}
								if *prefac {
									WRITE((*(nout)), *func() *[]byte {
										y := []byte(" %s, fact='%c', uplo='%c', N=%5d, equed='%c', type %1d, test(%1d)=%12.5f\n")
										return &y
									}(), *func() *[]byte {y := []byte("Dppsvx"); return &y }(), (*fact), (*uplo), (*n), (*equed), (*imat), (*k), (*result)[(*k)-1])
								} else {
									WRITE((*(nout)), *func() *[]byte {
										y := []byte(" %s, fact='%c', uplo='%c', N=%5d, type %1d, test(%1d)=%12.5f\n")
										return &y
									}(), *func() *[]byte {y := []byte("Dppsvx"); return &y }(), (*fact), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
								}
								(*nfail) = (*nfail) + 1
							}
							//Label80:
						}
						(*nrun) = (*nrun) + 7 - (*k1)
					Label90:
						;
					Label100:
					}
					//Label110:
				}
			Label120:
			}
		Label130:
		}
		//Label140:
	}
	//
	//     Print a summary of the results.
	//
	Alasvm(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Ddrvpp
	//
}
