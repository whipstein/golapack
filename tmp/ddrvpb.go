package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Ddrvpb tests the driver routines DPBSV and -SVX.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvpb( dotype, nn, nval, nrhs, thresh, tsterr, nmax,
//                          a, afac, asav, b, bsav, x, xact, s, work,
//                          rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nmax, nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), nval(*)
//       DOUBLE PRECISION   a(*), afac(*), asav(*), B(*),
//      $                   bsav(*), rwork(*), S(*), work(*),
//      $                   X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvpb tests the driver routines DPBSV and -SVX.
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
// \param[out] asav
// \verbatim
//          asav is DOUBLE PRECISION array, dimension (nmax*nmax)
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
func Ddrvpb(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, asav *[]float64, b *[]float64, bsav *[]float64, x *[]float64, xact *[]float64, s *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	nbw := new(int)
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
	i1 := new(int)
	i2 := new(int)
	iequed := new(int)
	ifact := new(int)
	ikd := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	ioff := new(int)
	iuplo := new(int)
	iw := new(int)
	izero := new(int)
	k := new(int)
	k1 := new(int)
	kd := new(int)
	kl := new(int)
	koff := new(int)
	ku := new(int)
	lda := new(int)
	ldab := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nbmin := new(int)
	nerrs := new(int)
	nfact := new(int)
	nfail := new(int)
	nimat := new(int)
	nkd := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	amax := new(float64)
	anorm := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
	roldc := new(float64)
	scond := new(float64)
	equeds := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	facts := func() *[]byte {
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
	kdval := func() *[]int {
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
	(*ntypes) = 8
	(*ntests) = 6
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
	(*facts)[0], (*facts)[1], (*facts)[2] = 'F', 'N', 'E'
	(*equeds)[0], (*equeds)[1] = 'N', 'Y'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y := []byte("PB"); return &y }()
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
	(*kdval)[0] = 0
	//
	//     Set the block size and minimum block size for testing.
	//
	(*nb) = 1
	(*nbmin) = 2
	Xlaenv(func() *int {y := 1; return &y }(), nb)
	Xlaenv(func() *int {y := 2; return &y }(), nbmin)
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
					(*packit) = 'Q'
					(*koff) = (MAX(1, (*kd)+2-(*n)))
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
						goto Label80
					}
					//
					//                 Skip types 2, 3, or 4 if the matrix size is too small.
					//
					(*zerot) = (*imat) >= 2 && (*imat) <= 4
					if (*zerot) && (*n) < (*imat)-1 {
						goto Label80
					}
					//
					if !(*zerot) || !(*(dotype))[0] {
						//
						//                    Set up parameters with Dlatb4 and generate a test
						//                    matrix with Dlatms.
						//
						Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
						//
						(*srnamt) = *func() *[]byte {y := []byte("Dlatms"); return &y }()
						Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kd, kd, packit, &((*(a))[(*koff)-1]), ldab, (work), info)
						//
						//                    Check error code from Dlatms.
						//
						if (*info) != 0 {
							Alaerh(path, func() *[]byte {y := []byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
							goto Label80
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
					//                 Save a copy of the matrix A in asav.
					//
					Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (*kd)+1, n, (a), ldab, (asav), ldab)
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
									goto Label60
								}
								(*rcondc) = (*zero)
								//
							} else if !blas.Lsame(fact, func() *byte {y := byte('N'); return &y }()) {
								//
								//                          Compute the condition number for comparison
								//                          with the value returned by Dpbsvx (fact =
								//                          'N' reuses the condition number from the
								//                          previous iteration with fact = 'F').
								//
								Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (*kd)+1, n, (asav), ldab, (afac), ldab)
								if (*equil) || (*iequed) > 1 {
									//
									//                             Compute row and column scale factors to
									//                             equilibrate the matrix A.
									//
									Dpbequ(uplo, n, kd, (afac), ldab, (s), scond, amax, info)
									if (*info) == 0 && (*n) > 0 {
										if (*iequed) > 1 {
											(*scond) = (*zero)
										}
										//
										//                                Equilibrate the matrix.
										//
										Dlaqsb(uplo, n, kd, (afac), ldab, (s), scond, amax, equed)
									}
								}
								//
								//                          Save the condition number of the
								//                          non-equilibrated system for use in Dget04.
								//
								if *equil {
									(*amax) = (*rcondc)
								}
								//
								//                          Compute the 1-norm of A.
								//
								(*anorm) = (*Dlansb(func() *byte {y := byte('1'); return &y }(), uplo, n, kd, (afac), ldab, (rwork)))
								//
								//                          Factor the matrix A.
								//
								Dpbtrf(uplo, n, kd, (afac), ldab, info)
								//
								//                          Form the inverse of A.
								//
								Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), n, n, zero, one, (a), lda)
								(*srnamt) = *func() *[]byte {y := []byte("DPBTRS"); return &y }()
								Dpbtrs(uplo, n, kd, n, (afac), ldab, (a), lda, info)
								//
								//                          Compute the 1-norm condition number of A.
								//
								(*ainvnm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), n, n, (a), lda, (rwork)))
								if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
									(*rcondc) = (*one)
								} else {
									(*rcondc) = ((*one) / (*anorm)) / (*ainvnm)
								}
							}
							//
							//                       Restore the matrix A.
							//
							Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (*kd)+1, n, (asav), ldab, (a), ldab)
							//
							//                       Form an exact solution and set the right hand
							//                       side.
							//
							(*srnamt) = *func() *[]byte {y := []byte("Dlarhs"); return &y }()
							Dlarhs(path, xtype, uplo, func() *byte {y := byte(' '); return &y }(), n, n, kd, kd, (nrhs), (a), ldab, (xact), lda, (b), lda, iseed, info)
							(*xtype) = 'C'
							Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (bsav), lda)
							//
							if *nofact {
								//
								//                          --- Test DPBSV  ---
								//
								//                          Compute the L*l' or U'*U factorization of the
								//                          matrix and solve the system.
								//
								Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (*kd)+1, n, (a), ldab, (afac), ldab)
								Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (x), lda)
								//
								(*srnamt) = *func() *[]byte {y := []byte("DPBSV "); return &y }()
								DPBSV(uplo, n, kd, (nrhs), (afac), ldab, (x), lda, info)
								//
								//                          Check error code from DPBSV .
								//
								if (*info) != (*izero) {
									Alaerh(path, func() *[]byte {y := []byte("DPBSV "); return &y }(), info, izero, uplo, n, n, kd, kd, (nrhs), imat, nfail, nerrs, (nout))
									goto Label40
								} else if (*info) != 0 {
									goto Label40
								}
								//
								//                          Re_construct matrix from factors and compute
								//                          residual.
								//
								Dpbt01(uplo, n, kd, (a), ldab, (afac), ldab, (rwork), &((*result)[0]))
								//
								//                          Compute residual of the computed solution.
								//
								Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (b), lda, (work), lda)
								Dpbt02(uplo, n, kd, (nrhs), (a), ldab, (x), lda, (work), lda, (rwork), &((*result)[1]))
								//
								//                          Check solution from generated exact solution.
								//
								Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
								(*nt) = 3
								//
								//                          Print information about the tests that did
								//                          not pass the threshold.
								//
								for (*k) = 1; (*k) <= (*nt); (*k)++ {
									if (*result)[(*k)-1] >= (*(thresh)) {
										if (*nfail) == 0 && (*nerrs) == 0 {
											Aladhd((nout), path)
										}
										WRITE((*(nout)), *func() *[]byte {
											y := []byte(" %s, uplo='%c', N =%5d, kd =%5d, type %1d, test(%1d)=%12.5f\n")
											return &y
										}(), *func() *[]byte {y := []byte("DPBSV "); return &y }(), (*uplo), (*n), (*kd), (*imat), (*k), (*result)[(*k)-1])
										(*nfail) = (*nfail) + 1
									}
									//Label30:
								}
								(*nrun) = (*nrun) + (*nt)
							Label40:
							}
							//
							//                       --- Test Dpbsvx ---
							//
							if !(*prefac) {
								Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (*kd)+1, n, zero, zero, (afac), ldab)
							}
							Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), zero, zero, (x), lda)
							if (*iequed) > 1 && (*n) > 0 {
								//
								//                          Equilibrate the matrix if fact='F' and
								//                          equed='Y'
								//
								Dlaqsb(uplo, n, kd, (a), ldab, (s), scond, amax, equed)
							}
							//
							//                       Solve the system and compute the condition
							//                       number and error bounds using Dpbsvx.
							//
							(*srnamt) = *func() *[]byte {y := []byte("Dpbsvx"); return &y }()
							Dpbsvx(fact, uplo, n, kd, (nrhs), (a), ldab, (afac), ldab, equed, (s), (b), lda, (x), lda, rcond, (rwork), &((*(rwork))[(*(nrhs))+0]), (work), (iwork), info)
							//
							//                       Check the error code from Dpbsvx.
							//
							if (*info) != (*izero) {
								Alaerh(path, func() *[]byte {y := []byte("Dpbsvx"); return &y }(), info, izero, append(append([]byte{}, (*fact)), (*uplo)), n, n, kd, kd, (nrhs), imat, nfail, nerrs, (nout))
								goto Label60
							}
							//
							if (*info) == 0 {
								if !(*prefac) {
									//
									//                             Re_construct matrix from factors and
									//                             compute residual.
									//
									Dpbt01(uplo, n, kd, (a), ldab, (afac), ldab, &((*(rwork))[2*(*(nrhs))+0]), &((*result)[0]))
									(*k1) = 1
								} else {
									(*k1) = 2
								}
								//
								//                          Compute residual of the computed solution.
								//
								Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, (nrhs), (bsav), lda, (work), lda)
								Dpbt02(uplo, n, kd, (nrhs), (asav), ldab, (x), lda, (work), lda, &((*(rwork))[2*(*(nrhs))+0]), &((*result)[1]))
								//
								//                          Check solution from generated exact solution.
								//
								if (*nofact) || ((*prefac) && blas.Lsame(equed, func() *byte {y := byte('N'); return &y }())) {
									Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
								} else {
									Dget04(n, (nrhs), (x), lda, (xact), lda, amax, &((*result)[2]))
								}
								//
								//                          Check the error bounds from iterative
								//                          refinement.
								//
								Dpbt05(uplo, n, kd, (nrhs), (asav), ldab, (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*(nrhs))+0]), &((*result)[3]))
							} else {
								(*k1) = 6
							}
							//
							//                       Compare rcond from Dpbsvx with the computed
							//                       value in rcondc.
							//
							(*result)[5] = (*Dget06(rcond, rcondc))
							//
							//                       Print information about the tests that did not
							//                       pass the threshold.
							//
							for (*k) = (*k1); (*k) <= 6; (*k)++ {
								if (*result)[(*k)-1] >= (*(thresh)) {
									if (*nfail) == 0 && (*nerrs) == 0 {
										Aladhd((nout), path)
									}
									if *prefac {
										WRITE((*(nout)), *func() *[]byte {
											y := []byte(" %s( '%c', '%c', %5d, %5d, ...), equed='%c', type %1d, test(%1d)=%12.5f\n")
											return &y
										}(), *func() *[]byte {y := []byte("Dpbsvx"); return &y }(), (*fact), (*uplo), (*n), (*kd), (*equed), (*imat), (*k), (*result)[(*k)-1])
									} else {
										WRITE((*(nout)), *func() *[]byte {
											y := []byte(" %s( '%c', '%c', %5d, %5d, ...), type %1d, test(%1d)=%12.5f\n")
											return &y
										}(), *func() *[]byte {y := []byte("Dpbsvx"); return &y }(), (*fact), (*uplo), (*n), (*kd), (*imat), (*k), (*result)[(*k)-1])
									}
									(*nfail) = (*nfail) + 1
								}
								//Label50:
							}
							(*nrun) = (*nrun) + 7 - (*k1)
						Label60:
						}
						//Label70:
					}
				Label80:
				}
				//Label90:
			}
			//Label100:
		}
		//Label110:
	}
	//
	//     Print a summary of the results.
	//
	Alasvm(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Ddrvpb
	//
}
