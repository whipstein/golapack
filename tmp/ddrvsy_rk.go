package goblas

import 

// DdrvsyRk tests the driver routines DsysvRk.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//     SUBROUTinE DdrvsyRk( dotype, nn, nval, nrhs, thresh, tsterr,
//    $                      nmax, a, afac, e, ainv, b, x, xact, work,
//    $                      rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nmax, nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), nval(*)
//       DOUBLE PRECISION   a(*), afac(*), ainv(*), B(*), E(*),
//      $                   rwork(*), work(*), X(*), xact(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
// DdrvsyRk tests the driver routines DsysvRk.
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
// \param[out] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (nmax)
// \endverbatim
//
// \param[out] ainv
// \verbatim
//          ainv is DOUBLE PRECISION array, dimension (nmax*nmax)
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
//          work is DOUBLE PRECISION array, dimension (nmax*max(2,nrhs))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (nmax+2*nrhs)
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
func DdrvsyRk(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, e *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	nfact := new(int)
	zerot := new(bool)
	dist := new(byte)
	fact := new(byte)
	_type := new(byte)
	uplo := new(byte)
	xtype := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	matpath := func() *[]byte {
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
	kl := new(int)
	ku := new(int)
	lda := new(int)
	lwork := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nbmin := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	cndnum := new(float64)
	rcondc := new(float64)
	facts := func() *[]byte {
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
		arr := make([]float64, 3)
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
	(*ntests) = 3
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
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
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
	(*uplos)[0], (*uplos)[1], (*facts)[0], (*facts)[1] = 'U', 'L', 'F', 'N'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	//     Test path
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y :=[]byte("SK"); return &y }()
	//
	//     Path to generate matrices
	//
	(*matpath)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*matpath)[1] = *func() *[]byte {y :=[]byte("SY"); return &y }()
	//
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
	//     Set the block size and minimum block size for which the block
	//     routine should be used, which will be later returned by ILAENV.
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
				(*uplo) = (*uplos)[(*iuplo)-1]
				//
				//              Begin generate the test matrix A.
				//
				//              Set up parameters with Dlatb4 for the matrix generator
				//              based on the type of matrix to be generated.
				//
				Dlatb4(matpath, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
				//
				//              Generate a matrix with Dlatms.
				//
				(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y }()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, uplo, (a), lda, (work), info)
				//
				//              Check error code from Dlatms and handle error.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					//
					//                 Skip all tests for this generated matrix
					//
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
							(*ioff) = ((*izero) - 1) * (*lda)
							for (*i) = 1; (*i) <= (*izero)-1; (*i)++ {
								(*(a))[(*ioff)+(*i)-1] = (*zero)
								//Label20:
							}
							(*ioff) = (*ioff) + (*izero)
							for (*i) = (*izero); (*i) <= (*n); (*i)++ {
								(*(a))[(*ioff)-1] = (*zero)
								(*ioff) = (*ioff) + (*lda)
								//Label30:
							}
						} else {
							(*ioff) = (*izero)
							for (*i) = 1; (*i) <= (*izero)-1; (*i)++ {
								(*(a))[(*ioff)-1] = (*zero)
								(*ioff) = (*ioff) + (*lda)
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
								(*ioff) = (*ioff) + (*lda)
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
								(*ioff) = (*ioff) + (*lda)
								//Label90:
							}
						}
					}
				} else {
					(*izero) = 0
				}
				//
				//              End generate the test matrix A.
				//
				for (*ifact) = 1; (*ifact) <= (*nfact); (*ifact)++ {
					//
					//                 Do first for fact = 'F', then for other values.
					//
					(*fact) = (*facts)[(*ifact)-1]
					//
					//                 Compute the condition number
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
						(*anorm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), uplo, n, (a), lda, (rwork)))
						//
						//                    Factor the matrix A.
						//
						Dlacpy(uplo, n, n, (a), lda, (afac), lda)
						DsytrfRk(uplo, n, (afac), lda, (e), (iwork), (work), lwork, info)
						//
						//                    Compute inv(a) and take its norm.
						//
						Dlacpy(uplo, n, n, (afac), lda, (ainv), lda)
						(*lwork) = ((*n) + (*nb) + 1) * ((*nb) + 3)
						//
						//                    We need to copute the invesrse to compute
						//                    rcondc that is used later in TEST3.
						//
						Dsytri3(uplo, n, (ainv), lda, (e), (iwork), (work), lwork, info)
						(*ainvnm) = (*Dlansy(func() *byte {y := byte('1'); return &y }(), uplo, n, (ainv), lda, (rwork)))
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
					(*srnamt) = *func() *[]byte {y :=[]byte("Dlarhs"); return &y }()
					Dlarhs(matpath, xtype, uplo, func() *byte {y := byte(' '); return &y }(), n, n, kl, ku, (nrhs), (a), lda, (xact), lda, (b), lda, iseed, info)
					(*xtype) = 'C'
					//
					//                 --- Test DsysvRk  ---
					//
					if (*ifact) == 2 {
						Dlacpy(uplo, n, n, (a), lda, (afac), lda)
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, (nrhs), (b), lda, (x), lda)
						//
						//                    Factor the matrix and solve the system using
						//                    DsysvRk.
						//
						(*srnamt) = *func() *[]byte {y :=[]byte("DsysvRk"); return &y }()
						DsysvRk(uplo, n, (nrhs), (afac), lda, (e), (iwork), (x), lda, (work), lwork, info)
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
						//                    Check error code from DsysvRk and handle error.
						//
						if (*info) != (*k) {
							Alaerh(path, func() *[]byte {y :=[]byte("DsysvRk"); return &y }(), info, k, uplo, n, n, -1, -1, (nrhs), imat, nfail, nerrs, (nout))
							goto Label120
						} else if (*info) != 0 {
							goto Label120
						}
						//
						//+    TEST 1      Re_construct matrix from factors and compute
						//                 residual.
						//
						Dsyt01_3(uplo, n, (a), lda, (afac), lda, (e), (iwork), (ainv), lda, (rwork), &((*result)[0]))
						//
						//+    TEST 2      Compute residual of the computed solution.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y }(), n, (nrhs), (b), lda, (work), lda)
						Dpot02(uplo, n, (nrhs), (a), lda, (x), lda, (work), lda, (rwork), &((*result)[1]))
						//
						//+    TEST 3
						//                 Check solution from generated exact solution.
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
								WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %s, uplo='%c', N =%5d, type %2d, test %2d, ratio =%12.5f\n"); return &y }(), *func() *[]byte {y :=[]byte("DsysvRk"); return &y }(), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
								(*nfail) = (*nfail) + 1
							}
							//Label110:
						}
						(*nrun) = (*nrun) + (*nt)
					Label120:
					}
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
	//     End of DdrvsyRk
	//
}
