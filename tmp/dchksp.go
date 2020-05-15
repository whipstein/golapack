package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dchksp tests Dsptrf, -tri, -TRS, -RFS, and -CON
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchksp( dotype, nn, nval, nns, nsval, thresh, tsterr,
//                          nmax, a, afac, ainv, b, x, xact, work, rwork,
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
// Dchksp tests Dsptrf, -tri, -TRS, -RFS, and -CON
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
//                      (nmax*max(2,NSMAX))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array,
//                                 dimension (nmax+2*nsmax)
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
func Dchksp(dotype *[]bool, nn *int, nval *[]int, nns *int, nsval *[]int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, ainv *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	trfcon := new(bool)
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
	imat := new(int)
	in := new(int)
	info := new(int)
	ioff := new(int)
	irhs := new(int)
	iuplo := new(int)
	izero := new(int)
	j := new(int)
	k := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	npp := new(int)
	nrhs := new(int)
	nrun := new(int)
	nt := new(int)
	anorm := new(float64)
	cndnum := new(float64)
	rcond := new(float64)
	rcondc := new(float64)
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
	(*zero) = 0.0e+0
	(*ntypes) = 10
	(*ntests) = 8
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
	(*uplos)[0], (*uplos)[1] = 'U', 'L'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y := []byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y := []byte("SP"); return &y }()
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
		Derrsy(path, (nout))
	}
	(*infot) = 0
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
		(*izero) = 0
		for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
			//
			//           Do the tests only if dotype( imat) is true.
			//
			if !(*(dotype))[(*imat)-1] {
				goto Label160
			}
			//
			//           Skip types 3, 4, 5, or 6 if the matrix size is too small.
			//
			(*zerot) = (*imat) >= 3 && (*imat) <= 6
			if (*zerot) && (*n) < (*imat)-2 {
				goto Label160
			}
			//
			//           Do first for uplo = 'U', then for uplo = 'L'
			//
			for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
				(*uplo) = (*uplos)[(*iuplo)-1]
				if blas.Lsame(uplo, func() *byte {y := byte('U'); return &y }()) {
					(*packit) = 'C'
				} else {
					(*packit) = 'R'
				}
				//
				//              Set up parameters with Dlatb4 and generate a test matrix
				//              with Dlatms.
				//
				Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cndnum, dist)
				//
				(*srnamt) = *func() *[]byte {y := []byte("Dlatms"); return &y }()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, kl, ku, packit, (a), lda, (work), info)
				//
				//              Check error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y := []byte("Dlatms"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					goto Label150
				}
				//
				//              For types 3-6, zero one or more rows and columns of
				//              the matrix to test that info is returned correctly.
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
				//              Compute the L*D*l' or U*D*U' factorization of the matrix.
				//
				(*npp) = (*n) * ((*n) + 1) / 2
				Dcopy(npp, (a), func() *int {y := 1; return &y }(), (afac), func() *int {y := 1; return &y }())
				(*srnamt) = *func() *[]byte {y := []byte("Dsptrf"); return &y }()
				Dsptrf(uplo, n, (afac), (iwork), info)
				//
				//              Adjust the expected value of info to account for
				//              pivoting.
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
				//              Check error code from Dsptrf.
				//
				if (*info) != (*k) {
					Alaerh(path, func() *[]byte {y := []byte("Dsptrf"); return &y }(), info, k, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
				}
				if (*info) != 0 {
					(*trfcon) = true
				} else {
					(*trfcon) = false
				}
				//
				//+    TEST 1
				//              Re_construct matrix from factors and compute residual.
				//
				Dspt01(uplo, n, (a), (afac), (iwork), (ainv), lda, (rwork), &((*result)[0]))
				(*nt) = 1
				//
				//+    TEST 2
				//              Form the inverse and compute the residual.
				//
				if !(*trfcon) {
					Dcopy(npp, (afac), func() *int {y := 1; return &y }(), (ainv), func() *int {y := 1; return &y }())
					(*srnamt) = *func() *[]byte {y := []byte("Dsptri"); return &y }()
					Dsptri(uplo, n, (ainv), (iwork), (work), info)
					//
					//              Check error code from Dsptri.
					//
					if (*info) != 0 {
						Alaerh(path, func() *[]byte {y := []byte("Dsptri"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
					}
					//
					Dppt03(uplo, n, (a), (ainv), (work), lda, (rwork), rcondc, &((*result)[1]))
					(*nt) = 2
				}
				//
				//              Print information about the tests that did not pass
				//              the threshold.
				//
				for (*k) = 1; (*k) <= (*nt); (*k)++ {
					if (*result)[(*k)-1] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Alahd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte {y := []byte(" uplo = '%c', N =%5d, type %2d, test %2d, ratio =%12.5f\n"); return &y }(), (*uplo), (*n), (*imat), (*k), (*result)[(*k)-1])
						(*nfail) = (*nfail) + 1
					}
					//Label110:
				}
				(*nrun) = (*nrun) + (*nt)
				//
				//              Do only the condition estimate if info is not 0.
				//
				if *trfcon {
					(*rcondc) = (*zero)
					goto Label140
				}
				//
				for (*irhs) = 1; (*irhs) <= (*(nns)); (*irhs)++ {
					(*nrhs) = (*(nsval))[(*irhs)-1]
					//
					//+    TEST 3
					//              Solve and compute residual for  A * X = B.
					//
					(*srnamt) = *func() *[]byte {y := []byte("Dlarhs"); return &y }()
					Dlarhs(path, xtype, uplo, func() *byte {y := byte(' '); return &y }(), n, n, kl, ku, nrhs, (a), lda, (xact), lda, (b), lda, iseed, info)
					Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, nrhs, (b), lda, (x), lda)
					//
					(*srnamt) = *func() *[]byte {y := []byte("DSPTRS"); return &y }()
					Dsptrs(uplo, n, nrhs, (afac), (iwork), (x), lda, info)
					//
					//              Check error code from DSPTRS.
					//
					if (*info) != 0 {
						Alaerh(path, func() *[]byte {y := []byte("DSPTRS"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, (nout))
					}
					//
					Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), n, nrhs, (b), lda, (work), lda)
					Dppt02(uplo, n, nrhs, (a), (x), lda, (work), lda, (rwork), &((*result)[2]))
					//
					//+    TEST 4
					//              Check solution from generated exact solution.
					//
					Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[3]))
					//
					//+    TESts 5, 6, and 7
					//              Use iterative refinement to improve the solution.
					//
					(*srnamt) = *func() *[]byte {y := []byte("DSPRFS"); return &y }()
					Dsprfs(uplo, n, nrhs, (a), (afac), (iwork), (b), lda, (x), lda, (rwork), &((*(rwork))[(*nrhs)+0]), (work), &((*(iwork))[(*n)+0]), info)
					//
					//              Check error code from DSPRFS.
					//
					if (*info) != 0 {
						Alaerh(path, func() *[]byte {y := []byte("DSPRFS"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, (nout))
					}
					//
					Dget04(n, nrhs, (x), lda, (xact), lda, rcondc, &((*result)[4]))
					Dppt05(uplo, n, nrhs, (a), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*nrhs)+0]), &((*result)[5]))
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = 3; (*k) <= 7; (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y := []byte(" uplo = '%c', N =%5d, nrhs=%3d, type %2d, test(%2d) =%12.5f\n")
								return &y
							}(), (*uplo), (*n), (*nrhs), (*imat), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label120:
					}
					(*nrun) = (*nrun) + 5
					//Label130:
				}
				//
				//+    TEST 8
				//              Get an estimate of rcond = 1/cndnum.
				//
			Label140:
				;
				(*anorm) = (*Dlansp(func() *byte {y := byte('1'); return &y }(), uplo, n, (a), (rwork)))
				(*srnamt) = *func() *[]byte {y := []byte("DSPCON"); return &y }()
				Dspcon(uplo, n, (afac), (iwork), anorm, rcond, (work), &((*(iwork))[(*n)+0]), info)
				//
				//              Check error code from DSPCON.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y := []byte("DSPCON"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
				}
				//
				(*result)[7] = (*Dget06(rcond, rcondc))
				//
				//              Print the test ratio if it is .GE. thresh.
				//
				if (*result)[7] >= (*(thresh)) {
					if (*nfail) == 0 && (*nerrs) == 0 {
						Alahd((nout), path)
					}
					WRITE((*(nout)), *func() *[]byte {y := []byte(" uplo = '%c', N =%5d, type %2d, test %2d, ratio =%12.5f\n"); return &y }(), (*uplo), (*n), (*imat), 8, (*result)[7])
					(*nfail) = (*nfail) + 1
				}
				(*nrun) = (*nrun) + 1
			Label150:
			}
		Label160:
		}
		//Label170:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchksp
	//
}
