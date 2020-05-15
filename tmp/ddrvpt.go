package goblas

import 

// Ddrvpt tests Dptsv and -SVX.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvpt( dotype, nn, nval, nrhs, thresh, tsterr, a, d,
//                          e, b, x, xact, work, rwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       inTEGER            nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       inTEGER            nval(*)
//       DOUBLE PRECISION   a(*), B(*), d(*), E(*), rwork(*),
//      $                   work(*), X(*), XACt(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvpt tests Dptsv and -SVX.
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
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (nmax*2)
// \endverbatim
//
// \param[out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (nmax*2)
// \endverbatim
//
// \param[out] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (nmax*2)
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
//                      (nmax*max(3,nrhs))
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension
//                      (max(nmax,2*nrhs))
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
func Ddrvpt(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, a *[]float64, d *[]float64, e *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	zerot := new(bool)
	dist := new(byte)
	fact := new(byte)
	_type := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	Ia := new(int)
	ifact := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	ix := new(int)
	izero := new(int)
	j := new(int)
	k := new(int)
	k1 := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	cond := new(float64)
	dmax := new(float64)
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
	result := func() *[]float64 {
		arr := make([]float64, 6)
		return &arr
	}()
	z := func() *[]float64 {
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
	(*ntypes) = 12
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
	//     .. Intrinsic Functions ..
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
	//     .. Data statements ..
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3] = 0, 0, 0, 1
	//     ..
	//     .. Executable Statements ..
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y}()
	(*path)[1] = *func() *[]byte {y :=[]byte("PT"); return &y}()
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
	for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
		//
		//        Do for each value of N in nval.
		//
		(*n) = (*(nval))[(*in)-1]
		(*lda) = (MAX(1, (*n)))
		(*nimat) = (*ntypes)
		if (*n) <= 0 {
			(*nimat) = 1
		}
		//
		for (*imat) = 1; (*imat) <= (*nimat); (*imat)++ {
			//
			//           Do the tests only if dotype( imat) is true.
			//
			if (*n) > 0 && !(*(dotype))[(*imat)-1] {
				goto Label110
			}
			//
			//           Set up parameters with Dlatb4.
			//
			Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cond, dist)
			//
			(*zerot) = (*imat) >= 8 && (*imat) <= 10
			if (*imat) <= 6 {
				//
				//              Type 1-6:  generate a symmetric tridiagonal matrix of
				//              known condition number in lower triangular band storage.
				//
				(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y}()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cond, anorm, kl, ku, func() *byte {y := byte('B'); return &y}(), (a), func() *int {y := 2; return &y}(), (work), info)
				//
				//              Check the error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y}(), info, func() *int {y := 0; return &y}(), func() *byte {y := byte(' '); return &y}(), n, n, kl, ku, -1, imat, nfail, nerrs, (nout))
					goto Label110
				}
				(*izero) = 0
				//
				//              Copy the matrix to D and E.
				//
				(*ia) = 1
				for (*i) = 1; (*i) <= (*n)-1; (*i)++ {
					(*(d))[(*i)-1] = (*(a))[(*ia)-1]
					(*(e))[(*i)-1] = (*(a))[(*ia)+0]
					(*ia) = (*ia) + 2
					//Label20:
				}
				if (*n) > 0 {
					(*(d))[(*n)-1] = (*(a))[(*ia)-1]
				}
			} else {
				//
				//              Type 7-12:  generate a diagonally dominant matrix with
				//              unknown condition number in the vectors D and E.
				//
				if !(*zerot) || !(*(dotype))[6] {
					//
					//                 Let D and E have values from[-1,1].
					//
					Dlarnv(func() *int {y := 2; return &y}(), iseed, n, (d))
					Dlarnv(func() *int {y := 2; return &y}(), iseed, (*n)-1, (e))
					//
					//                 Make the tridiagonal matrix diagonally dominant.
					//
					if (*n) == 1 {
						(*(d))[0] = (ABS(((*(d))[0])))
					} else {
						(*(d))[0] = ABS(((*(d))[0])) + ABS(((*(e))[0]))
						(*(d))[(*n)-1] = ABS(((*(d))[(*n)-1])) + ABS(((*(e))[(*n)-0]))
						for (*i) = 2; (*i) <= (*n)-1; (*i)++ {
							(*(d))[(*i)-1] = ABS(((*(d))[(*i)-1])) + ABS(((*(e))[(*i)-1])) + ABS(((*(e))[(*i)-0]))
							//Label30:
						}
					}
					//
					//                 Scale D and E so the maximum element is anorm.
					//
					(*ix) = (*Idamax(n, (d), func() *int {y := 1; return &y}()))
					(*dmax) = (*(d))[(*ix)-1]
					Dscal(n, (*anorm)/(*dmax), (d), func() *int {y := 1; return &y}())
					if (*n) > 1 {
						Dscal((*n)-1, (*anorm)/(*dmax), (e), func() *int {y := 1; return &y}())
					}
					//
				} else if (*izero) > 0 {
					//
					//                 Reuse the last matrix by copying back the zeroed out
					//                 elements.
					//
					if (*izero) == 1 {
						(*(d))[0] = (*z)[1]
						if (*n) > 1 {
							(*(e))[0] = (*z)[2]
						}
					} else if (*izero) == (*n) {
						(*(e))[(*n)-0] = (*z)[0]
						(*(d))[(*n)-1] = (*z)[1]
					} else {
						(*(e))[(*izero)-0] = (*z)[0]
						(*(d))[(*izero)-1] = (*z)[1]
						(*(e))[(*izero)-1] = (*z)[2]
					}
				}
				//
				//              For types 8-10, set one row and column of the matrix to
				//              zero.
				//
				(*izero) = 0
				if (*imat) == 8 {
					(*izero) = 1
					(*z)[1] = (*(d))[0]
					(*(d))[0] = (*zero)
					if (*n) > 1 {
						(*z)[2] = (*(e))[0]
						(*(e))[0] = (*zero)
					}
				} else if (*imat) == 9 {
					(*izero) = (*n)
					if (*n) > 1 {
						(*z)[0] = (*(e))[(*n)-0]
						(*(e))[(*n)-0] = (*zero)
					}
					(*z)[1] = (*(d))[(*n)-1]
					(*(d))[(*n)-1] = (*zero)
				} else if (*imat) == 10 {
					(*izero) = ((*n) + 1) / 2
					if (*izero) > 1 {
						(*z)[0] = (*(e))[(*izero)-0]
						(*z)[2] = (*(e))[(*izero)-1]
						(*(e))[(*izero)-0] = (*zero)
						(*(e))[(*izero)-1] = (*zero)
					}
					(*z)[1] = (*(d))[(*izero)-1]
					(*(d))[(*izero)-1] = (*zero)
				}
			}
			//
			//           Generate nrhs random solution vectors.
			//
			(*ix) = 1
			for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
				Dlarnv(func() *int {y := 2; return &y}(), iseed, n, &((*(xact))[(*ix)-1]))
				(*ix) = (*ix) + (*lda)
				//Label40:
			}
			//
			//           Set the right hand side.
			//
			Dlaptm(n, (nrhs), one, (d), (e), (xact), lda, zero, (b), lda)
			//
			for (*ifact) = 1; (*ifact) <= 2; (*ifact)++ {
				if (*ifact) == 1 {
					(*fact) = 'F'
				} else {
					(*fact) = 'N'
				}
				//
				//              Compute the condition number for comparison with
				//              the value returned by Dptsvx.
				//
				if *zerot {
					if (*ifact) == 1 {
						goto Label100
					}
					(*rcondc) = (*zero)
					//
				} else if (*ifact) == 1 {
					//
					//                 Compute the 1-norm of A.
					//
					(*anorm) = (*Dlanst(func() *byte {y := byte('1'); return &y}(), n, (d), (e)))
					//
					Dcopy(n, (d), func() *int {y := 1; return &y}(), &((*(d))[(*n)+0]), func() *int {y := 1; return &y}())
					if (*n) > 1 {
						Dcopy((*n)-1, (e), func() *int {y := 1; return &y}(), &((*(e))[(*n)+0]), func() *int {y := 1; return &y}())
					}
					//
					//                 Factor the matrix A.
					//
					Dpttrf(n, &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), info)
					//
					//                 Use Dpttrs to solve for one column at a time of
					//                 inv(a), computing the maximum column sum as we go.
					//
					(*ainvnm) = (*zero)
					for (*i) = 1; (*i) <= (*n); (*i)++ {
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							(*(x))[(*j)-1] = (*zero)
							//Label50:
						}
						(*(x))[(*i)-1] = (*one)
						Dpttrs(n, func() *int {y := 1; return &y}(), &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), (x), lda, info)
						(*ainvnm) = (MAX((*ainvnm), Dasum(n, (x), func() *int {y := 1; return &y}())))
						//Label60:
					}
					//
					//                 Compute the 1-norm condition number of A.
					//
					if (*anorm) <= (*zero) || (*ainvnm) <= (*zero) {
						(*rcondc) = (*one)
					} else {
						(*rcondc) = ((*one) / (*anorm)) / (*ainvnm)
					}
				}
				//
				if (*ifact) == 2 {
					//
					//                 --- Test Dptsv --
					//
					Dcopy(n, (d), func() *int {y := 1; return &y}(), &((*(d))[(*n)+0]), func() *int {y := 1; return &y}())
					if (*n) > 1 {
						Dcopy((*n)-1, (e), func() *int {y := 1; return &y}(), &((*(e))[(*n)+0]), func() *int {y := 1; return &y}())
					}
					Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (x), lda)
					//
					//                 Factor A as L*D*L' and solve the system A*X = B.
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dptsv "); return &y}()
					Dptsv(n, (nrhs), &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), (x), lda, info)
					//
					//                 Check error code from Dptsv .
					//
					if (*info) != (*izero) {
						Alaerh(path, func() *[]byte {y :=[]byte("Dptsv "); return &y}(), info, izero, func() *byte {y := byte(' '); return &y}(), n, n, func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), (nrhs), imat, nfail, nerrs, (nout))
					}
					(*nt) = 0
					if (*izero) == 0 {
						//
						//                    Check the factorization by computing the ratio
						//                       norm(L*D*L' - A) / (n * norm(a) * eps)
						//
						Dptt01(n, (d), (e), &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), (work), &((*result)[0]))
						//
						//                    Compute the residual in the solution.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
						Dptt02(n, (nrhs), (d), (e), (x), lda, (work), lda, &((*result)[1]))
						//
						//                    Check solution from generated exact solution.
						//
						Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
						(*nt) = 3
					}
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = 1; (*k) <= (*nt); (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Aladhd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %s, N =%5d, type %2d, test %2d, ratio = %12.5f\n"); return &y}(), *func() *[]byte {y :=[]byte("Dptsv "); return &y}(), (*n), (*imat), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label70:
					}
					(*nrun) = (*nrun) + (*nt)
				}
				//
				//              --- Test Dptsvx ---
				//
				if (*ifact) > 1 {
					//
					//                 Initialize d( N+1:2*n) and E( N+1:2*n) to zero.
					//
					for (*i) = 1; (*i) <= (*n)-1; (*i)++ {
						(*(d))[(*n)+(*i)-1] = (*zero)
						(*(e))[(*n)+(*i)-1] = (*zero)
						//Label80:
					}
					if (*n) > 0 {
						(*(d))[(*n)+(*n)-1] = (*zero)
					}
				}
				//
				Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), zero, zero, (x), lda)
				//
				//              Solve the system and compute the condition number and
				//              error bounds using Dptsvx.
				//
				(*srnamt) = *func() *[]byte {y :=[]byte("Dptsvx"); return &y}()
				Dptsvx(fact, n, (nrhs), (d), (e), &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), (b), lda, (x), lda, rcond, (rwork), &((*(rwork))[(*(nrhs))+0]), (work), info)
				//
				//              Check the error code from Dptsvx.
				//
				if (*info) != (*izero) {
					Alaerh(path, func() *[]byte {y :=[]byte("Dptsvx"); return &y}(), info, izero, fact, n, n, func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), (nrhs), imat, nfail, nerrs, (nout))
				}
				if (*izero) == 0 {
					if (*ifact) == 2 {
						//
						//                    Check the factorization by computing the ratio
						//                       norm(L*D*L' - A) / (n * norm(a) * eps)
						//
						(*k1) = 1
						Dptt01(n, (d), (e), &((*(d))[(*n)+0]), &((*(e))[(*n)+0]), (work), &((*result)[0]))
					} else {
						(*k1) = 2
					}
					//
					//                 Compute the residual in the solution.
					//
					Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
					Dptt02(n, (nrhs), (d), (e), (x), lda, (work), lda, &((*result)[1]))
					//
					//                 Check solution from generated exact solution.
					//
					Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
					//
					//                 Check error bounds from iterative refinement.
					//
					Dptt05(n, (nrhs), (d), (e), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*(nrhs))+0]), &((*result)[3]))
				} else {
					(*k1) = 6
				}
				//
				//              Check the reciprocal of the condition number.
				//
				(*result)[5] = (*Dget06(rcond, rcondc))
				//
				//              Print information about the tests that did not pass
				//              the threshold.
				//
				for (*k) = (*k1); (*k) <= 6; (*k)++ {
					if (*result)[(*k)-1] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Aladhd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %s, fact='%c', N =%5d, type %2d, test %2d, ratio = %12.5f\n"); return &y}(), *func() *[]byte {y :=[]byte("Dptsvx"); return &y}(), (*fact), (*n), (*imat), (*k), (*result)[(*k)-1])
						(*nfail) = (*nfail) + 1
					}
					//Label90:
				}
				(*nrun) = (*nrun) + 7 - (*k1)
			Label100:
			}
		Label110:
		}
		//Label120:
	}
	//
	//     Print a summary of the results.
	//
	Alasvm(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Ddrvpt
	//
}
