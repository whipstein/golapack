package goblas

import 

// Ddrvgt tests Dgtsv and -SVX.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Ddrvgt( dotype, nn, nval, nrhs, thresh, tsterr, a, af,
//                          b, x, xact, work, rwork, iwork, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       inTEGER            nn, nout, nrhs
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       inTEGER            iwork(*), nval(*)
//       DOUBLE PRECISION   a(*), af(*), B(*), rwork(*), work(*),
//      $                   X(*), XACt(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ddrvgt tests Dgtsv and -SVX.
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
//          The number of right hand sides, nrhs >= 0.
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
//          A is DOUBLE PRECISION array, dimension (nmax*4)
// \endverbatim
//
// \param[out] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (nmax*4)
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
func Ddrvgt(dotype *[]bool, nn *int, nval *[]int, nrhs *int, thresh *float64, tsterr *bool, a *[]float64, af *[]float64, b *[]float64, x *[]float64, xact *[]float64, work *[]float64, rwork *[]float64, iwork *[]int, nout *int) {
	one := new(float64)
	zero := new(float64)
	ntypes := new(int)
	ntests := new(int)
	trfcon := new(bool)
	zerot := new(bool)
	dist := new(byte)
	fact := new(byte)
	trans := new(byte)
	_type := new(byte)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	ifact := new(int)
	imat := new(int)
	in := new(int)
	info := new(int)
	itran := new(int)
	ix := new(int)
	izero := new(int)
	j := new(int)
	k := new(int)
	k1 := new(int)
	kl := new(int)
	trfcon := new(int)
	ku := new(int)
	lda := new(int)
	m := new(int)
	mode := new(int)
	n := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nrun := new(int)
	nt := new(int)
	ainvnm := new(float64)
	anorm := new(float64)
	anormi := new(float64)
	anormo := new(float64)
	cond := new(float64)
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
	(*iseedy)[0], (*iseedy)[1], (*iseedy)[2], (*iseedy)[3], (*transs)[0], (*transs)[1], (*transs)[2] = 0, 0, 0, 1, 'N', 'T', 'C'
	//     ..
	//     .. Executable Statements ..
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y}()
	(*path)[1] = *func() *[]byte {y :=[]byte("GT"); return &y}()
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
		(*m) = (MAX((*n)-1, 0))
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
			if !(*(dotype))[(*imat)-1] {
				goto Label130
			}
			//
			//           Set up parameters with Dlatb4.
			//
			Dlatb4(path, imat, n, n, _type, kl, ku, anorm, mode, cond, dist)
			//
			(*zerot) = (*imat) >= 8 && (*imat) <= 10
			if (*imat) <= 6 {
				//
				//              Types 1-6:  generate matrices of known condition number.
				//
				(*koff) = (MAX(2-(*ku), 3-MAX(1, (*n))))
				(*srnamt) = *func() *[]byte {y :=[]byte("Dlatms"); return &y}()
				Dlatms(n, n, dist, iseed, _type, (rwork), mode, cond, anorm, kl, ku, func() *byte {y := byte('Z'); return &y}(), &((*(af))[(*koff)-1]), func() *int {y := 3; return &y}(), (work), info)
				//
				//              Check the error code from Dlatms.
				//
				if (*info) != 0 {
					Alaerh(path, func() *[]byte {y :=[]byte("Dlatms"); return &y}(), info, func() *int {y := 0; return &y}(), func() *byte {y := byte(' '); return &y}(), n, n, kl, ku, -1, imat, nfail, nerrs, (nout))
					goto Label130
				}
				(*izero) = 0
				//
				if (*n) > 1 {
					Dcopy((*n)-1, &((*(af))[3]), func() *int {y := 3; return &y}(), (a), func() *int {y := 1; return &y}())
					Dcopy((*n)-1, &((*(af))[2]), func() *int {y := 3; return &y}(), &((*(a))[(*n)+(*m)+0]), func() *int {y := 1; return &y}())
				}
				Dcopy(n, &((*(af))[1]), func() *int {y := 3; return &y}(), &((*(a))[(*m)+0]), func() *int {y := 1; return &y}())
			} else {
				//
				//              Types 7-12:  generate tridiagonal matrices with
				//              unknown condition numbers.
				//
				if !(*zerot) || !(*(dotype))[6] {
					//
					//                 Generate a matrix with elements from[-1,1].
					//
					Dlarnv(func() *int {y := 2; return &y}(), iseed, (*n)+2*(*m), (a))
					if (*anorm) != (*one) {
						Dscal((*n)+2*(*m), anorm, (a), func() *int {y := 1; return &y}())
					}
				} else if (*izero) > 0 {
					//
					//                 Reuse the last matrix by copying back the zeroed out
					//                 elements.
					//
					if (*izero) == 1 {
						(*(a))[(*n)-1] = (*z)[1]
						if (*n) > 1 {
							(*(a))[0] = (*z)[2]
						}
					} else if (*izero) == (*n) {
						(*(a))[3*(*n)-1] = (*z)[0]
						(*(a))[2*(*n)-0] = (*z)[1]
					} else {
						(*(a))[2*(*n)-2+(*izero)-1] = (*z)[0]
						(*(a))[(*n)-1+(*izero)-1] = (*z)[1]
						(*(a))[(*izero)-1] = (*z)[2]
					}
				}
				//
				//              If imat > 7, set one column of the matrix to 0.
				//
				if !(*zerot) {
					(*izero) = 0
				} else if (*imat) == 8 {
					(*izero) = 1
					(*z)[1] = (*(a))[(*n)-1]
					(*(a))[(*n)-1] = (*zero)
					if (*n) > 1 {
						(*z)[2] = (*(a))[0]
						(*(a))[0] = (*zero)
					}
				} else if (*imat) == 9 {
					(*izero) = (*n)
					(*z)[0] = (*(a))[3*(*n)-1]
					(*z)[1] = (*(a))[2*(*n)-0]
					(*(a))[3*(*n)-1] = (*zero)
					(*(a))[2*(*n)-0] = (*zero)
				} else {
					(*izero) = ((*n) + 1) / 2
					for (*i) = (*izero); (*i) <= (*n)-1; (*i)++ {
						(*(a))[2*(*n)-2+(*i)-1] = (*zero)
						(*(a))[(*n)-1+(*i)-1] = (*zero)
						(*(a))[(*i)-1] = (*zero)
						//Label20:
					}
					(*(a))[3*(*n)-1] = (*zero)
					(*(a))[2*(*n)-0] = (*zero)
				}
			}
			//
			for (*ifact) = 1; (*ifact) <= 2; (*ifact)++ {
				if (*ifact) == 1 {
					(*fact) = 'F'
				} else {
					(*fact) = 'N'
				}
				//
				//              Compute the condition number for comparison with
				//              the value returned by Dgtsvx.
				//
				if *zerot {
					if (*ifact) == 1 {
						goto Label120
					}
					(*rcondo) = (*zero)
					(*rcondi) = (*zero)
					//
				} else if (*ifact) == 1 {
					Dcopy((*n)+2*(*m), (a), func() *int {y := 1; return &y}(), (af), func() *int {y := 1; return &y}())
					//
					//                 Compute the 1-norm and infinity-norm of A.
					//
					(*anormo) = (*Dlangt(func() *byte {y := byte('1'); return &y}(), n, (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0])))
					(*anormi) = (*Dlangt(func() *byte {y := byte('I'); return &y}(), n, (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0])))
					//
					//                 Factor the matrix A.
					//
					Dgttrf(n, (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), &((*(af))[(*n)+2*(*m)+0]), (iwork), info)
					//
					//                 Use Dgttrs to solve for one column at a time of
					//                 inv(a), computing the maximum column sum as we go.
					//
					(*ainvnm) = (*zero)
					for (*i) = 1; (*i) <= (*n); (*i)++ {
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							(*(x))[(*j)-1] = (*zero)
							//Label30:
						}
						(*(x))[(*i)-1] = (*one)
						Dgttrs(func() *[]byte {y :=[]byte("No transpose"); return &y}(), n, func() *int {y := 1; return &y}(), (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), &((*(af))[(*n)+2*(*m)+0]), (iwork), (x), lda, info)
						(*ainvnm) = (MAX((*ainvnm), Dasum(n, (x), func() *int {y := 1; return &y}())))
						//Label40:
					}
					//
					//                 Compute the 1-norm condition number of A.
					//
					if (*anormo) <= (*zero) || (*ainvnm) <= (*zero) {
						(*rcondo) = (*one)
					} else {
						(*rcondo) = ((*one) / (*anormo)) / (*ainvnm)
					}
					//
					//                 Use Dgttrs to solve for one column at a time of
					//                 inv(A'), computing the maximum column sum as we go.
					//
					(*ainvnm) = (*zero)
					for (*i) = 1; (*i) <= (*n); (*i)++ {
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							(*(x))[(*j)-1] = (*zero)
							//Label50:
						}
						(*(x))[(*i)-1] = (*one)
						Dgttrs(func() *[]byte {y :=[]byte("Transpose"); return &y}(), n, func() *int {y := 1; return &y}(), (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), &((*(af))[(*n)+2*(*m)+0]), (iwork), (x), lda, info)
						(*ainvnm) = (MAX((*ainvnm), Dasum(n, (x), func() *int {y := 1; return &y}())))
						//Label60:
					}
					//
					//                 Compute the infinity-norm condition number of A.
					//
					if (*anormi) <= (*zero) || (*ainvnm) <= (*zero) {
						(*rcondi) = (*one)
					} else {
						(*rcondi) = ((*one) / (*anormi)) / (*ainvnm)
					}
				}
				//
				for (*itran) = 1; (*itran) <= 3; (*itran)++ {
					(*trans) = (*transs)[(*itran)-1]
					if (*itran) == 1 {
						(*rcondc) = (*rcondo)
					} else {
						(*rcondc) = (*rcondi)
					}
					//
					//                 Generate nrhs random solution vectors.
					//
					(*ix) = 1
					for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
						Dlarnv(func() *int {y := 2; return &y}(), iseed, n, &((*(xact))[(*ix)-1]))
						(*ix) = (*ix) + (*lda)
						//Label70:
					}
					//
					//                 Set the right hand side.
					//
					Dlagtm(trans, n, (nrhs), one, (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (xact), lda, zero, (b), lda)
					//
					if (*ifact) == 2 && (*itran) == 1 {
						//
						//                    --- Test Dgtsv  ---
						//
						//                    Solve the system using Gaussian elimination with
						//                    partial pivoting.
						//
						Dcopy((*n)+2*(*m), (a), func() *int {y := 1; return &y}(), (af), func() *int {y := 1; return &y}())
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (x), lda)
						//
						(*srnamt) = *func() *[]byte {y :=[]byte("Dgtsv "); return &y}()
						Dgtsv(n, (nrhs), (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), (x), lda, info)
						//
						//                    Check error code from Dgtsv .
						//
						if (*info) != (*izero) {
							Alaerh(path, func() *[]byte {y :=[]byte("Dgtsv "); return &y}(), info, izero, func() *byte {y := byte(' '); return &y}(), n, n, func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), (nrhs), imat, nfail, nerrs, (nout))
						}
						(*nt) = 1
						if (*izero) == 0 {
							//
							//                       Check residual of computed solution.
							//
							Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
							Dgtt02(trans, n, (nrhs), (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (x), lda, (work), lda, &((*result)[1]))
							//
							//                       Check solution from generated exact solution.
							//
							Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
							(*nt) = 3
						}
						//
						//                    Print information about the tests that did not pass
						//                    the threshold.
						//
						for (*k) = 2; (*k) <= (*nt); (*k)++ {
							if (*result)[(*k)-1] >= (*(thresh)) {
								if (*nfail) == 0 && (*nerrs) == 0 {
									Aladhd((nout), path)
								}
								WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %s, N =%5d, type %2d, test %2d, ratio = %12.5f\n"); return &y}(), *func() *[]byte {y :=[]byte("Dgtsv "); return &y}(), (*n), (*imat), (*k), (*result)[(*k)-1])
								(*nfail) = (*nfail) + 1
							}
							//Label80:
						}
						(*nrun) = (*nrun) + (*nt) - 1
					}
					//
					//                 --- Test Dgtsvx ---
					//
					if (*ifact) > 1 {
						//
						//                    Initialize af to zero.
						//
						for (*i) = 1; (*i) <= 3*(*n)-2; (*i)++ {
							(*(af))[(*i)-1] = (*zero)
							//Label90:
						}
					}
					Dlaset(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), zero, zero, (x), lda)
					//
					//                 Solve the system and compute the condition number and
					//                 error bounds using Dgtsvx.
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dgtsvx"); return &y}()
					Dgtsvx(fact, trans, n, (nrhs), (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), &((*(af))[(*n)+2*(*m)+0]), (iwork), (b), lda, (x), lda, rcond, (rwork), &((*(rwork))[(*(nrhs))+0]), (work), &((*(iwork))[(*n)+0]), info)
					//
					//                 Check the error code from Dgtsvx.
					//
					if (*info) != (*izero) {
						Alaerh(path, func() *[]byte {y :=[]byte("Dgtsvx"); return &y}(), info, izero, append(append([]byte{}, (*fact)), (*trans)), n, n, func() *int {y := 1; return &y}(), func() *int {y := 1; return &y}(), (nrhs), imat, nfail, nerrs, (nout))
					}
					//
					if (*ifact) >= 2 {
						//
						//                    Re_construct matrix from factors and compute
						//                    residual.
						//
						Dgtt01(n, (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (af), &((*(af))[(*m)+0]), &((*(af))[(*n)+(*m)+0]), &((*(af))[(*n)+2*(*m)+0]), (iwork), (work), lda, (rwork), &((*result)[0]))
						(*k1) = 1
					} else {
						(*k1) = 2
					}
					//
					if (*info) == 0 {
						(*trfcon) = false
						//
						//                    Check residual of computed solution.
						//
						Dlacpy(func() *[]byte {y :=[]byte("Full"); return &y}(), n, (nrhs), (b), lda, (work), lda)
						Dgtt02(trans, n, (nrhs), (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (x), lda, (work), lda, &((*result)[1]))
						//
						//                    Check solution from generated exact solution.
						//
						Dget04(n, (nrhs), (x), lda, (xact), lda, rcondc, &((*result)[2]))
						//
						//                    Check the error bounds from iterative refinement.
						//
						Dgtt05(trans, n, (nrhs), (a), &((*(a))[(*m)+0]), &((*(a))[(*n)+(*m)+0]), (b), lda, (x), lda, (xact), lda, (rwork), &((*(rwork))[(*(nrhs))+0]), &((*result)[3]))
						(*nt) = 5
					}
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = (*k1); (*k) <= (*nt); (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Aladhd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" %s, fact='%c', trans='%c', N =%5d, type %2d, test %2d, ratio = %12.5f\n")
								return &y
							}(), *func() *[]byte {y :=[]byte("Dgtsvx"); return &y}(), (*fact), (*trans), (*n), (*imat), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label100:
					}
					//
					//                 Check the reciprocal of the condition number.
					//
					(*result)[5] = (*Dget06(rcond, rcondc))
					if (*result)[5] >= (*(thresh)) {
						if (*nfail) == 0 && (*nerrs) == 0 {
							Aladhd((nout), path)
						}
						WRITE((*(nout)), *func() *[]byte {
							y :=[]byte(" %s, fact='%c', trans='%c', N =%5d, type %2d, test %2d, ratio = %12.5f\n")
							return &y
						}(), *func() *[]byte {y :=[]byte("Dgtsvx"); return &y}(), (*fact), (*trans), (*n), (*imat), (*k), (*result)[(*k)-1])
						(*nfail) = (*nfail) + 1
					}
					(*nrun) = (*nrun) + (*nt) - (*k1) + 2
					//
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
	//     End of Ddrvgt
	//
}
