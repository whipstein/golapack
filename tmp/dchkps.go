package goblas

import 

// Dchkps tests DPSTRF.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkps( dotype, nn, nval, nnb, nbval, nrank, rankval,
//                          thresh, tsterr, nmax, a, afac, perm, piv, work,
//                          rwork, nout)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   thresh
//       intEGER            nmax, nn, nnb, nout, nrank
//       LOGICAL            tsterr
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a(*), afac(*), perm(*), rwork(*),
//      $                   work(*)
//       intEGER            nbval(*), nval(*), piv(*), rankval(*)
//       LOGICAL            dotype(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkps tests DPSTRF.
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
//          The values of the block size nb.
// \endverbatim
//
// \param[in] nrank
// \verbatim
//          nrank is intEGER
//          The number of values of rank contained in the vector rankval.
// \endverbatim
//
// \param[in] rankval
// \verbatim
//          rankval is intEGER array, dimension (nbval)
//          The values of the block size nb.
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
// \param[out] perm
// \verbatim
//          perm is DOUBLE PRECISION array, dimension (nmax*nmax)
// \endverbatim
//
// \param[out] piv
// \verbatim
//          piv is intEGER array, dimension (nmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (nmax*3)
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (nmax)
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
func Dchkps(dotype *[]bool, nn *int, nval *[]int, nnb *int, nbval *[]int, nrank *int, rankval *[]int, thresh *float64, tsterr *bool, nmax *int, a *[]float64, afac *[]float64, perm *[]float64, piv *[]int, work *[]float64, rwork *[]float64, nout *int) {
	one := new(float64)
	ntypes := new(int)
	anorm := new(float64)
	cndnum := new(float64)
	result := new(float64)
	tol := new(float64)
	comprank := new(int)
	i := new(int)
	imat := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	irank := new(int)
	iuplo := new(int)
	izero := new(int)
	kl := new(int)
	ku := new(int)
	lda := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nimat := new(int)
	nrun := new(int)
	rank := new(int)
	rankdiff := new(int)
	dist := new(byte)
	_type := new(byte)
	uplo := new(byte)
	path := func() *[]byte {
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
	uplos := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	infot := new(int)
	nunit := new(int)
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
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
	(*ntypes) = 9
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
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
	(*uplos)[0], (*uplos)[1] = 'U', 'L'
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize _constants and the random number seed.
	//
	(*path)[0] = *func() *[]byte {y :=[]byte("Double precision"); return &y }()
	(*path)[1] = *func() *[]byte {y :=[]byte("PS"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-(1)] = (*iseedy)[(*i)-(1)]
		//Label100:
	}
	//
	//     Test the error exits
	//
	if *(tsterr) {
		Derrps(path, (nout))
	}
	(*infot) = 0
	Xlaenv(func() *int {y := 2; return &y }(), func() *int {y := 2; return &y }())
	//
	//     Do for each value of N in nval
	//
	for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
		(*n) = (*(nval))[(*in)-(1)]
		(*lda) = (MAX((*n), 1))
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
			if !(*(dotype))[(*imat)-(1)] {
				goto Label140
			}
			//
			//              Do for each value of rank in rankval
			//
			for (*irank) = 1; (*irank) <= (*(nrank)); (*irank)++ {
				//
				//              Only repeat test 3 to 5 for different ranks
				//              Other tests use full rank
				//
				if ((*imat) < 3 || (*imat) > 5) && (*irank) > 1 {
					goto Label130
				}
				//
				(*rank) = (*CEILING(((*n) * DBLE(((*(rankval))[(*irank)-(1)]))) / 100.e+0))
				//
				//
				//           Do first for uplo = 'U', then for uplo = 'L'
				//
				for (*iuplo) = 1; (*iuplo) <= 2; (*iuplo)++ {
					(*uplo) = (*uplos)[(*iuplo)-(1)]
					//
					//              Set up parameters with Dlatb5 and generate a test matrix
					//              with DLAtmT.
					//
					dlatb5(path, imat, n, _type, kl, ku, anorm, mode, cndnum, dist)
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("DLAtmT"); return &y }()
					dlatmt(n, n, dist, iseed, _type, (rwork), mode, cndnum, anorm, rank, kl, ku, uplo, (a), lda, (work), info)
					//
					//              Check error code from DLAtmT.
					//
					if (*info) != 0 {
						Alaerh(path, func() *[]byte {y :=[]byte("DLAtmT"); return &y }(), info, func() *int {y := 0; return &y }(), uplo, n, n, -1, -1, -1, imat, nfail, nerrs, (nout))
						goto Label120
					}
					//
					//              Do for each value of nb in nbval
					//
					for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
						(*nb) = (*(nbval))[(*inb)-(1)]
						Xlaenv(func() *int {y := 1; return &y }(), nb)
						//
						//                 Compute the pivoted L*l' or U'*U factorization
						//                 of the matrix.
						//
						Dlacpy(uplo, n, n, (a), lda, (afac), lda)
						(*srnamt) = *func() *[]byte {y :=[]byte("DPSTRF"); return &y }()
						//
						//                 Use default tolerance
						//
						(*tol) = -(*one)
						Dpstrf(uplo, n, (afac), lda, (piv), comprank, tol, (work), info)
						//
						//                 Check error code from DPSTRF.
						//
						if ((*info) < (*izero)) || ((*info) != (*izero) && (*rank) == (*n)) || ((*info) <= (*izero) && (*rank) < (*n)) {
							Alaerh(path, func() *[]byte {y :=[]byte("DPSTRF"); return &y }(), info, izero, uplo, n, n, -1, -1, nb, imat, nfail, nerrs, (nout))
							goto Label110
						}
						//
						//                 Skip the test if info is not 0.
						//
						if (*info) != 0 {
							goto Label110
						}
						//
						//                 Re_construct matrix from factors and compute residual.
						//
						//                 perm holds permuted L*l^T or U^T*U
						//
						Dpst01(uplo, n, (a), lda, (afac), lda, (perm), lda, (piv), (rwork), result, comprank)
						//
						//                 Print information about the tests that did not pass
						//                 the threshold or where computed rank was not rank.
						//
						if (*n) == 0 {
							(*comprank) = 0
						}
						(*rankdiff) = (*rank) - (*comprank)
						if (*result) >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" uplo = '%c', N =%5d, rank =%3d, Diff =%5d, nb =%4d, type %2d, ratio =%12.5f\n")
								return &y
							}(), (*uplo), (*n), (*rank), (*rankdiff), (*nb), (*imat), (*result))
							(*nfail) = (*nfail) + 1
						}
						(*nrun) = (*nrun) + 1
					Label110:
					}
					//
				Label120:
				}
			Label130:
			}
		Label140:
		}
		//Label150:
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchkps
	//
}
