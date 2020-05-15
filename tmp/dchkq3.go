package goblas

import 

// Dchkq3 tests Dgeqp3.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchkq3( dotype, nm, mval, nn, nval, nnb, nbval, nxval,
//                          thresh, a, copya, s, tau, work, iwork,
//                          nout)
//
//       .. Scalar Arguments ..
//       intEGER            nm, nn, nnb, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       intEGER            iwork(*), mval(*), nbval(*), nval(*),
//      $                   nxval(*)
//       DOUBLE PRECISION   a(*), copya(*), S(*),
//      $                   tau(*), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkq3 tests Dgeqp3.
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
// \param[in] nm
// \verbatim
//          nm is intEGER
//          The number of values of M contained in the vector mval.
// \endverbatim
//
// \param[in] mval
// \verbatim
//          mval is intEGER array, dimension (nm)
//          The values of the matrix row dimension M.
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
// \param[in] nnb
// \verbatim
//          nnb is intEGER
//          The number of values of nb and nx contained in the
//          vectors nbval and nxval.  The blocking parameters are used
//          in pairs (nb,nx).
// \endverbatim
//
// \param[in] nbval
// \verbatim
//          nbval is intEGER array, dimension (nnb)
//          The values of the blocksize nb.
// \endverbatim
//
// \param[in] nxval
// \verbatim
//          nxval is intEGER array, dimension (nnb)
//          The values of the crossover point nx.
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
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (mmax*nmax)
//          where mmax is the maximum value of M in mval and nmax is the
//          maximum value of N in nval.
// \endverbatim
//
// \param[out] copya
// \verbatim
//          copya is DOUBLE PRECISION array, dimension (mmax*nmax)
// \endverbatim
//
// \param[out] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension
//                      (min(mmax,nmax))
// \endverbatim
//
// \param[out] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (mmax)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension
//                      (mmax*nmax + 4*nmax + mmax)
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
func Dchkq3(dotype *[]bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nxval *[]int, thresh *float64, a *[]float64, copya *[]float64, s *[]float64, tau *[]float64, work *[]float64, iwork *[]int, nout *int) {
	ntypes := new(int)
	ntests := new(int)
	one := new(float64)
	zero := new(float64)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	ihigh := new(int)
	ilow := new(int)
	im := new(int)
	imode := new(int)
	in := new(int)
	inb := new(int)
	info := new(int)
	istep := new(int)
	k := new(int)
	lda := new(int)
	lw := new(int)
	lwork := new(int)
	m := new(int)
	mnmin := new(int)
	mode := new(int)
	n := new(int)
	nb := new(int)
	nerrs := new(int)
	nfail := new(int)
	nrun := new(int)
	nx := new(int)
	eps := new(float64)
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
	(*ntypes) = 6
	(*ntests) = 3
	(*one) = 1.0
	(*zero) = 0.0
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
	iounit = common.infoc.iounit
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
	(*path)[1] = *func() *[]byte {y :=[]byte("Q3"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	for (*i) = 1; (*i) <= 4; (*i)++ {
		(*iseed)[(*i)-1] = (*iseedy)[(*i)-1]
		//Label10:
	}
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*infot) = 0
	//
	for (*im) = 1; (*im) <= (*(nm)); (*im)++ {
		//
		//        Do for each value of M in mval.
		//
		(*m) = (*(mval))[(*im)-1]
		(*lda) = (MAX(1, (*m)))
		//
		for (*in) = 1; (*in) <= (*(nn)); (*in)++ {
			//
			//           Do for each value of N in nval.
			//
			(*n) = (*(nval))[(*in)-1]
			(*mnmin) = (Min((*m), (*n)))
			(*lwork) = (*MAX(func() *int {y := 1; return &y }(), (*m)*MAX((*m), (*n))+4*(*mnmin)+MAX((*m), (*n)), (*m)*(*n)+2*(*mnmin)+4*(*n)))
			//
			for (*imode) = 1; (*imode) <= (*ntypes); (*imode)++ {
				if !(*(dotype))[(*imode)-1] {
					goto Label70
				}
				//
				//              Do for each type of matrix
				//                 1:  zero matrix
				//                 2:  one small singular value
				//                 3:  geometric distribution of singular values
				//                 4:  first n/2 columns fixed
				//                 5:  last n/2 columns fixed
				//                 6:  every second column fixed
				//
				(*mode) = (*imode)
				if (*imode) > 3 {
					(*mode) = 1
				}
				//
				//              Generate test matrix of size m by n using
				//              singular value distribution indicated by `mode'.
				//
				for (*i) = 1; (*i) <= (*n); (*i)++ {
					(*(iwork))[(*i)-1] = 0
					//Label20:
				}
				if (*imode) == 1 {
					Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), m, n, zero, zero, (copya), lda)
					for (*i) = 1; (*i) <= (*mnmin); (*i)++ {
						(*(s))[(*i)-1] = (*zero)
						//Label30:
					}
				} else {
					Dlatms(m, n, func() *[]byte {y :=[]byte("Uniform"); return &y }(), iseed, func() *[]byte {y :=[]byte("Nonsymm"); return &y }(), (s), mode, (*one)/(*eps), one, m, n, func() *[]byte {y :=[]byte("No packing"); return &y }(), (copya), lda, (work), info)
					if (*imode) >= 4 {
						if (*imode) == 4 {
							(*ilow) = 1
							(*istep) = 1
							(*ihigh) = (MAX(1, (*n)/2))
						} else if (*imode) == 5 {
							(*ilow) = (MAX(1, (*n)/2))
							(*istep) = 1
							(*ihigh) = (*n)
						} else if (*imode) == 6 {
							(*ilow) = 1
							(*istep) = 2
							(*ihigh) = (*n)
						}
						for (*i) = (*ilow); (*i) <= (*ihigh); (*i) += (*istep) {
							(*(iwork))[(*i)-1] = 1
							//Label40:
						}
					}
					Dlaord(func() *[]byte {y :=[]byte("Decreasing"); return &y }(), mnmin, (s), func() *int {y := 1; return &y }())
				}
				//
				for (*inb) = 1; (*inb) <= (*(nnb)); (*inb)++ {
					//
					//                 Do for each pair of values (nb,nx) in nbval and nxval.
					//
					(*nb) = (*(nbval))[(*inb)-1]
					Xlaenv(func() *int {y := 1; return &y }(), nb)
					(*nx) = (*(nxval))[(*inb)-1]
					Xlaenv(func() *int {y := 3; return &y }(), nx)
					//
					//                 Get a working copy of copya into A and a copy of
					//                 vector iwork.
					//
					Dlacpy(func() *[]byte {y :=[]byte("All"); return &y }(), m, n, (copya), lda, (a), lda)
					icopy(n, &((*(iwork))[0]), func() *int {y := 1; return &y }(), &((*(iwork))[(*n)+0]), func() *int {y := 1; return &y }())
					//
					//                 Compute the QR factorization with pivoting of A
					//
					(*lw) = (MAX(1, 2*(*n)+(*nb)*((*n)+1)))
					//
					//                 Compute the QP3 factorization of A
					//
					(*srnamt) = *func() *[]byte {y :=[]byte("Dgeqp3"); return &y }()
					Dgeqp3(m, n, (a), lda, &((*(iwork))[(*n)+0]), (tau), (work), LW, info)
					//
					//                 Compute norm(svd(a) - svd(r))
					//
					(*result)[0] = (*Dqrt12(m, n, (a), lda, (s), (work), lwork))
					//
					//                 Compute norm( A*P - Q*R)
					//
					(*result)[1] = (*Dqpt01(m, n, mnmin, (copya), (a), lda, (tau), &((*(iwork))[(*n)+0]), (work), lwork))
					//
					//                 Compute Q'*Q
					//
					(*result)[2] = (*Dqrt11(m, mnmin, (a), lda, (tau), (work), lwork))
					//
					//                 Print information about the tests that did not pass
					//                 the threshold.
					//
					for (*k) = 1; (*k) <= (*ntests); (*k)++ {
						if (*result)[(*k)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {
								y :=[]byte(" %s M =%5d, N =%5d, nb =%4d, type %2d, test %2d, ratio =%12.5f\n")
								return &y
							}(), *func() *[]byte {y :=[]byte("Dgeqp3"); return &y }(), (*m), (*n), (*nb), (*imode), (*k), (*result)[(*k)-1])
							(*nfail) = (*nfail) + 1
						}
						//Label50:
					}
					(*nrun) = (*nrun) + (*ntests)
					//
					//Label60:
				}
			Label70:
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
	//
	//     End of Dchkq3
	//
}
