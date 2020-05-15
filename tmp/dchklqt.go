package goblas

import 

// Dchklqt tests DGELQT and DGEMLQT.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dchklqt( thresh, tsterr, nm, mval, nn, nval, nnb,
//                           nbval, nout)
//
//       .. Scalar Arguments ..
//       LOGICAL            tsterr
//       intEGER            nm, nn, nnb, nout
//       DOUBLE PRECISION   thresh
//       ..
//       .. Array Arguments ..
//       intEGER            mval(*), nbval(*), nval(*)
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchklqt tests DGELQT and DGEMLQT.
// \endverbatim
//
//  Arguments:
//  ==========
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
//          The number of values of nb contained in the vector nbval.
// \endverbatim
//
// \param[in] nbval
// \verbatim
//          nbval is intEGER array, dimension (nbval)
//          The values of the blocksize nb.
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
func Dchklqt(thresh *float64, tsterr *bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nout *int) {
	ntests := new(int)
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	i := new(int)
	j := new(int)
	k := new(int)
	t := new(int)
	m := new(int)
	n := new(int)
	nb := new(int)
	nfail := new(int)
	nerrs := new(int)
	nrun := new(int)
	minmn := new(int)
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
	(*ntests) = 6
	//     ..
	//     .. Local Scalars ..
	//
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
	//     .. Executable Statements ..
	//
	//     Initialize _constants
	//
	(*path)[0] = 'D'
	(*path)[1] = *func() *[]byte {y :=[]byte("TQ"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	//
	//     Test the error exits
	//
	if *(tsterr) {
		Derrlqt(path, (nout))
	}
	(*infot) = 0
	//
	//     Do for each value of M in mval.
	//
	for (*i) = 1; (*i) <= (*(nm)); (*i)++ {
		(*m) = (*(mval))[(*i)-1]
		//
		//        Do for each value of N in nval.
		//
		for (*j) = 1; (*j) <= (*(nn)); (*j)++ {
			(*n) = (*(nval))[(*j)-1]
			//
			//        Do for each possible value of nb
			//
			(*minmn) = (Min((*m), (*n)))
			for (*k) = 1; (*k) <= (*(nnb)); (*k)++ {
				(*nb) = (*(nbval))[(*k)-1]
				//
				//              Test DGELQT and DGEMLQT
				//
				if ((*nb) <= (*minmn)) && ((*nb) > 0) {
					Dlqt04(m, n, nb, result)
					//
					//                 Print information about the tests that did not
					//                 pass the threshold.
					//
					for (*t) = 1; (*t) <= (*ntests); (*t)++ {
						if (*result)[(*t)-1] >= (*(thresh)) {
							if (*nfail) == 0 && (*nerrs) == 0 {
								Alahd((nout), path)
							}
							WRITE((*(nout)), *func() *[]byte {y :=[]byte(" M=%5d, N=%5d, nb=%4d test(%2d)=%12.5f\n"); return &y }(), (*m), (*n), (*nb), (*t), (*result)[(*t)-1])
							(*nfail) = (*nfail) + 1
						}
					}
					(*nrun) = (*nrun) + (*ntests)
				}
			}
		}
	}
	//
	//     Print a summary of the results.
	//
	Alasum(path, (nout), nfail, nrun, nerrs)
	//
	return
	//
	//     End of Dchklqt
	//
}
