package goblas

import 

// DchkorhrCol tests DorhrCol using DLAtsQR and Dgemqrt. Therefore, DLAtsQR
// (used in DGEQR) and Dgemqrt (used in DGEMQR) have to be tested
// before this test.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE DchkorhrCol( thresh, tsterr, nm, mval, nn, nval, nnb,
//                                nbval, nout)
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
// DchkorhrCol tests DorhrCol using DLAtsQR and Dgemqrt. Therefore, DLAtsQR
// (used in DGEQR) and Dgemqrt (used in DGEMQR) have to be tested
// before this test.
//
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
// \date November 2019
//
// \ingroup double_lin
//
//  =====================================================================
func DchkorhrCol(thresh *float64, tsterr *bool, nm *int, mval *[]int, nn *int, nval *[]int, nnb *int, nbval *[]int, nout *int) {
	ntests := new(int)
	_len := func() *[]byte {
		arr := make([]byte, -1)
		return &arr
	}()
	i := new(int)
	imb1 := new(int)
	inb1 := new(int)
	inb2 := new(int)
	j := new(int)
	t := new(int)
	m := new(int)
	n := new(int)
	mb1 := new(int)
	nb1 := new(int)
	nb2 := new(int)
	nfail := new(int)
	nerrs := new(int)
	nrun := new(int)
	result := func() *[]float64 {
		arr := make([]float64, 6)
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	_len := func() *[]byte {
		arr := make([]byte, -1)
		return &arr
	}()
	infot := new(int)
	nunit := new(int)
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.nunit = new(float64)
	common.infoc.infot = new(int)
	srnamt := new(int)
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
	//     .. Executable Statements ..
	//
	//     Initialize _constants
	//
	path ( 1 : 1) = 'D'
	path ( 2 : 3) = *func() *[]byte{y := []byte("HH"); return &y }()
	(*nrun) = 0
	(*nfail) = 0
	(*nerrs) = 0
	//
	//     Test the error exits
	//
	if *(tsterr) {
		DerrorhrCol(path, (nout))
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
			//           Only for M >= N
			//
			if Min((*m), (*n)) > 0 && (*m) >= (*n) {
				//
				//              Do for each possible value of mb1
				//
				for (*imB1) = 1; (*imB1) <= (*(nnb)); (*imB1)++ {
					(*mB1) = (*(nbval))[(*imB1)-1]
					//
					//                 Only for mb1 > N
					//
					if (*mB1) > (*n) {
						//
						//                    Do for each possible value of nb1
						//
						for (*inb1) = 1; (*inb1) <= (*(nnb)); (*inb1)++ {
							(*nb1) = (*(nbval))[(*inb1)-1]
							//
							//                       Do for each possible value of nb2
							//
							for (*inb2) = 1; (*inb2) <= (*(nnb)); (*inb2)++ {
								(*nb2) = (*(nbval))[(*inb2)-1]
								//
								if (*nb1) > 0 && (*nb2) > 0 {
									//
									//                             Test DorhrCol
									//
									DorhrCol01(m, n, mb1, nb1, nb2, result)
									//
									//                             Print information about the tests that did
									//                             not pass the threshold.
									//
									for  (*t) = 1;  (*t) <= (*ntests);  (*t)++ {
										if (*result)[(*t)-1] >= (*(thresh)) {
											if (*nfail) == 0 && (*nerrs) == 0 {
												Alahd((nout), path)
											}
											WRITE((*(nout)), *func() *[]byte{y := []byte("M=%5d, N=%5d, mb1=%5d, nb1=%5d, nb2=%5d test(%2d)=%12.5f\n"); return &y }(), (*m), (*n), (*mB1), (*nb1), (*nb2),  (*t), (*result)[(*t)-1])
											(*nfail) = (*nfail) + 1
										}
									}
									(*nrun) = (*nrun) + (*ntests)
								}
							}
						}
					}
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
	//     End of DchkorhrCol
	//
}
