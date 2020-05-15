package goblas

import 

// Ilaenv returns problem-dependent parameters for the local
// environment.  See ispec for a description of the parameters.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       intEGER FUNCTION Ilaenv( ispec, name, opts, n1, n2, n3,
//                        n4)
//
//       .. Scalar Arguments ..
//       CHARACTER*(*)    name, opts
//       intEGER            ispec, n1, n2, n3, n4
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ilaenv returns problem-dependent parameters for the local
// environment.  See ispec for a description of the parameters.
//
// In this version, the problem-dependent parameters are contained in
// the integer array iparms in the common block claenv and the value
// with index ispec is copied to Ilaenv.  This version of Ilaenv is
// to be used in conjunction with Xlaenv in TESTinG and TIMinG.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] ispec
// \verbatim
//          ispec is intEGER
//          Specifies the parameter to be returned as the value of
//          Ilaenv.
//          = 1: the optimal blocksize; if this value is 1, an unblocked
//               algorithm will give the best performance.
//          = 2: the minimum block size for which the block routine
//               should be used; if the usable block size is less than
//               this value, an unblocked routine should be used.
//          = 3: the crossover point (in a block routine, for N less
//               than this value, an unblocked routine should be used)
//          = 4: the number of shifts, used in the nonsymmetric
//               eigenvalue routines
//          = 5: the minimum column dimension for blocking to be used;
//               rectangular blocks must have dimension at least k by m,
//               where k is given by Ilaenv(2,...) and m by Ilaenv(5,...)
//          = 6: the crossover point for the SVD (when reducing an m by n
//               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
//               this value, a QR factorization is used first to reduce
//               the matrix to a triangular form.)
//          = 7: the number of processors
//          = 8: the crossover point for the multishift QR and QZ methods
//               for nonsymmetric eigenvalue problems.
//          = 9: maximum size of the subproblems at the bottom of the
//               computation tree in the divide-and-conquer algorithm
//          =10: ieee NaN arithmetic can be trusted not to trap
//          =11: infinity arithmetic can be trusted not to trap
//
//          Other specifications (up to 100) can be added later.
// \endverbatim
//
// \param[in] name
// \verbatim
//          name is CHARACTER*(*)
//          The name of the calling subroutine.
// \endverbatim
//
// \param[in] opts
// \verbatim
//          opts is CHARACTER*(*)
//          The character options to the subroutine name, concatenated
//          into a single character string.  For example, uplo = 'U',
//          trans = 'T', and diag = 'N' for a triangular routine would
//          be specified as opts = 'UTN'.
// \endverbatim
//
// \param[in] n1
// \verbatim
//          n1 is intEGER
// \endverbatim
//
// \param[in] n2
// \verbatim
//          n2 is intEGER
// \endverbatim
//
// \param[in] n3
// \verbatim
//          n3 is intEGER
// \endverbatim
//
// \param[in] n4
// \verbatim
//          n4 is intEGER
//
//          Problem dimensions for the subroutine name; these may not all
//          be required.
// \endverbatim
//
// \return Ilaenv
// \verbatim
//          Ilaenv is intEGER
//          >= 0: the value of the parameter specified by ispec
//          < 0:  if Ilaenv = -k, the k-th argument had an illegal value.
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
// \date November 2017
//
// \ingroup aux_lin
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  The following conventions have been used when calling Ilaenv from the
//  lapACK routines:
//  1)  opts is a concatenation of all of the character options to
//      subroutine name, in the same order that they appear in the
//      argument list for name, even if they are not used in determining
//      the value of the parameter specified by ispec.
//  2)  The problem dimensions n1, n2, n3, n4 are specified in the order
//      that they appear in the argument list for name.  n1 is used
//      first, n2 second, and so on, and unused problem dimensions are
//      passed a value of -1.
//  3)  The parameter value returned by Ilaenv is checked for validity in
//      the calling subroutine.  For example, Ilaenv is used to retrieve
//      the optimal blocksize for STRtri as follows:
//
//      nb = Ilaenv( 1, 'STRtri', uplo // diag, n, -1, -1, -1)
//      IF( nb.LE.1) nb = MAX( 1, N)
// \endverbatim
//
//  =====================================================================
func Ilaenv(ispec *int, name *[]byte, opts *[]byte, n1 *int, n2 *int, n3 *int, n4 *int) (ilaenvReturn *int) {
	iparms := func() *[]int {
		arr := make([]int, 100)
		return &arr
	}()
	common.claenv.iparms = new([]int)
	//
	//  -- lapACK test routine (version 3.8.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     November 2017
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Arrays in common ..
	//     ..
	//     .. common blocks ..
	iparms = common.claenv.iparms
	//     ..
	//     .. Save statement ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*(ispec)) >= 1 && (*(ispec)) <= 5 {
		//
		//        Return a value from the common block.
		//
		if (*(name))[1] == *func() *[]byte {y :=[]byte("GEQR "); return &y }() {
			if (*(n3)) == 2 {
				(*(ilaenvReturn)) = (*iparms)[1]
			} else {
				(*(ilaenvReturn)) = (*iparms)[0]
			}
		} else if (*(name))[1] == *func() *[]byte {y :=[]byte("GELQ "); return &y }() {
			if (*(n3)) == 2 {
				(*(ilaenvReturn)) = (*iparms)[1]
			} else {
				(*(ilaenvReturn)) = (*iparms)[0]
			}
		} else {
			(*(ilaenvReturn)) = (*iparms)[(*(ispec))-1]
		}
		//
	} else if (*(ispec)) == 6 {
		//
		//        Compute SVD crossover point.
		//
		(*(ilaenvReturn)) = (*int(real(Min((*(n1)), (*(n2)))) * 1.6))
		//
	} else if (*(ispec)) >= 7 && (*(ispec)) <= 9 {
		//
		//        Return a value from the common block.
		//
		(*(ilaenvReturn)) = (*iparms)[(*(ispec))-1]
		//
	} else if (*(ispec)) == 10 {
		//
		//        IEEE NaN arithmetic can be trusted not to trap
		//
		//C        Ilaenv = 0
		(*(ilaenvReturn)) = 1
		if (*(ilaenvReturn)) == 1 {
			(*(ilaenvReturn)) = (*ieeeck(func() *int {y := 1; return &y }(), func() *float64 {y := 0.0; return &y }(), func() *float64 {y := 1.0; return &y }()))
		}
		//
	} else if (*(ispec)) == 11 {
		//
		//        Infinity arithmetic can be trusted not to trap
		//
		//C        Ilaenv = 0
		(*(ilaenvReturn)) = 1
		if (*(ilaenvReturn)) == 1 {
			(*(ilaenvReturn)) = (*ieeeck(func() *int {y := 0; return &y }(), func() *float64 {y := 0.0; return &y }(), func() *float64 {y := 1.0; return &y }()))
		}
		//
	} else {
		//
		//        Invalid value for ispec
		//
		(*(ilaenvReturn)) = -1
	}
	//
	return
	//
	//     End of Ilaenv
	//
}

// Ilaenv2stage ...
func Ilaenv2stage(ispec *int, name *[]byte, opts *[]byte, n1 *int, n2 *int, n3 *int, n4 *int) (ilaenv2stageReturn *int) {
	ilaenv2stageReturn = new(int)
	iispec := new(int)
	iparms := func() *[]int {
		arr := make([]int, 100)
		return &arr
	}()
	common.claenv.iparms = new([]int)
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Local variables ..
	//     .. External Functions ..
	//     ..
	//     .. Arrays in common ..
	//     ..
	//     .. common blocks ..
	iparms = common.claenv.iparms
	//     ..
	//     .. Save statement ..
	//     ..
	//     .. Executable Statements ..
	//
	if ((*(ispec)) >= 1) && ((*(ispec)) <= 5) {
		//
		//     1 <= ispec <= 5: 2stage eigenvalues SVD routines.
		//
		if (*(ispec)) == 1 {
			(*(ilaenv2stageReturn)) = (*iparms)[0]
		} else {
			(*iispec) = 16 + (*(ispec))
			(*(ilaenv2stageReturn)) = (*Iparam2stage(iispec, (name), (opts), (n1), (n2), (n3), (n4)))
		}
		//
	} else {
		//
		//        Invalid value for ispec
		//
		(*(ilaenv2stageReturn)) = -1
	}
	//
	return
}
