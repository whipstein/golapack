package goblas

// Xlaenv sets certain machine- and problem-dependent quantities
// which will later be retrieved by ILAENV.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Xlaenv( ispec, nvalue)
//
//       .. Scalar Arguments ..
//       intEGER            ispec, nvalue
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Xlaenv sets certain machine- and problem-dependent quantities
// which will later be retrieved by ILAENV.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] ispec
// \verbatim
//          ispec is intEGER
//          Specifies the parameter to be set in the common array iparms.
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
//               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
//          = 6: the crossover point for the SVD (when reducing an m by n
//               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
//               this value, a QR factorization is used first to reduce
//               the matrix to a triangular form)
//          = 7: the number of processors
//          = 8: another crossover point, for the multishift QR and QZ
//               methods for nonsymmetric eigenvalue problems.
//          = 9: maximum size of the subproblems at the bottom of the
//               computation tree in the divide-and-conquer algorithm
//               (used by xGELSD and xGESDD)
//          =10: ieee NaN arithmetic can be trusted not to trap
//          =11: infinity arithmetic can be trusted not to trap
// \endverbatim
//
// \param[in] nvalue
// \verbatim
//          nvalue is intEGER
//          The value of the parameter specified by ispec.
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
// \ingroup aux_lin
//
//  =====================================================================
func Xlaenv(ispec *int, nvalue *int) {
	iparms := func() *[]int {
		arr := make([]int, 100)
		return &arr
	}()
	common.claenv.iparms = new([]int)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Arrays in common ..
	//     ..
	//     .. common blocks ..
	iparms = common.claenv.iparms
	//     ..
	//     .. Save statement ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*(ispec)) >= 1 && (*(ispec)) <= 9 {
		(*iparms)[(*(ispec))-(1)] = (*(nvalue))
	}
	//
	return
	//
	//     End of Xlaenv
	//
}
