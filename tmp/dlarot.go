package goblas

//    Dlarot applies a (Givens) rotation to two adjacent rows or
//    columns, where one element of the first and/or last column/row
//    for use on matrices stored in some format other than GE, so
//    that elements of the matrix may be used or modified for which
//    no array element is provided.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlarot( lrows, lleft, lright, nl, c, s, a, lda, xleft,
//                          xright)
//
//       .. Scalar Arguments ..
//       LOGICAL            lleft, lright, lrows
//       intEGER            lda, nl
//       DOUBLE PRECISION   c, s, xleft, xright
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Dlarot applies a (Givens) rotation to two adjacent rows or
//    columns, where one element of the first and/or last column/row
//    for use on matrices stored in some format other than GE, so
//    that elements of the matrix may be used or modified for which
//    no array element is provided.
//
//    one example is a symmetric matrix in SB format (bandwidth=4), for
//    which uplo='L':  Two adjacent rows will have the format:
//
//    row j:     C> C> C> C> C> .  .  .  .
//    row j+1:      C> C> C> C> C> .  .  .  .
//
//    '*' indicates elements for which storage is provided,
//    '.' indicates elements for which no storage is provided, but
//    are not necessarily zero; their values are determined by
//    symmetry.  ' ' indicates elements which are necessarily zero,
//     and have no storage provided.
//
//    Those columns which have two '*'s can be handled by DROT.
//    Those columns which have no '*'s can be ignored, since as long
//    as the Givens rotations are carefully applied to preserve
//    symmetry, their values are determined.
//    Those columns which have one '*' have to be handled separately,
//    by using separate variables "p" and "q":
//
//    row j:     C> C> C> C> C> p  .  .  .
//    row j+1:   q  C> C> C> C> C> .  .  .  .
//
//    The element p would have to be set correctly, then that column
//    is rotated, setting p to its new value.  The next call to
//    Dlarot would rotate columns j and j+1, using p, and restore
//    symmetry.  The element q would start out being zero, and be
//    made non-zero by the rotation.  Later, rotations would presumably
//    be chosen to zero q out.
//
//    Typical Calling Sequences: rotating the i-th and (i+1)-st rows.
//    ------- ------- ---------
//
//      General dense matrix:
//
//              CALL Dlarot(.TRUE.,.FALSE.,.FALSE., n, c,S,
//                      a(i,1),lda, dummy, dummy)
//
//      General banded matrix in GB format:
//
//              j = MAX(1, i-kl)
//              nl = Min( n, i+ku+1) + 1-j
//              CALL Dlarot( .TRUE., i-kl.GE.1, i+ku.LT.N, nl, c,S,
//                      a(ku+i+1-j,j),lda-1, xleft, xright)
//
//             [note that i+1-j is just Min(i,kl+1)]
//
//      Symmetric banded matrix in SY format, bandwidth k,
//      lower triangle only:
//
//              j = MAX(1, i-K)
//              nl = Min( K+1, i) + 1
//              CALL Dlarot( .TRUE., i-K.GE.1, .TRUE., nl, c,S,
//                      a(i,j), lda, xleft, xright)
//
//      Same, but upper triangle only:
//
//              nl = Min( K+1, N-i) + 1
//              CALL Dlarot( .TRUE., .TRUE., i+K.LT.N, nl, c,S,
//                      a(i,i), lda, xleft, xright)
//
//      Symmetric banded matrix in SB format, bandwidth k,
//      lower triangle only:
//
//             [same as for SY, except:]
//                  . . . .
//                      a(i+1-j,j), lda-1, xleft, xright)
//
//             [note that i+1-j is just Min(i,K+1)]
//
//      Same, but upper triangle only:
//                   . . .
//                      a(K+1,i), lda-1, xleft, xright)
//
//      Rotating columns is just the transpose of rotating rows, except
//      for GB and SB: (rotating columns i and i+1)
//
//      GB:
//              j = MAX(1, i-ku)
//              nl = Min( n, i+kl+1) + 1-j
//              CALL Dlarot( .TRUE., i-ku.GE.1, i+kl.LT.N, nl, c,S,
//                      a(ku+j+1-i,i),lda-1, xtOP, XBOTtm)
//
//             [note that ku+j+1-i is just MAX(1,ku+2-i)]
//
//      SB: (upper triangle)
//
//                   . . . . . .
//                      a(K+j+1-i,i),lda-1, xtOP, XBOTtm)
//
//      SB: (lower triangle)
//
//                   . . . . . .
//                      a(1,i),lda-1, xtOP, XBOTtm)
// \endverbatim
//
//  Arguments:
//  ==========
//
// \verbatim
//  lrows  - LOGICAL
//           If .TRUE., then Dlarot will rotate two rows.  If .FALSE.,
//           then it will rotate two columns.
//           Not modified.
//
//  lleft  - LOGICAL
//           If .TRUE., then xleft will be used instead of the
//           corresponding element of A for the first element in the
//           second row (if lrows=.FALSE.) or column (if lrows=.TRUE.)
//           If .FALSE., then the corresponding element of A will be
//           used.
//           Not modified.
//
//  lright - LOGICAL
//           If .TRUE., then xright will be used instead of the
//           corresponding element of A for the last element in the
//           first row (if lrows=.FALSE.) or column (if lrows=.TRUE.) If
//           .FALSE., then the corresponding element of A will be used.
//           Not modified.
//
//  nl     - intEGER
//           The length of the rows (if lrows=.TRUE.) or columns (if
//           lrows=.FALSE.) to be rotated.  If xleft and/or xright are
//           used, the columns/rows they are in should be included in
//           nl, e.g., if lleft = lright = .TRUE., then nl must be at
//           least 2.  The number of rows/columns to be rotated
//           exclusive of those involving xleft and/or xright may
//           not be negative, i.e., nl minus how many of lleft and
//           lright are .TRUE. must be at least zero; if not, Xerbla
//           will be called.
//           Not modified.
//
//  c, S   - DOUBLE PRECISION
//           Specify the Givens rotation to be applied.  If lrows is
//           true, then the matrix ( c  s)
//                                 (-s  c)  is applied from the left;
//           if false, then the transpose thereof is applied from the
//           right.  For a Givens rotation, C**2 + S**2 should be 1,
//           but this is not checked.
//           Not modified.
//
//  A      - DOUBLE PRECISION array.
//           The array containing the rows/columns to be rotated.  The
//           first element of A should be the upper left element to
//           be rotated.
//           Read and modified.
//
//  lda    - intEGER
//           The "effective" leading dimension of A.  If A contains
//           a matrix stored in GE or SY format, then this is just
//           the leading dimension of A as dimensioned in the calling
//           routine.  If A contains a matrix stored in band (GB or SB)
//           format, then this should be *one less* than the leading
//           dimension used in the calling routine.  Thus, if
//           A were dimensioned a(lda,*) in Dlarot, then a(1,j) would
//           be the j-th element in the first of the two rows
//           to be rotated, and a(2,j) would be the j-th in the second,
//           regardless of how the array may be stored in the calling
//           routine. [A cannot, however, actually be dimensioned thus,
//           since for band format, the row number may exceed lda, which
//           is not legal FORtran.]
//           If lrows=.TRUE., then lda must be at least 1, otherwise
//           it must be at least nl minus the number of .TRUE. values
//           in xleft and xright.
//           Not modified.
//
//  xleft  - DOUBLE PRECISION
//           If lleft is .TRUE., then xleft will be used and modified
//           instead of a(2,1) (if lrows=.TRUE.) or a(1,2)
//           (if lrows=.FALSE.).
//           Read and modified.
//
//  xright - DOUBLE PRECISION
//           If lright is .TRUE., then xright will be used and modified
//           instead of a(1,nl) (if lrows=.TRUE.) or a(nl,1)
//           (if lrows=.FALSE.).
//           Read and modified.
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
// \ingroup double_matgen
//
//  =====================================================================
func Dlarot(lrows *bool, lleft *bool, lright *bool, nl *int, c *float64, s *float64, a *[]float64, lda *int, xleft *float64, xright *float64) {
	iinc := new(int)
	inEXt := new(int)
	ix := new(int)
	iy := new(int)
	iyt := new(int)
	nt := new(int)
	Xt := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	yt := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
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
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Set up indices, arrays for ends
	//
	if *(lrows) {
		(*iinc) = (*(lda))
		(*inext) = 1
	} else {
		(*iinc) = 1
		(*inext) = (*(lda))
	}
	//
	if *(lleft) {
		(*nt) = 1
		(*ix) = 1 + (*iinc)
		(*iy) = 2 + (*(lda))
		(*xt)[0] = (*(a))[0]
		(*yt)[0] = (*(xleft))
	} else {
		(*nt) = 0
		(*ix) = 1
		(*iy) = 1 + (*inext)
	}
	//
	if *(lright) {
		(*iyt) = 1 + (*inext) + ((*(nl))-1)*(*iinc)
		(*nt) = (*nt) + 1
		(*xt)[(*nt)-(1)] = (*(xright))
		(*yt)[(*nt)-(1)] = (*(a))[(*iyt)-(1)]
	}
	//
	//     Check for errors
	//
	if (*(nl)) < (*nt) {
		Xerbla(func() *[]byte {y := []byte("Dlarot"); return &y }(), func() *int {y := 4; return &y }())
		return
	}
	if (*(lda)) <= 0 || (!(*(lrows)) && (*(lda)) < (*(nl))-(*nt)) {
		Xerbla(func() *[]byte {y := []byte("Dlarot"); return &y }(), func() *int {y := 8; return &y }())
		return
	}
	//
	//     Rotate
	//
	Drot((*(nl))-(*nt), &((*(a))[(*ix)-(1)]), iinc, &((*(a))[(*iy)-(1)]), iinc, (c), (s))
	Drot(NT, xt, func() *int {y := 1; return &y }(), yt, func() *int {y := 1; return &y }(), (c), (s))
	//
	//     Stuff values back into xleft, xright, etc.
	//
	if *(lleft) {
		(*(a))[0] = (*xt)[0]
		(*(xleft)) = (*yt)[0]
	}
	//
	if *(lright) {
		(*(xright)) = (*xt)[(*nt)-(1)]
		(*(a))[(*iyt)-(1)] = (*yt)[(*nt)-(1)]
	}
	//
	return
	//
	//     End of Dlarot
	//
}
