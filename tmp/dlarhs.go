package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dlarhs chooses a set of nrhs random solution vectors and sets
// up the right hand sides for the linear system
//    op( A) * X = b,
// where op( A) may be A or A' (transpose of A).
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlarhs( path, xtype, uplo, trans, m, n, kl, ku, nrhs,
//                          a, lda, x, ldx, b, ldb, iseed, info)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans, uplo, xtype
//       CHARACTER*3        path
//       inTEGER            info, kl, ku, lda, ldb, ldx, m, n, nrhs
//       ..
//       .. Array Arguments ..
//       inTEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dlarhs chooses a set of nrhs random solution vectors and sets
// up the right hand sides for the linear system
//    op( A) * X = b,
// where op( A) may be A or A' (transpose of A).
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The type of the real matrix A.  path may be given in any
//          combination of upper and lower case.  Valid types include
//             xGE:  General m x n matrix
//             xGB:  General banded matrix
//             xPO:  Symmetric positive definite, 2-D storage
//             xPP:  Symmetric positive definite packed
//             xPB:  Symmetric positive definite banded
//             xSY:  Symmetric indefinite, 2-D storage
//             xSP:  Symmetric indefinite packed
//             xSB:  Symmetric indefinite banded
//             xTR:  Triangular
//             xTP:  Triangular packed
//             xTB:  Triangular banded
//             xQR:  General m x n matrix
//             xLQ:  General m x n matrix
//             xQL:  General m x n matrix
//             xRQ:  General m x n matrix
//          where the leading character indicates the precision.
// \endverbatim
//
// \param[in] xtype
// \verbatim
//          xtype is CHARACTER*1
//          Specifies how the exact solution X will be determined:
//          = 'N':  New solution; generate a random X.
//          = 'C':  Computed; use value of X on entry.
// \endverbatim
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the upper or lower triangular part of the
//          matrix A is stored, if A is symmetric.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the operation applied to the matrix A.
//          = 'N':  System is  A * x = b
//          = 'T':  System is  A'* x = b
//          = 'C':  System is  A'* x = b
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is inTEGER
//          The number or rows of the matrix A.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is inTEGER
//          The number of columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is inTEGER
//          Used only if A is a band matrix; specifies the number of
//          subdiagonals of A if A is a general band matrix or if A is
//          symmetric or triangular and uplo = 'L'; specifies the number
//          of superdiagonals of A if A is symmetric or triangular and
//          uplo = 'U'.  0 <= kl <= M-1.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is inTEGER
//          Used only if A is a general band matrix or if A is
//          triangular.
//
//          If path = xGB, specifies the number of superdiagonals of a,
//          and 0 <= ku <= N-1.
//
//          If path = xTR, xTP, or xTB, specifies whether or not the
//          matrix has unit diagonal:
//          = 1:  matrix has non-unit diagonal (default)
//          = 2:  matrix has unit diagonal
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is inTEGER
//          The number of right hand side vectors in the system A*X = B.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The test matrix whose type is given by path.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is inTEGER
//          The leading dimension of the array A.
//          If path = xGB, lda >= kl+ku+1.
//          If path = xPB, xSB, xHB, or xTB, lda >= kl+1.
//          Otherwise, lda >= max(1,M).
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is or output) DOUBLE PRECISION array, dimension(ldx,nrhs)
//          On entry, if xtype = 'C' (for 'Computed'), then X contains
//          the exact solution to the system of linear equations.
//          On exit, if xtype = 'N' (for 'New'), then X is initialized
//          with random values.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is inTEGER
//          The leading dimension of the array X.  If trans = 'N',
//          ldx >= max(1,N); if trans = 'T', ldx >= max(1,M).
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          The right hand side vector(s) for the system of equations,
//          computed from B = op(a) * x, where op(a) is determined by
//          trans.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is inTEGER
//          The leading dimension of the array B.  If trans = 'N',
//          ldb >= max(1,M); if trans = 'T', ldb >= max(1,N).
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is inTEGER array, dimension (4)
//          The seed vector for the random number generator (used in
//          Dlatms).  Modified on exit.
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is inTEGER
//          = 0: successful exit
//          < 0: if info = -i, the i-th argument had an illegal value
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
func Dlarhs(path *[]byte, xtype *byte, uplo *byte, trans *byte, m *int, n *int, kl *int, ku *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, b *[][]float64, ldb *int, iseed *[]int, info *int) {
	one := new(float64)
	zero := new(float64)
	band := new(bool)
	gen := new(bool)
	notran := new(bool)
	qrs := new(bool)
	sym := new(bool)
	tran := new(bool)
	tri := new(bool)
	c1 := new(byte)
	diag := new(byte)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	j := new(int)
	mb := new(int)
	nx := new(int)
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
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input parameters.
	//
	(*(info)) = 0
	(*c1) = (*(path))[0]
	(*c2) = (*(path))[1]
	(*tran) = blas.Lsame((trans), func() *byte {y := byte('T'); return &y }()) || blas.Lsame((trans), func() *byte {y := byte('C'); return &y }())
	(*notran) = !(*tran)
	(*gen) = (*blas.Lsame(&((*(path))[1]), func() *byte {y := byte('G'); return &y }()))
	(*qrs) = blas.Lsame(&((*(path))[1]), func() *byte {y := byte('Q'); return &y }()) || blas.Lsame(&((*(path))[2]), func() *byte {y := byte('Q'); return &y }())
	(*sym) = blas.Lsame(&((*(path))[1]), func() *byte {y := byte('P'); return &y }()) || blas.Lsame(&((*(path))[1]), func() *byte {y := byte('S'); return &y }())
	(*tri) = (*blas.Lsame(&((*(path))[1]), func() *byte {y := byte('T'); return &y }()))
	(*band) = (*blas.Lsame(&((*(path))[2]), func() *byte {y := byte('B'); return &y }()))
	if !blas.Lsame(c1, func() *[]byte {y := []byte("Double precision"); return &y }()) {
		(*(info)) = -1
	} else if !(blas.Lsame((xtype), func() *byte {y := byte('N'); return &y }()) || blas.Lsame((xtype), func() *byte {y := byte('C'); return &y }())) {
		(*(info)) = -2
	} else if ((*sym) || (*tri)) && !(blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) || blas.Lsame((uplo), func() *byte {y := byte('L'); return &y }())) {
		(*(info)) = -3
	} else if ((*gen) || (*qrs)) && !((*tran) || blas.Lsame((trans), func() *byte {y := byte('N'); return &y }())) {
		(*(info)) = -4
	} else if (*(m)) < 0 {
		(*(info)) = -5
	} else if (*(n)) < 0 {
		(*(info)) = -6
	} else if (*band) && (*(kl)) < 0 {
		(*(info)) = -7
	} else if (*band) && (*(ku)) < 0 {
		(*(info)) = -8
	} else if (*(nrhs)) < 0 {
		(*(info)) = -9
	} else if (!(*band) && (*(lda)) < MAX(1, (*(m)))) || ((*band) && ((*sym) || (*tri)) && (*(lda)) < (*(kl))+1) || ((*band) && (*gen) && (*(lda)) < (*(kl))+(*(ku))+1) {
		(*(info)) = -11
	} else if ((*notran) && (*(ldx)) < MAX(1, (*(n)))) || ((*tran) && (*(ldx)) < MAX(1, (*(m)))) {
		(*(info)) = -13
	} else if ((*notran) && (*(ldb)) < MAX(1, (*(m)))) || ((*tran) && (*(ldb)) < MAX(1, (*(n)))) {
		(*(info)) = -15
	}
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("Dlarhs"); return &y }(), -(*(info)))
		return
	}
	//
	//     Initialize X to nrhs random vectors unless xtype = 'C'.
	//
	if *tran {
		(*nx) = (*(m))
		(*mb) = (*(n))
	} else {
		(*nx) = (*(n))
		(*mb) = (*(m))
	}
	if !blas.Lsame((xtype), func() *byte {y := byte('C'); return &y }()) {
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), (iseed), (n), &((*(x))[0][(*j)-1]))
			//Label10:
		}
	}
	//
	//     Multiply X by op( A) using an appropriate
	//     matrix multiply routine.
	//
	if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("GE"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QR"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("LQ"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QL"); return &y }()) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("RQ"); return &y }())) {
		//
		//        General matrix
		//
		Dgemm((trans), func() *byte {y := byte('N'); return &y }(), mb, (nrhs), nx, one, (a), (lda), (x), (ldx), zero, (b), (ldb))
		//
	} else if (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PO"); return &y }())) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SY"); return &y }())) {
		//
		//        Symmetric matrix, 2-D storage
		//
		Dsymm(func() *[]byte {y := []byte("Left"); return &y }(), (uplo), (n), (nrhs), one, (a), (lda), (x), (ldx), zero, (b), (ldb))
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("GB"); return &y }()) {
		//
		//        General matrix, band storage
		//
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			Dgbmv((trans), mb, nx, (kl), (ku), one, (a), (lda), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), zero, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label20:
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PB"); return &y }()) {
		//
		//        Symmetric matrix, band storage
		//
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			Dsbmv((uplo), (n), (kl), one, (a), (lda), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), zero, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label30:
		}
		//
	} else if (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PP"); return &y }())) || (*Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SP"); return &y }())) {
		//
		//        Symmetric matrix, packed storage
		//
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			Dspmv((uplo), (n), one, (a), &((*(x))[0][(*j)-1]), func() *int {y := 1; return &y }(), zero, &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label40:
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TR"); return &y }()) {
		//
		//        Triangular matrix.  Note that for triangular matrices,
		//           ku = 1 => non-unit triangular
		//           ku = 2 => unit triangular
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (nrhs), (x), (ldx), (b), (ldb))
		if (*(ku)) == 2 {
			(*diag) = 'U'
		} else {
			(*diag) = 'N'
		}
		Dtrmm(func() *[]byte {y := []byte("Left"); return &y }(), (uplo), (trans), diag, (n), (nrhs), one, (a), (lda), (b), (ldb))
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TP"); return &y }()) {
		//
		//        Triangular matrix, packed storage
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (nrhs), (x), (ldx), (b), (ldb))
		if (*(ku)) == 2 {
			(*diag) = 'U'
		} else {
			(*diag) = 'N'
		}
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			Dtpmv((uplo), (trans), diag, (n), (a), &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label50:
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TB"); return &y }()) {
		//
		//        Triangular matrix, banded storage
		//
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (n), (nrhs), (x), (ldx), (b), (ldb))
		if (*(ku)) == 2 {
			(*diag) = 'U'
		} else {
			(*diag) = 'N'
		}
		for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
			dtbmv((uplo), (trans), diag, (n), (kl), (a), (lda), &((*(b))[0][(*j)-1]), func() *int {y := 1; return &y }())
			//Label60:
		}
		//
	} else {
		//
		//        If path is none of the above, return with an error code.
		//
		(*(info)) = -1
		Xerbla(func() *[]byte {y := []byte("Dlarhs"); return &y }(), -(*(info)))
	}
	//
	return
	//
	//     End of Dlarhs
	//
}
