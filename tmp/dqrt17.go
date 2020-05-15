package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dqrt17 computes the ratio
//
//    || R'*op(a) ||/(||A||*alpha*max(m,N,nrhs)*eps)
//
// where R = op(a)*X - b, op(a) is A or A', and
//
//    alpha = ||B|| if iresid = 1 (zero-residual problem)
//    alpha = ||R|| if iresid = 2 (otherwise).
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dqrt17( trans, iresid, m, n, nrhs, a,
//                        lda, x, ldx, b, ldb, c, work, lwork)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       intEGER            iresid, lda, ldb, ldx, lwork, m, n, nrhs
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), c( ldb, *),
//      $                   work( lwork), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt17 computes the ratio
//
//    || R'*op(a) ||/(||A||*alpha*max(m,N,nrhs)*eps)
//
// where R = op(a)*X - b, op(a) is A or A', and
//
//    alpha = ||B|| if iresid = 1 (zero-residual problem)
//    alpha = ||R|| if iresid = 2 (otherwise).
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies whether or not the transpose of A is used.
//          = 'N':  No transpose, op(a) = A.
//          = 'T':  Transpose, op(a) = A'.
// \endverbatim
//
// \param[in] iresid
// \verbatim
//          iresid is intEGER
//          iresid = 1 indicates zero-residual problem.
//          iresid = 2 indicates non-zero residual.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix A.
//          If trans = 'N', the number of rows of the matrix B.
//          If trans = 'T', the number of rows of the matrix X.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix  A.
//          If trans = 'N', the number of rows of the matrix X.
//          If trans = 'T', the number of rows of the matrix B.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of columns of the matrices X and B.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The m-by-n matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A. lda >= M.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          If trans = 'N', the n-by-nrhs matrix X.
//          If trans = 'T', the m-by-nrhs matrix X.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is intEGER
//          The leading dimension of the array X.
//          If trans = 'N', ldx >= N.
//          If trans = 'T', ldx >= M.
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          If trans = 'N', the m-by-nrhs matrix B.
//          If trans = 'T', the n-by-nrhs matrix B.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B.
//          If trans = 'N', ldb >= M.
//          If trans = 'T', ldb >= N.
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION array, dimension (ldb,nrhs)
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (lwork)
// \endverbatim
//
// \param[in] lwork
// \verbatim
//          lwork is intEGER
//          The length of the array work.  lwork >= nrhs*(M+N).
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
func Dqrt17(trans *byte, iresid *int, m *int, n *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, b *[][]float64, ldb *int, c *[][]float64, work *[]float64, lwork *int) (dqrt17Return *float64) {
	dqrt17Return = new(float64)
	zero := new(float64)
	one := new(float64)
	info := new(int)
	iscl := new(int)
	ncols := new(int)
	nrows := new(int)
	bignum := new(float64)
	err := new(float64)
	norma := new(float64)
	normb := new(float64)
	normrs := new(float64)
	smlnum := new(float64)
	rwork := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
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
	(*zero) = 0.0
	(*one) = 1.0
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
	//     .. Executable Statements ..
	//
	(*(dqrt17Return)) = (*zero)
	//
	if blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		(*nROWS) = (*(m))
		(*ncols) = (*(n))
	} else if blas.Lsame((trans), func() *byte {y := byte('T'); return &y }()) {
		(*nROWS) = (*(n))
		(*ncols) = (*(m))
	} else {
		Xerbla(func() *[]byte {y := []byte("Dqrt17"); return &y }(), func() *int {y := 1; return &y }())
		return
	}
	//
	if (*(lwork)) < (*ncols)*(*(nrhs)) {
		Xerbla(func() *[]byte {y := []byte("Dqrt17"); return &y }(), func() *int {y := 13; return &y }())
		return
	}
	//
	if (*(m)) <= 0 || (*(n)) <= 0 || (*(nrhs)) <= 0 {
		return
	}
	//
	(*norma) = (*Dlange(func() *[]byte {y := []byte("one-norm"); return &y }(), (m), (n), (a), (lda), rwork))
	(*smlnum) = Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()) / Dlamch(func() *[]byte {y := []byte("Precision"); return &y }())
	(*bignum) = (*one) / (*smlnum)
	(*iscl) = 0
	//
	//     compute residual and scale it
	//
	Dlacpy(func() *[]byte {y := []byte("All"); return &y }(), nrows, (nrhs), (b), (ldb), (c), (ldb))
	Dgemm((trans), func() *[]byte {y := []byte("No transpose"); return &y }(), nrows, (nrhs), ncols, -(*one), (a), (lda), (x), (ldx), one, (c), (ldb))
	(*normrs) = (*Dlange(func() *[]byte {y := []byte("Max"); return &y }(), nrows, (nrhs), (c), (ldb), rwork))
	if (*normrs) > (*smlnum) {
		(*iscl) = 1
		dlascl(func() *[]byte {y := []byte("General"); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), normrs, one, nrows, (nrhs), (c), (ldb), info)
	}
	//
	//     compute R'*a
	//
	Dgemm(func() *[]byte {y := []byte("Transpose"); return &y }(), (trans), (nrhs), ncols, nrows, one, (c), (ldb), (a), (lda), zero, (work), (nrhs))
	//
	//     compute and properly scale error
	//
	(*err) = (*Dlange(func() *[]byte {y := []byte("one-norm"); return &y }(), (nrhs), ncols, (work), (nrhs), rwork))
	if (*norma) != (*zero) {
		(*err) = (*err) / (*norma)
	}
	//
	if (*iscl) == 1 {
		(*err) = (*err) * (*normrs)
	}
	//
	if (*(iresid)) == 1 {
		(*normb) = (*Dlange(func() *[]byte {y := []byte("one-norm"); return &y }(), nrows, (nrhs), (b), (ldb), rwork))
		if (*normb) != (*zero) {
			(*err) = (*err) / (*normb)
		}
	} else {
		if (*normrs) != (*zero) {
			(*err) = (*err) / (*normrs)
		}
	}
	//
	(*(dqrt17Return)) = (*err) / (Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()) * DBLE(MAX((m), (n), (nrhs))))
	return
	//
	//     End of Dqrt17
	//
}
