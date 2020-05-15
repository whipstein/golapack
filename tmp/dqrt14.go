package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dqrt14 checks whether X is in the row space of A or A'.  It does so
// by scaling both X and A such that their norms are in the range
//[sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of[A,X]
// (if trans = 'T') or an LQ factorization of[A',X]' (if trans = 'N'),
// and returning the norm of the trailing triangle, scaled by
// MAX(m,N,nrhs)*eps.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dqrt14( trans, m, n, nrhs, a, lda, x,
//                        ldx, work, lwork)
//
//       .. Scalar Arguments ..
//       CHARACTER          trans
//       intEGER            lda, ldx, lwork, m, n, nrhs
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), work( lwork), X( ldx, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt14 checks whether X is in the row space of A or A'.  It does so
// by scaling both X and A such that their norms are in the range
//[sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of[A,X]
// (if trans = 'T') or an LQ factorization of[A',X]' (if trans = 'N'),
// and returning the norm of the trailing triangle, scaled by
// MAX(m,N,nrhs)*eps.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          = 'N':  No transpose, check for X in the row space of A
//          = 'T':  Transpose, check for X in the row space of A'.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix A.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix A.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of columns
//          of X.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The M-by-N matrix A.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.
// \endverbatim
//
// \param[in] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension (ldx,nrhs)
//          If trans = 'N', the N-by-nrhs matrix X.
//          IF trans = 'T', the M-by-nrhs matrix X.
// \endverbatim
//
// \param[in] ldx
// \verbatim
//          ldx is intEGER
//          The leading dimension of the array X.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array dimension (lwork)
// \endverbatim
//
// \param[in] lwork
// \verbatim
//          lwork is intEGER
//          length of workspace array required
//          If trans = 'N', lwork >= (M+nrhs)*(N+2);
//          if trans = 'T', lwork >= (N+nrhs)*(M+2).
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
func Dqrt14(trans *byte, m *int, n *int, nrhs *int, a *[][]float64, lda *int, x *[][]float64, ldx *int, work *[]float64, lwork *int) (dqrt14Return *float64) {
	dqrt14Return = new(float64)
	zero := new(float64)
	one := new(float64)
	tpsd := new(bool)
	i := new(int)
	info := new(int)
	j := new(int)
	ldwork := new(int)
	anrm := new(float64)
	err := new(float64)
	xnrm := new(float64)
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
	(*(dqrt14Return)) = (*zero)
	if blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		(*ldwork) = (*(m)) + (*(nrhs))
		(*tpsd) = false
		if (*(lwork)) < ((*(m))+(*(nrhs)))*((*(n))+2) {
			Xerbla(func() *[]byte {y := []byte("Dqrt14"); return &y }(), func() *int {y := 10; return &y }())
			return
		} else if (*(n)) <= 0 || (*(nrhs)) <= 0 {
			return
		}
	} else if blas.Lsame((trans), func() *byte {y := byte('T'); return &y }()) {
		(*ldwork) = (*(m))
		(*tpsd) = true
		if (*(lwork)) < ((*(n))+(*(nrhs)))*((*(m))+2) {
			Xerbla(func() *[]byte {y := []byte("Dqrt14"); return &y }(), func() *int {y := 10; return &y }())
			return
		} else if (*(m)) <= 0 || (*(nrhs)) <= 0 {
			return
		}
	} else {
		Xerbla(func() *[]byte {y := []byte("Dqrt14"); return &y }(), func() *int {y := 1; return &y }())
		return
	}
	//
	//     Copy and scale A
	//
	Dlacpy(func() *[]byte {y := []byte("All"); return &y }(), (m), (n), (a), (lda), (work), ldwork)
	(*anrm) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (m), (n), (work), ldwork, rwork))
	if (*anrm) != (*zero) {
		dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), anrm, one, (m), (n), (work), ldwork, info)
	}
	//
	//     Copy X or X' into the right place and scale it
	//
	if *tpsd {
		//
		//        Copy X into columns n+1:n+nrhs of work
		//
		Dlacpy(func() *[]byte {y := []byte("All"); return &y }(), (m), (nrhs), (x), (ldx), &((*(work))[(*(n))*(*ldwork)+0]), ldwork)
		(*xnrm) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (m), (nrhs), &((*(work))[(*(n))*(*ldwork)+0]), ldwork, rwork))
		if (*xnrm) != (*zero) {
			dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), xnrm, one, (m), (nrhs), &((*(work))[(*(n))*(*ldwork)+0]), ldwork, info)
		}
		(*anrm) = (*Dlange(func() *[]byte {y := []byte("one-norm"); return &y }(), (m), (*(n))+(*(nrhs)), (work), ldwork, rwork))
		//
		//        Compute QR factorization of X
		//
		Dgeqr2((m), (*(n))+(*(nrhs)), (work), ldwork, &((*(work))[(*ldwork)*((*(n))+(*(nrhs)))+0]), &((*(work))[(*ldwork)*((*(n))+(*(nrhs)))+Min((*(m)), (*(n))+(*(nrhs)))+0]), info)
		//
		//        Compute largest entry in upper triangle of
		//        work(n+1:m,n+1:n+nrhs)
		//
		(*err) = (*zero)
		for (*j) = (*(n)) + 1; (*j) <= (*(n))+(*(nrhs)); (*j)++ {
			for (*i) = (*(n)) + 1; (*i) <= (Min((*(m)), (*j))); (*i)++ {
				(*err) = (MAX((*err), ABS(((*(work))[(*i)+((*j)-1)*(*(m))-1]))))
				//Label10:
			}
			//Label20:
		}
		//
	} else {
		//
		//        Copy X' into rows m+1:m+nrhs of work
		//
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
				(*(work))[(*(m))+(*j)+((*i)-1)*(*ldwork)-1] = (*(x))[(*i)-1][(*j)-1]
				//Label30:
			}
			//Label40:
		}
		//
		(*xnrm) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (nrhs), (n), &((*(work))[(*(m))+0]), ldwork, rwork))
		if (*xnrm) != (*zero) {
			dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), xnrm, one, (nrhs), (n), &((*(work))[(*(m))+0]), ldwork, info)
		}
		//
		//        Compute LQ factorization of work
		//
		Dgelq2(ldwork, (n), (work), ldwork, &((*(work))[(*ldwork)*(*(n))+0]), &((*(work))[(*ldwork)*((*(n))+1)+0]), info)
		//
		//        Compute largest entry in lower triangle in
		//        work(m+1:m+nrhs,m+1:n)
		//
		(*err) = (*zero)
		for (*j) = (*(m)) + 1; (*j) <= (*(n)); (*j)++ {
			for (*i) = (*j); (*i) <= (*ldwork); (*i)++ {
				(*err) = (MAX((*err), ABS(((*(work))[(*i)+((*j)-1)*(*ldwork)-1]))))
				//Label50:
			}
			//Label60:
		}
		//
	}
	//
	(*(dqrt14Return)) = (*err) / (DBLE(MAX((m), (n), (nrhs))) * Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	//
	return
	//
	//     End of Dqrt14
	//
}
