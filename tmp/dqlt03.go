package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Dqlt03 tests Dormql, which computes Q*C, Q'*C, C*Q or C*Q'.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqlt03( m, n, k, af, c, CC, q, lda, tau, work, lwork,
//                          rwork, result)
//
//       .. Scalar Arguments ..
//       intEGER            k, lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   af( lda, *), c( lda, *), cc( lda, *),
//      $                   q( lda, *), result(*), rwork(*), tau(*),
//      $                   work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqlt03 tests Dormql, which computes Q*C, Q'*C, C*Q or C*Q'.
//
// Dqlt03 compares the results of a call to Dormql with the results of
// forming Q explicitly by a call to Dorgql and then performing matrix
// multiplication by a call to Dgemm.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The order of the orthogonal matrix Q.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of rows or columns of the matrix C; C is m-by-n if
//          Q is applied from the left, or n-by-m if Q is applied from
//          the right.  N >= 0.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is intEGER
//          The number of elementary reflectors whose product defines the
//          orthogonal matrix Q.  M >= K >= 0.
// \endverbatim
//
// \param[in] af
// \verbatim
//          af is DOUBLE PRECISION array, dimension (lda,N)
//          Details of the QL factorization of an m-by-n matrix, as
//          returned by Dgeqlf. See SGEQLF for further details.
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[out] CC
// \verbatim
//          CC is DOUBLE PRECISION array, dimension (lda,N)
// \endverbatim
//
// \param[out] Q
// \verbatim
//          Q is DOUBLE PRECISION array, dimension (lda,M)
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the arrays af, c, CC, and Q.
// \endverbatim
//
// \param[in] tau
// \verbatim
//          tau is DOUBLE PRECISION array, dimension (min(m,N))
//          The scalar factors of the elementary reflectors corresponding
//          to the QL factorization in af.
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
//          The length of work.  lwork must be at least m, and should be
//          M*nb, where nb is the blocksize for this environment.
// \endverbatim
//
// \param[out] rwork
// \verbatim
//          rwork is DOUBLE PRECISION array, dimension (m)
// \endverbatim
//
// \param[out] result
// \verbatim
//          result is DOUBLE PRECISION array, dimension (4)
//          The test ratios compare two techniques for multiplying a
//          random matrix C by an m-by-m orthogonal matrix Q.
//          result1 = norm( Q*C - Q*C)  / ( m * norm(c) * eps)
//          result(2) = norm( C*Q - C*Q)  / ( m * norm(c) * eps)
//          result(3) = norm( Q'*C - Q'*C)/ ( m * norm(c) * eps)
//          result(4) = norm( C*Q' - C*Q')/ ( m * norm(c) * eps)
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
func Dqlt03(m *int, n *int, k *int, af *[][]float64, c *[][]float64, cc *[][]float64, q *[][]float64, lda *int, tau *[]float64, work *[]float64, lwork *int, rwork *[]float64, result *[]float64) {
	zero := new(float64)
	one := new(float64)
	rogue := new(float64)
	side := new(byte)
	trans := new(byte)
	info := new(int)
	iside := new(int)
	itrans := new(int)
	j := new(int)
	mc := new(int)
	minmn := new(int)
	nc := new(int)
	cnorm := new(float64)
	eps := new(float64)
	resid := new(float64)
	iseed := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
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
	(*zero) = 0.0
	(*one) = 1.0
	(*rogue) = -1.0e+10
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Scalars in common ..
	//     ..
	//     .. common blocks ..
	srnamt = common.srnamc.srnamt
	//     ..
	//     .. Data statements ..
	(*iseed)[0], (*iseed)[1], (*iseed)[2], (*iseed)[3] = 1988, 1989, 1990, 1991
	//     ..
	//     .. Executable Statements ..
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	(*minmn) = (Min((*(m)), (*(n))))
	//
	//     Quick return if possible
	//
	if (*minmn) == 0 {
		(*(result))[0] = (*zero)
		(*(result))[1] = (*zero)
		(*(result))[2] = (*zero)
		(*(result))[3] = (*zero)
		return
	}
	//
	//     Copy the last k columns of the factorization to the array Q
	//
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), (m), (m), rogue, rogue, (q), (lda))
	if (*(k)) > 0 && (*(m)) > (*(k)) {
		Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), (*(m))-(*(k)), (k), &((*(af))[0][(*(n))-(*(k))+0]), (lda), &((*(q))[0][(*(m))-(*(k))+0]), (lda))
	}
	if (*(k)) > 1 {
		Dlacpy(func() *[]byte {y := []byte("Upper"); return &y }(), (*(k))-1, (*(k))-1, &((*(af))[(*(m))-(*(k))+0][(*(n))-(*(k))+1]), (lda), &((*(q))[(*(m))-(*(k))+0][(*(m))-(*(k))+1]), (lda))
	}
	//
	//     Generate the m-by-m matrix Q
	//
	(*srnamt) = *func() *[]byte {y := []byte("Dorgql"); return &y }()
	Dorgql((m), (m), (k), (q), (lda), &((*(tau))[(*minmn)-(*(k))+0]), (work), (lwork), info)
	//
	for (*iside) = 1; (*iside) <= 2; (*iside)++ {
		if (*iside) == 1 {
			(*side) = 'L'
			(*mc) = (*(m))
			(*nc) = (*(n))
		} else {
			(*side) = 'R'
			(*mc) = (*(n))
			(*nc) = (*(m))
		}
		//
		//        Generate mc by nc matrix C
		//
		for (*j) = 1; (*j) <= (*nc); (*j)++ {
			Dlarnv(func() *int {y := 2; return &y }(), iseed, mc, &((*(c))[0][(*j)-1]))
			//Label10:
		}
		(*cnorm) = (*Dlange(func() *byte {y := byte('1'); return &y }(), mc, nc, (c), (lda), (rwork)))
		if (*cnorm) == 0.0 {
			(*cnorm) = (*one)
		}
		//
		for (*itrans) = 1; (*itrans) <= 2; (*itrans)++ {
			if (*itrans) == 1 {
				(*trans) = 'N'
			} else {
				(*trans) = 'T'
			}
			//
			//           Copy C
			//
			Dlacpy(func() *[]byte {y := []byte("Full"); return &y }(), mc, nc, (c), (lda), (cc), (lda))
			//
			//           Apply Q or Q' to C
			//
			(*srnamt) = *func() *[]byte {y := []byte("Dormql"); return &y }()
			if (*(k)) > 0 {
				Dormql(side, trans, mc, nc, (k), &((*(af))[0][(*(n))-(*(k))+0]), (lda), &((*(tau))[(*minmn)-(*(k))+0]), (cc), (lda), (work), (lwork), info)
			}
			//
			//           Form explicit product and subtract
			//
			if blas.Lsame(side, func() *byte {y := byte('L'); return &y }()) {
				Dgemm(trans, func() *[]byte {y := []byte("No transpose"); return &y }(), mc, nc, mc, -(*one), (q), (lda), (c), (lda), one, (cc), (lda))
			} else {
				Dgemm(func() *[]byte {y := []byte("No transpose"); return &y }(), trans, mc, nc, nc, -(*one), (c), (lda), (q), (lda), one, (cc), (lda))
			}
			//
			//           Compute error in the difference
			//
			(*resid) = (*Dlange(func() *byte {y := byte('1'); return &y }(), mc, nc, (cc), (lda), (rwork)))
			(*(result))[((*iside)-1)*2+(*itrans)-1] = (*resid) / (DBLE(MAX(1, (*(m)))) * (*cnorm) * (*eps))
			//
			//Label20:
		}
		//Label30:
	}
	//
	return
	//
	//     End of Dqlt03
	//
}
