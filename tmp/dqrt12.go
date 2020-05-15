package goblas

import 

// Dqrt12 computes the singular values `svlues' of the upper trapezoid
// of a(1:M,1:N) and returns the ratio
//
//      || s - svlues||/(||svlues||*eps*max(m,N))
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       DOUBLE PRECISION FUNCTION Dqrt12( m, n, a, lda, s, work, lwork)
//
//       .. Scalar Arguments ..
//       intEGER            lda, lwork, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), S(*), work( lwork)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt12 computes the singular values `svlues' of the upper trapezoid
// of a(1:M,1:N) and returns the ratio
//
//      || s - svlues||/(||svlues||*eps*max(m,N))
// \endverbatim
//
//  Arguments:
//  ==========
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
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The M-by-N matrix A. Only the upper trapezoid is referenced.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of the array A.
// \endverbatim
//
// \param[in] S
// \verbatim
//          S is DOUBLE PRECISION array, dimension (min(m,N))
//          The singular values of the matrix A.
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
//          The length of the array work. lwork >= max(M*n + 4*min(m,N) +
//          max(m,N), M*n+2*min( m, N)+4*n).
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
func Dqrt12(m *int, n *int, a *[][]float64, lda *int, s *[]float64, work *[]float64, lwork *int) (dqrt12Return *float64) {
	dqrt12Return = new(float64)
	zero := new(float64)
	one := new(float64)
	i := new(int)
	info := new(int)
	iscl := new(int)
	j := new(int)
	mn := new(int)
	anrm := new(float64)
	bignum := new(float64)
	nrmsvl := new(float64)
	smlnum := new(float64)
	dummy := func() *[]float64 {
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
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Executable Statements ..
	//
	(*(dqrt12Return)) = (*zero)
	//
	//     Test that enough workspace is supplied
	//
	if (*(lwork)) < (MAX((*(m))*(*(n))+4*Min((*(m)), (*(n)))+MAX((*(m)), (*(n))), (*(m))*(*(n))+2*Min((*(m)), (*(n)))+4*(*(n)))) {
		Xerbla(func() *[]byte {y :=[]byte("Dqrt12"); return &y }(), func() *int {y := 7; return &y }())
		return
	}
	//
	//     Quick return if possible
	//
	(*mn) = (Min((*(m)), (*(n))))
	if (*mn) <= (*zero) {
		return
	}
	//
	(*nrmsvl) = (*Dnrm2(mn, (s), func() *int {y := 1; return &y }()))
	//
	//     Copy upper triangle of A into work
	//
	Dlaset(func() *[]byte {y :=[]byte("Full"); return &y }(), (m), (n), zero, zero, (work), (m))
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		for (*i) = 1; (*i) <= (Min((*j), (*(m)))); (*i)++ {
			(*(work))[((*j)-1)*(*(m))+(*i)-(1)] = (*(a))[(*i)-(1)][(*j)-(1)]
			//Label10:
		}
		//Label20:
	}
	//
	//     Get machine parameters
	//
	(*smlnum) = Dlamch(func() *byte {y := byte('S'); return &y }()) / Dlamch(func() *byte {y := byte('P'); return &y }())
	(*bignum) = (*one) / (*smlnum)
	Dlabad(smlnum, bignum)
	//
	//     Scale work if max entry outside range[smlnum,bignum]
	//
	(*anrm) = (*Dlange(func() *byte {y := byte('M'); return &y }(), (m), (n), (work), (m), dummy))
	(*iscl) = 0
	if (*anrm) > (*zero) && (*anrm) < (*smlnum) {
		//
		//        Scale matrix norm up to smlnum
		//
		dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), anrm, smlnum, (m), (n), (work), (m), info)
		(*iscl) = 1
	} else if (*anrm) > (*bignum) {
		//
		//        Scale matrix norm down to bignum
		//
		dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), anrm, bignum, (m), (n), (work), (m), info)
		(*iscl) = 1
	}
	//
	if (*anrm) != (*zero) {
		//
		//        Compute SVD of work
		//
		Dgebd2((m), (n), (work), (m), &((*(work))[(*(m))*(*(n))+0]), &((*(work))[(*(m))*(*(n))+(*mn)+0]), &((*(work))[(*(m))*(*(n))+2*(*mn)+0]), &((*(work))[(*(m))*(*(n))+3*(*mn)+0]), &((*(work))[(*(m))*(*(n))+4*(*mn)+0]), info)
		Dbdsqr(func() *[]byte {y :=[]byte("Upper"); return &y }(), mn, func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), &((*(work))[(*(m))*(*(n))+0]), &((*(work))[(*(m))*(*(n))+(*mn)+0]), dummy, mn, dummy, func() *int {y := 1; return &y }(), dummy, mn, &((*(work))[(*(m))*(*(n))+2*(*mn)+0]), info)
		//
		if (*iscl) == 1 {
			if (*anrm) > (*bignum) {
				dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), bignum, anrm, mn, func() *int {y := 1; return &y }(), &((*(work))[(*(m))*(*(n))+0]), mn, info)
			}
			if (*anrm) < (*smlnum) {
				dlascl(func() *byte {y := byte('G'); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), smlnum, anrm, mn, func() *int {y := 1; return &y }(), &((*(work))[(*(m))*(*(n))+0]), mn, info)
			}
		}
		//
	} else {
		//
		for (*i) = 1; (*i) <= (*mn); (*i)++ {
			(*(work))[(*(m))*(*(n))+(*i)-(1)] = (*zero)
			//Label30:
		}
	}
	//
	//     Compare s and singular values of work
	//
	Daxpy(mn, -(*one), (s), func() *int {y := 1; return &y }(), &((*(work))[(*(m))*(*(n))+0]), func() *int {y := 1; return &y }())
	(*(dqrt12Return)) = Dasum(mn, &((*(work))[(*(m))*(*(n))+0]), func() *int {y := 1; return &y }()) / (Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()) * DBLE(MAX((*(m)), (*(n)))))
	if (*nrmsvl) != (*zero) {
		(*(dqrt12Return)) = (*(dqrt12Return)) / (*nrmsvl)
	}
	//
	return
	//
	//     End of Dqrt12
	//
}
