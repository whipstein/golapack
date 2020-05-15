package goblas

// Dqrt13 generates a full-rank matrix that may be scaled to have large
// or small norm.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dqrt13( scale, m, n, a, lda, norma, iseed)
//
//       .. Scalar Arguments ..
//       intEGER            lda, m, n, scale
//       DOUBLE PRECISION   norma
//       ..
//       .. Array Arguments ..
//       intEGER            iseed( 4)
//       DOUBLE PRECISION   a( lda, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dqrt13 generates a full-rank matrix that may be scaled to have large
// or small norm.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] scale
// \verbatim
//          scale is intEGER
//          scale = 1: normally scaled matrix
//          scale = 2: matrix scaled up
//          scale = 3: matrix scaled down
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
//          The number of columns of A.
// \endverbatim
//
// \param[out] A
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
// \param[out] norma
// \verbatim
//          norma is DOUBLE PRECISION
//          The one-norm of A.
// \endverbatim
//
// \param[in,out] iseed
// \verbatim
//          iseed is integer array, dimension (4)
//          Seed for random number generator
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
func Dqrt13(scale *int, m *int, n *int, a *[][]float64, lda *int, norma *float64, iseed *[]int) {
	one := new(float64)
	info := new(int)
	j := new(int)
	bignum := new(float64)
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
	if (*(m)) <= 0 || (*(n)) <= 0 {
		return
	}
	//
	//     benign matrix
	//
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		Dlarnv(func() *int {y := 2; return &y }(), (iseed), (m), &((*(a))[0][(*j)-(1)]))
		if (*j) <= (*(m)) {
			(*(a))[(*j)-(1)][(*j)-(1)] = (*(a))[(*j)-(1)][(*j)-(1)] + SIGN(Dasum((m), &((*(a))[0][(*j)-(1)]), func() *int {y := 1; return &y }()), &((*(a))[(*j)-(1)][(*j)-(1)]))
		}
		//Label10:
	}
	//
	//     scaled versions
	//
	if (*(scale)) != 1 {
		(*(norma)) = (*Dlange(func() *[]byte {y := []byte("Max"); return &y }(), (m), (n), (a), (lda), dummy))
		(*smlnum) = (*Dlamch(func() *[]byte {y := []byte("Safe minimum"); return &y }()))
		(*bignum) = (*one) / (*smlnum)
		Dlabad(smlnum, bignum)
		(*smlnum) = (*smlnum) / Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }())
		(*bignum) = (*one) / (*smlnum)
		//
		if (*(scale)) == 2 {
			//
			//           matrix scaled up
			//
			dlascl(func() *[]byte {y := []byte("General"); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), (norma), bignum, (m), (n), (a), (lda), info)
		} else if (*(scale)) == 3 {
			//
			//           matrix scaled down
			//
			dlascl(func() *[]byte {y := []byte("General"); return &y }(), func() *int {y := 0; return &y }(), func() *int {y := 0; return &y }(), (norma), smlnum, (m), (n), (a), (lda), info)
		}
	}
	//
	(*(norma)) = (*Dlange(func() *[]byte {y := []byte("one-norm"); return &y }(), (m), (n), (a), (lda), dummy))
	return
	//
	//     End of Dqrt13
	//
}
