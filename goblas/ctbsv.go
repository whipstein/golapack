package goblas

import (
	"math/cmplx"
)

// Ctbsv ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ctbsv(uplo,trans,diag,n,k,a,lda,x,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,k,lda,n
//       CHARACTER diag,trans,uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX a(lda,*),x(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ctbsv  solves one of the systems of equations
//
//    a*x = b,   or   a**T*x = b,   or   a**H*x = b,
//
// where b and x are n element vectors and a is an n by n unit, or
// non-unit, upper or lower triangular band matrix, with ( k + 1 )
// diagonals.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              uplo = 'U' or 'u'   a is an upper triangular matrix.
//
//              uplo = 'L' or 'l'   a is a lower triangular matrix.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry, trans specifies the equations to be solved as
//           follows:
//
//              trans = 'N' or 'n'   a*x = b.
//
//              trans = 'T' or 't'   a**T*x = b.
//
//              trans = 'C' or 'c'   a**H*x = b.
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//           On entry, diag specifies whether or not a is unit
//           triangular as follows:
//
//              diag = 'U' or 'u'   a is assumed to be unit triangular.
//
//              diag = 'N' or 'n'   a is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the order of the matrix a.
//           n must be at least zero.
// \endverbatim
//
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry with uplo = 'U' or 'u', k specifies the number of
//           super-diagonals of the matrix a.
//           On entry with uplo = 'L' or 'l', k specifies the number of
//           sub-diagonals of the matrix a.
//           k must satisfy  0 .le. k.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX array, dimension ( lda, n )
//           Before entry with uplo = 'U' or 'u', the leading ( k + 1 )
//           by n part of the array a must contain the upper triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1 ) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array a is not referenced.
//           The following program segment will transfer an upper
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, n
//                    m = k + 1 - j
//                    DO 10, I = MAX( 1, j - k ), j
//                       a( m + I, j ) = matrix( I, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with uplo = 'L' or 'l', the leading ( k + 1 )
//           by n part of the array a must contain the lower triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array a is not referenced.
//           The following program segment will transfer a lower
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, j = 1, n
//                    m = 1 - j
//                    DO 10, I = j, MIN( n, j + k )
//                       a( m + I, j ) = matrix( I, j )
//              10    CONTINUE
//              20 CONTINUE
//
//           Note that when diag = 'U' or 'u' the elements of the array a
//           corresponding to the diagonal elements of the matrix are not
//           referenced, but are assumed to be unity.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. lda must be at least
//           ( k + 1 ).
// \endverbatim
//
// \param[in,out] x
// \verbatim
//          x is COMPLEX array, dimension at least
//           ( 1 + ( n - 1 )*abs( incx ) ).
//           Before entry, the incremented array x must contain the n
//           element right-hand side vector b. On exit, x is overwritten
//           with the solution vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           x. incx must not be zero.
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG ld.
//
// \date December 2016
//
// \ingroup complex_blas_level2
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  lvel 2 Blas routine.
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National lb.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National lbs.
// \endverbatim
//
//  =====================================================================
func Ctbsv(major, uplo, trans, diag *byte, n, k *int, a *[][]complex64, lda *int, x *[]complex64, incx *int) {
	var temp complex64
	var i, info, ix, j, jx, kplus1, kx, l int
	var noconj, nounit bool
	//
	//  -- Reference BLAS level2 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
	//     December 2016
	//
	//     Test the input parameters.
	//
	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !Lsame(uplo, func() *byte { y := byte('U'); return &y }()) && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 2
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 3
	} else if !Lsame(diag, func() *byte { y := byte('U'); return &y }()) && !Lsame(diag, func() *byte { y := byte('N'); return &y }()) {
		info = 4
	} else if *n < 0 {
		info = 5
	} else if *k < 0 {
		info = 6
	} else if *lda < (*k)+1 {
		info = 8
	} else if *incx == 0 {
		info = 10
	}
	if info != 0 {
		name := "Ctbsv"
		if common.infoc.test {
			xerblaTest(&name, &info)
			return
		}
		Xerbla(&name, &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *n == 0 {
		return
	}
	//
	noconj = Lsame(trans, func() *byte { y := byte('T'); return &y }())
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	//
	//     Set up the start point in x if the increment is not unity. This
	//     will be  ( n - 1 )*incx  too small for descending loops.
	//
	if *incx <= 0 {
		kx = 1 - ((*n)-1)*(*incx)
	} else if *incx != 1 {
		kx = 1
	}
	//
	//     Start the operations. In this version the elements of a are
	//     accessed by sequentially with one pass through a.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  x := inv( a )*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					if (*x)[j-1] != 0.0 {
						l = kplus1 - j
						if nounit {
							(*x)[j-1] = (*x)[j-1] / (*a)[kplus1-1][j-1]
						}
						temp = (*x)[j-1]
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							(*x)[i-1] -= temp * (*a)[l+i-1][j-1]
						}
					}
				}
			} else {
				kx += ((*n) - 1) * (*incx)
				jx = kx
				for j = *n; j >= 1; j-- {
					kx -= *incx
					if (*x)[jx-1] != 0.0 {
						ix = kx
						l = kplus1 - j
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] / (*a)[kplus1-1][j-1]
						}
						temp = (*x)[jx-1]
						for i = j - 1; i >= max(1, j-(*k)); i-- {
							(*x)[ix-1] -= temp * (*a)[l+i-1][j-1]
							ix -= *incx
						}
					}
					jx -= *incx
				}
			}
		} else {
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					if (*x)[j-1] != 0.0 {
						l = 1 - j
						if nounit {
							(*x)[j-1] = (*x)[j-1] / (*a)[0][j-1]
						}
						temp = (*x)[j-1]
						for i = j + 1; i <= min(*n, j+(*k)); i++ {
							(*x)[i-1] -= temp * (*a)[l+i-1][j-1]
						}
					}
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					kx += *incx
					if (*x)[jx-1] != 0.0 {
						ix = kx
						l = 1 - j
						if nounit {
							(*x)[jx-1] = (*x)[jx-1] / (*a)[0][j-1]
						}
						temp = (*x)[jx-1]
						for i = j + 1; i <= min(*n, j+(*k)); i++ {
							(*x)[ix-1] -= temp * (*a)[l+i-1][j-1]
							ix += *incx
						}
					}
					jx += *incx
				}
			}
		}
	} else {
		//
		//        Form  x := inv( a**T )*x  or  x := inv( a**H )*x.
		//
		if Lsame(uplo, func() *byte { y := byte('U'); return &y }()) {
			kplus1 = (*k) + 1
			if *incx == 1 {
				for j = 1; j <= *n; j++ {
					temp = (*x)[j-1]
					l = kplus1 - j
					if noconj {
						for i = max(1, j-(*k)); i <= j-1; i++ {
							temp -= (*a)[l+i-1][j-1] * (*x)[i-1]
						}
						if nounit {
							temp = temp / (*a)[kplus1-1][j-1]
						}
					} else {
						for i = max(1, j-(*k)); i <= j-1; i++ {
							temp -= complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
						}
						if nounit {
							temp = temp / complex64(cmplx.Conj(complex128((*a)[kplus1-1][j-1])))
						}
					}
					(*x)[j-1] = temp
				}
			} else {
				jx = kx
				for j = 1; j <= *n; j++ {
					temp = (*x)[jx-1]
					ix = kx
					l = kplus1 - j
					if noconj {
						for i = max(1, j-(*k)); i <= j-1; i++ {
							temp -= (*a)[l+i-1][j-1] * (*x)[ix-1]
							ix += *incx
						}
						if nounit {
							temp = temp / (*a)[kplus1-1][j-1]
						}
					} else {
						for i = max(1, j-(*k)); i <= j-1; i++ {
							temp -= complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
							ix += *incx
						}
						if nounit {
							temp = temp / complex64(cmplx.Conj(complex128((*a)[kplus1-1][j-1])))
						}
					}
					(*x)[jx-1] = temp
					jx += *incx
					if j > (*k) {
						kx += *incx
					}
				}
			}
		} else {
			if *incx == 1 {
				for j = *n; j >= 1; j-- {
					temp = (*x)[j-1]
					l = 1 - j
					if noconj {
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							temp -= (*a)[l+i-1][j-1] * (*x)[i-1]
						}
						if nounit {
							temp = temp / (*a)[0][j-1]
						}
					} else {
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							temp -= complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[i-1]
						}
						if nounit {
							temp = temp / complex64(cmplx.Conj(complex128((*a)[0][j-1])))
						}
					}
					(*x)[j-1] = temp
				}
			} else {
				kx += ((*n) - 1) * (*incx)
				jx = kx
				for j = *n; j >= 1; j-- {
					temp = (*x)[jx-1]
					ix = kx
					l = 1 - j
					if noconj {
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							temp -= (*a)[l+i-1][j-1] * (*x)[ix-1]
							ix -= *incx
						}
						if nounit {
							temp = temp / (*a)[0][j-1]
						}
					} else {
						for i = min(*n, j+(*k)); i >= j+1; i-- {
							temp -= complex64(cmplx.Conj(complex128((*a)[l+i-1][j-1]))) * (*x)[ix-1]
							ix -= *incx
						}
						if nounit {
							temp = temp / complex64(cmplx.Conj(complex128((*a)[0][j-1])))
						}
					}
					(*x)[jx-1] = temp
					jx -= *incx
					if ((*n) - j) >= *k {
						kx -= *incx
					}
				}
			}
		}
	}
}
