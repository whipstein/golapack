package goblas

import (
	"math/cmplx"
)

// Chemm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Chemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha,beta
//       INTEGER lda,ldb,ldc,m,n
//       CHARACTER side,uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX a(lda,*),b(ldb,*),c(ldc,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Chemm  performs one of the matrix-matrix operations
//
//    c := alpha*a*b + beta*c,
//
// or
//
//    c := alpha*b*a + beta*c,
//
// where alpha and beta are scalars, a is an hermitian matrix and  b and
// c are m by n matrices.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] side
// \verbatim
//          side is CHARACTER*1
//           On entry,  side  specifies whether  the  hermitian matrix  a
//           appears on the  left or right  in the  operation as follows:
//
//              side = 'L' or 'l'   c := alpha*a*b + beta*c,
//
//              side = 'R' or 'r'   c := alpha*b*a + beta*c,
// \endverbatim
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On  entry,   uplo  specifies  whether  the  upper  or  lower
//           triangular  part  of  the  hermitian  matrix   a  is  to  be
//           referenced as follows:
//
//              uplo = 'U' or 'u'   Only the upper triangular part of the
//                                  hermitian matrix is to be referenced.
//
//              uplo = 'L' or 'l'   Only the lower triangular part of the
//                                  hermitian matrix is to be referenced.
// \endverbatim
//
// \param[in] m
// \verbatim
//          m is INTEGER
//           On entry,  m  specifies the number of rows of the matrix  c.
//           m  must be at least zero.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the number of columns of the matrix c.
//           n  must be at least zero.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is COMPLEX
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX array, dimension ( lda, ka ), where ka is
//           m  when  side = 'L' or 'l'  and is n  otherwise.
//           Before entry  with  side = 'L' or 'l',  the  m by m  part of
//           the array  a  must contain the  hermitian matrix,  such that
//           when  uplo = 'U' or 'u', the leading m by m upper triangular
//           part of the array  a  must contain the upper triangular part
//           of the  hermitian matrix and the  strictly  lower triangular
//           part of  a  is not referenced,  and when  uplo = 'L' or 'l',
//           the leading  m by m  lower triangular part  of the  array  a
//           must  contain  the  lower triangular part  of the  hermitian
//           matrix and the  strictly upper triangular part of  a  is not
//           referenced.
//           Before entry  with  side = 'R' or 'r',  the  n by n  part of
//           the array  a  must contain the  hermitian matrix,  such that
//           when  uplo = 'U' or 'u', the leading n by n upper triangular
//           part of the array  a  must contain the upper triangular part
//           of the  hermitian matrix and the  strictly  lower triangular
//           part of  a  is not referenced,  and when  uplo = 'L' or 'l',
//           the leading  n by n  lower triangular part  of the  array  a
//           must  contain  the  lower triangular part  of the  hermitian
//           matrix and the  strictly upper triangular part of  a  is not
//           referenced.
//           Note that the imaginary parts  of the diagonal elements need
//           not be set, they are assumed to be zero.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the  calling (sub) program. When  side = 'L' or 'l'  then
//           lda must be at least  max( 1, m ), otherwise  lda must be at
//           least max( 1, n ).
// \endverbatim
//
// \param[in] b
// \verbatim
//          b is COMPLEX array, dimension ( ldb, n )
//           Before entry, the leading  m by n part of the array  b  must
//           contain the matrix b.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is INTEGER
//           On entry, ldb specifies the first dimension of b as declared
//           in  the  calling  (sub)  program.   ldb  must  be  at  least
//           max( 1, m ).
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is COMPLEX
//           On entry,  beta  specifies the scalar  beta.  When  beta  is
//           supplied as zero then c need not be set on input.
// \endverbatim
//
// \param[in,out] c
// \verbatim
//          c is COMPLEX array, dimension ( ldc, n )
//           Before entry, the leading  m by n  part of the array  c must
//           contain the matrix  c,  except when  beta  is zero, in which
//           case c need not be set on entry.
//           On exit, the array  c  is overwritten by the  m by n updated
//           matrix.
// \endverbatim
//
// \param[in] ldc
// \verbatim
//          ldc is INTEGER
//           On entry, ldc specifies the first dimension of c as declared
//           in  the  calling  (sub)  program.   ldc  must  be  at  least
//           max( 1, m ).
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
// \ingroup complex_blas_level3
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  lvel 3 Blas routine.
//
//  -- Written on 8-February-1989.
//     Jack Dongarra, Argonne National lboratory.
//     Iain Duff, AERE Harwell.
//     Jeremy Du Croz, Numerical Algorithms Group ld.
//     Sven Hammarling, Numerical Algorithms Group ld.
// \endverbatim
//
//  =====================================================================
func Chemm(major, side, uplo *byte, m, n *int, alpha *complex64, a *[][]complex64, lda *int, b *[][]complex64, ldb *int, beta *complex64, c *[][]complex64, ldc *int) {
	var temp1, temp2 complex64
	var i, info, j, k, nrowa int
	var upper bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
	//     December 2016
	//
	//     Set nrowa as the number of rows of a.
	//
	if Lsame(side, func() *byte { y := byte('L'); return &y }()) {
		nrowa = *m
	} else {
		nrowa = *n
	}
	upper = Lsame(uplo, func() *byte { y := byte('U'); return &y }())
	//
	//     Test the input parameters.
	//
	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !Lsame(side, func() *byte { y := byte('L'); return &y }()) && !Lsame(side, func() *byte { y := byte('R'); return &y }()) {
		info = 2
	} else if !upper && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 3
	} else if *m < 0 {
		info = 4
	} else if *n < 0 {
		info = 5
	} else if *lda < max(1, nrowa) {
		info = 8
	} else if *ldb < max(1, *m) {
		info = 10
	} else if *ldc < max(1, *m) {
		info = 13
	}
	if info != 0 {
		Xerbla(func() *string { y := "Chemm"; return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *m == 0 || *n == 0 || (*alpha == 0.0 && *beta == 1.0) {
		return
	}
	//
	//     And when  alpha.eq.zero.
	//
	if *alpha == 0.0 {
		if *beta == 0.0 {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] = 0.0
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
				}
			}
		}
		return
	}
	//
	//     Start the operations.
	//
	if Lsame(side, func() *byte { y := byte('L'); return &y }()) {
		//
		//        Form  c := alpha*a*b + beta*c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp1 = (*alpha) * (*b)[i-1][j-1]
					temp2 = 0.0
					for k = 1; k <= i-1; k++ {
						(*c)[k-1][j-1] += temp1 * (*a)[k-1][i-1]
						temp2 += (*b)[k-1][j-1] * complex64(cmplx.Conj(complex128((*a)[k-1][i-1])))
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = temp1*complex(real((*a)[i-1][i-1]), 0.0) + (*alpha)*temp2
					} else {
						(*c)[i-1][j-1] = (*beta)*(*c)[i-1][j-1] + temp1*complex(real((*a)[i-1][i-1]), 0.0) + (*alpha)*temp2
					}
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				for i = *m; i >= 1; i-- {
					temp1 = (*alpha) * (*b)[i-1][j-1]
					temp2 = 0.0
					for k = i + 1; k <= *m; k++ {
						(*c)[k-1][j-1] += temp1 * (*a)[k-1][i-1]
						temp2 += (*b)[k-1][j-1] * complex64(cmplx.Conj(complex128((*a)[k-1][i-1])))
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = temp1*complex(real((*a)[i-1][i-1]), 0.0) + (*alpha)*temp2
					} else {
						(*c)[i-1][j-1] = (*beta)*(*c)[i-1][j-1] + temp1*complex(real((*a)[i-1][i-1]), 0.0) + (*alpha)*temp2
					}
				}
			}
		}
	} else {
		//
		//        Form  c := alpha*b*a + beta*c.
		//
		for j = 1; j <= *n; j++ {
			temp1 = (*alpha) * complex(real((*a)[j-1][j-1]), 0.0)
			if *beta == 0.0 {
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] = temp1 * (*b)[i-1][j-1]
				}
			} else {
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] = (*beta)*(*c)[i-1][j-1] + temp1*(*b)[i-1][j-1]
				}
			}
			for k = 1; k <= j-1; k++ {
				if upper {
					temp1 = (*alpha) * (*a)[k-1][j-1]
				} else {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*a)[j-1][k-1])))
				}
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] += temp1 * (*b)[i-1][k-1]
				}
			}
			for k = j + 1; k <= *n; k++ {
				if upper {
					temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*a)[j-1][k-1])))
				} else {
					temp1 = (*alpha) * (*a)[k-1][j-1]
				}
				for i = 1; i <= *m; i++ {
					(*c)[i-1][j-1] += temp1 * (*b)[i-1][k-1]
				}
			}
		}
	}
}
