package goblas

import (
	"math/cmplx"
)

// Cgemm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha,beta
//       INTEGER k,lda,ldb,ldc,m,n
//       CHARACTER transa,transb
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
// Cgemm  performs one of the matrix-matrix operations
//
//    c := alpha*op( a )*op( b ) + beta*c,
//
// where  op( x ) is one of
//
//    op( x ) = x   or   op( x ) = x**T   or   op( x ) = x**H,
//
// alpha and beta are scalars, and a, b and c are matrices, with op( a )
// an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] transa
// \verbatim
//          transa is CHARACTER*1
//           On entry, transa specifies the form of op( a ) to be used in
//           the matrix multiplication as follows:
//
//              transa = 'N' or 'n',  op( a ) = a.
//
//              transa = 'T' or 't',  op( a ) = a**T.
//
//              transa = 'C' or 'c',  op( a ) = a**H.
// \endverbatim
//
// \param[in] transb
// \verbatim
//          transb is CHARACTER*1
//           On entry, transb specifies the form of op( b ) to be used in
//           the matrix multiplication as follows:
//
//              transb = 'N' or 'n',  op( b ) = b.
//
//              transb = 'T' or 't',  op( b ) = b**T.
//
//              transb = 'C' or 'c',  op( b ) = b**H.
// \endverbatim
//
// \param[in] m
// \verbatim
//          m is INTEGER
//           On entry,  m  specifies  the number  of rows  of the  matrix
//           op( a )  and of the  matrix  c.  m  must  be at least  zero.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry,  n  specifies the number  of columns of the matrix
//           op( b ) and the number of columns of the matrix c. n must be
//           at least zero.
// \endverbatim
//
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry,  k  specifies  the number of columns of the matrix
//           op( a ) and the number of rows of the matrix op( b ). k must
//           be at least  zero.
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
//           k  when  transa = 'N' or 'n',  and is  m  otherwise.
//           Before entry with  transa = 'N' or 'n',  the leading  m by k
//           part of the array  a  must contain the matrix  a,  otherwise
//           the leading  k by m  part of the array  a  must contain  the
//           matrix a.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program. When  transa = 'N' or 'n' then
//           lda must be at least  max( 1, m ), otherwise  lda must be at
//           least  max( 1, k ).
// \endverbatim
//
// \param[in] b
// \verbatim
//          b is COMPLEX array, dimension ( ldb, kb ), where kb is
//           n  when  transb = 'N' or 'n',  and is  k  otherwise.
//           Before entry with  transb = 'N' or 'n',  the leading  k by n
//           part of the array  b  must contain the matrix  b,  otherwise
//           the leading  n by k  part of the array  b  must contain  the
//           matrix b.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is INTEGER
//           On entry, ldb specifies the first dimension of b as declared
//           in the calling (sub) program. When  transb = 'N' or 'n' then
//           ldb must be at least  max( 1, k ), otherwise  ldb must be at
//           least  max( 1, n ).
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
//           On exit, the array  c  is overwritten by the  m by n  matrix
//           ( alpha*op( a )*op( b ) + beta*c ).
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
// \author NAG Ltd.
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
//  Level 3 Blas routine.
//
//  -- Written on 8-February-1989.
//     Jack Dongarra, Argonne National Laboratory.
//     Iain Duff, AERE Harwell.
//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
//     Sven Hammarling, Numerical Algorithms Group Ltd.
// \endverbatim
//
//  =====================================================================
func Cgemm(major, transa, transb *byte, m, n, k *int, alpha *complex64, a *[][]complex64, lda *int, b *[][]complex64, ldb *int, beta *complex64, c *[][]complex64, ldc *int) {
	var temp complex64
	var i, info, j, l, nrowa, nrowb int
	var conja, conjb, nota, notb bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     Set  nota  and  notb  as  true if  a  and  b  respectively are not
	//     conjugated or transposed, set  conja and conjb  as true if  a  and
	//     b  respectively are to be  transposed but  not conjugated  and set
	//     nrowa, ncola and  nrowb  as the number of rows and  columns  of  a
	//     and the number of rows of  b  respectively.
	//
	nota = Lsame(transa, func() *byte { y := byte('N'); return &y }())
	notb = Lsame(transb, func() *byte { y := byte('N'); return &y }())
	conja = Lsame(transa, func() *byte { y := byte('C'); return &y }())
	conjb = Lsame(transb, func() *byte { y := byte('C'); return &y }())
	if nota {
		nrowa = *m
	} else {
		nrowa = *k
	}
	if notb {
		nrowb = *k
	} else {
		nrowb = *n
	}
	//
	//     Test the input parameters.
	//
	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !nota && !conja && !Lsame(transa, func() *byte { y := byte('T'); return &y }()) {
		info = 2
	} else if !notb && !conjb && !Lsame(transb, func() *byte { y := byte('T'); return &y }()) {
		info = 3
	} else if *m < 0 {
		info = 4
	} else if *n < 0 {
		info = 5
	} else if *k < 0 {
		info = 6
	} else if *lda < max(1, nrowa) {
		info = 9
	} else if *ldb < max(1, nrowb) {
		info = 11
	} else if *ldc < max(1, *m) {
		info = 14
	}
	if info != 0 {
		Xerbla(func() *string { y := "Cgemm"; return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if *m == 0 || *n == 0 || ((*alpha == 0.0 || *k == 0) && *beta == 1.0) {
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
	if notb {
		if nota {
			//
			//           Form  c := alpha*a*b + beta*c.
			//
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
				}
				for l = 1; l <= *k; l++ {
					temp = (*alpha) * (*b)[l-1][j-1]
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
					}
				}
			}
		} else if conja {
			//
			//           Form  c := alpha*a**H*b + beta*c.
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += complex64(cmplx.Conj(complex128((*a)[l-1][i-1]))) * (*b)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		} else {
			//
			//           Form  c := alpha*a**T*b + beta*c
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += (*a)[l-1][i-1] * (*b)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		}
	} else if nota {
		if conjb {
			//
			//           Form  c := alpha*a*b**H + beta*c.
			//
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
				}
				for l = 1; l <= *k; l++ {
					temp = (*alpha) * complex64(cmplx.Conj(complex128((*b)[j-1][l-1])))
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = (*c)[i-1][j-1] + temp*(*a)[i-1][l-1]
					}
				}
			}
		} else {
			//
			//           Form  c := alpha*a*b**T + beta*c
			//
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
				}
				for l = 1; l <= *k; l++ {
					temp = (*alpha) * (*b)[j-1][l-1]
					for i = 1; i <= *m; i++ {
						(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
					}
				}
			}
		}
	} else if conja {
		if conjb {
			//
			//           Form  c := alpha*a**H*b**H + beta*c.
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += complex64(cmplx.Conj(complex128((*a)[l-1][i-1]))) * complex64(cmplx.Conj(complex128((*b)[j-1][l-1])))
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		} else {
			//
			//           Form  c := alpha*a**H*b**T + beta*c
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += complex64(cmplx.Conj(complex128((*a)[l-1][i-1]))) * (*b)[j-1][l-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		}
	} else {
		if conjb {
			//
			//           Form  c := alpha*a**T*b**H + beta*c
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += (*a)[l-1][i-1] * complex64(cmplx.Conj(complex128((*b)[j-1][l-1])))
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		} else {
			//
			//           Form  c := alpha*a**T*b**T + beta*c
			//
			for j = 1; j <= *n; j++ {
				for i = 1; i <= *m; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += (*a)[l-1][i-1] * (*b)[j-1][l-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha) * temp
					} else {
						(*c)[i-1][j-1] = (*alpha)*temp + (*beta)*(*c)[i-1][j-1]
					}
				}
			}
		}
	}
}
