package goblas

import (
	"math/cmplx"
)

// Cher2k ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha
//       REAL beta
//       INTEGER k,lda,ldb,ldc,n
//       CHARACTER trans,uplo
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
// Cher2k  performs one of the hermitian rank 2k operations
//
//    c := alpha*a*b**H + conjg( alpha )*b*a**H + beta*c,
//
// or
//
//    c := alpha*a**H*b + conjg( alpha )*b**H*a + beta*c,
//
// where  alpha and beta  are scalars with  beta  real,  c is an  n by n
// hermitian matrix and  a and b  are  n by k matrices in the first case
// and  k by n  matrices in the second case.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On  entry,   uplo  specifies  whether  the  upper  or  lower
//           triangular  part  of the  array  c  is to be  referenced  as
//           follows:
//
//              uplo = 'U' or 'u'   Only the  upper triangular part of  c
//                                  is to be referenced.
//
//              uplo = 'L' or 'l'   Only the  lower triangular part of  c
//                                  is to be referenced.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry,  trans  specifies the operation to be performed as
//           follows:
//
//              trans = 'N' or 'n'    c := alpha*a*b**H          +
//                                         conjg( alpha )*b*a**H +
//                                         beta*c.
//
//              trans = 'C' or 'c'    c := alpha*a**H*b          +
//                                         conjg( alpha )*b**H*a +
//                                         beta*c.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry,  n specifies the order of the matrix c.  n must be
//           at least zero.
// \endverbatim
//
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry with  trans = 'N' or 'n',  k  specifies  the number
//           of  columns  of the  matrices  a and b,  and on  entry  with
//           trans = 'C' or 'c',  k  specifies  the number of rows of the
//           matrices  a and b.  k must be at least zero.
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
//           k  when  trans = 'N' or 'n',  and is  n  otherwise.
//           Before entry with  trans = 'N' or 'n',  the  leading  n by k
//           part of the array  a  must contain the matrix  a,  otherwise
//           the leading  k by n  part of the array  a  must contain  the
//           matrix a.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in  the  calling  (sub)  program.   When  trans = 'N' or 'n'
//           then  lda must be at least  max( 1, n ), otherwise  lda must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] b
// \verbatim
//          b is COMPLEX array, dimension ( ldb, kb ), where kb is
//           k  when  trans = 'N' or 'n',  and is  n  otherwise.
//           Before entry with  trans = 'N' or 'n',  the  leading  n by k
//           part of the array  b  must contain the matrix  b,  otherwise
//           the leading  k by n  part of the array  b  must contain  the
//           matrix b.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is INTEGER
//           On entry, ldb specifies the first dimension of b as declared
//           in  the  calling  (sub)  program.   When  trans = 'N' or 'n'
//           then  ldb must be at least  max( 1, n ), otherwise  ldb must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is REAL
//           On entry, beta specifies the scalar beta.
// \endverbatim
//
// \param[in,out] c
// \verbatim
//          c is COMPLEX array, dimension ( ldc, n )
//           Before entry  with  uplo = 'U' or 'u',  the leading  n by n
//           upper triangular part of the array c must contain the upper
//           triangular part  of the  hermitian matrix  and the strictly
//           lower triangular part of c is not referenced.  On exit, the
//           upper triangular part of the array  c is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry  with  uplo = 'L' or 'l',  the leading  n by n
//           lower triangular part of the array c must contain the lower
//           triangular part  of the  hermitian matrix  and the strictly
//           upper triangular part of c is not referenced.  On exit, the
//           lower triangular part of the array  c is overwritten by the
//           lower triangular part of the updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set,  they are assumed to be zero,  and on exit they
//           are set to zero.
// \endverbatim
//
// \param[in] ldc
// \verbatim
//          ldc is INTEGER
//           On entry, ldc specifies the first dimension of c as declared
//           in  the  calling  (sub)  program.   ldc  must  be  at  least
//           max( 1, n ).
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
//
//  -- Modified 8-Nov-93 to set c(j,j) to REAL( c(j,j) ) when beta = 1.
//     Ed Anderson, Cray Research Inc.
// \endverbatim
//
//  =====================================================================
func Cher2k(major, uplo, trans *byte, n, k *int, alpha *complex64, a *[][]complex64, lda *int, b *[][]complex64, ldb *int, beta *float32, c *[][]complex64, ldc *int) {
	var temp1, temp2 complex64
	var i, info, j, l, nrowa int
	var upper bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
	//     December 2016
	//
	//     Test the input parameters.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		nrowa = *n
	} else {
		nrowa = *k
	}
	upper = Lsame(uplo, func() *byte { y := byte('U'); return &y }())

	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !upper && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 2
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('C'); return &y }()) {
		info = 3
	} else if *n < 0 {
		info = 4
	} else if *k < 0 {
		info = 5
	} else if *lda < max(1, nrowa) {
		info = 8
	} else if *ldb < max(1, nrowa) {
		info = 10
	} else if *ldc < max(1, *n) {
		info = 13
	}
	if info != 0 {
		name := "Cher2k"
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
	if *n == 0 || ((*alpha == 0.0 || *k == 0) && *beta == 1.0) {
		return
	}
	//
	//     And when  alpha.eq.zero.
	//
	if *alpha == 0.0 {
		if upper {
			if *beta == 0.0 {
				for j = 1; j <= *n; j++ {
					for i = 1; i <= j; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				}
			} else {
				for j = 1; j <= *n; j++ {
					for i = 1; i <= j-1; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
					(*c)[j-1][j-1] = complex(*beta, 0.0) * complex(real((*c)[j-1][j-1]), 0.0)
				}
			}
		} else {
			if *beta == 0.0 {
				for j = 1; j <= *n; j++ {
					for i = j; i <= *n; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				}
			} else {
				for j = 1; j <= *n; j++ {
					(*c)[j-1][j-1] = complex(*beta, 0.0) * complex(real((*c)[j-1][j-1]), 0.0)
					for i = j + 1; i <= *n; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
				}
			}
		}
		return
	}
	//
	//     Start the operations.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		//
		//        Form  c := alpha*a*b**H + conjg( alpha )*b*a**H +
		//                   c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = 1; i <= j; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = 1; i <= j-1; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
					(*c)[j-1][j-1] = complex(*beta, 0.0) * complex(real((*c)[j-1][j-1]), 0.0)
				} else {
					(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 || (*b)[j-1][l-1] != 0.0 {
						temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*b)[j-1][l-1])))
						temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*a)[j-1][l-1])))
						for i = 1; i <= j-1; i++ {
							(*c)[i-1][j-1] += (*a)[i-1][l-1]*temp1 + (*b)[i-1][l-1]*temp2
						}
						(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0) + complex(real((*a)[j-1][l-1]*temp1+(*b)[j-1][l-1]*temp2), 0.0)
					}
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = j; i <= *n; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = j + 1; i <= *n; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
					(*c)[j-1][j-1] = complex(*beta, 0.0) * complex(real((*c)[j-1][j-1]), 0.0)
				} else {
					(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 || (*b)[j-1][l-1] != 0.0 {
						temp1 = (*alpha) * complex64(cmplx.Conj(complex128((*b)[j-1][l-1])))
						temp2 = complex64(cmplx.Conj(complex128((*alpha) * (*a)[j-1][l-1])))
						for i = j + 1; i <= *n; i++ {
							(*c)[i-1][j-1] += (*a)[i-1][l-1]*temp1 + (*b)[i-1][l-1]*temp2
						}
						(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0) + complex(real((*a)[j-1][l-1]*temp1+(*b)[j-1][l-1]*temp2), 0.0)
					}
				}
			}
		}
	} else {
		//
		//        Form  c := alpha*a**H*b + conjg( alpha )*b**H*a +
		//                   c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= j; i++ {
					temp1 = 0.0
					temp2 = 0.0
					for l = 1; l <= *k; l++ {
						temp1 += complex64(cmplx.Conj(complex128((*a)[l-1][i-1]))) * (*b)[l-1][j-1]
						temp2 += complex64(cmplx.Conj(complex128((*b)[l-1][i-1]))) * (*a)[l-1][j-1]
					}
					if i == j {
						if *beta == 0.0 {
							(*c)[j-1][j-1] = complex(real((*alpha)*temp1+complex64(cmplx.Conj(complex128(*alpha)))*temp2), 0.0)
						} else {
							(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1]), 0.0) + complex(real((*alpha)*temp1+complex64(cmplx.Conj(complex128(*alpha)))*temp2), 0.0)
						}
					} else {
						if *beta == 0.0 {
							(*c)[i-1][j-1] = (*alpha)*temp1 + complex64(cmplx.Conj(complex128(*alpha)))*temp2
						} else {
							(*c)[i-1][j-1] = complex(*beta, 0.0)*(*c)[i-1][j-1] + (*alpha)*temp1 + complex64(cmplx.Conj(complex128(*alpha)))*temp2
						}
					}
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				for i = j; i <= *n; i++ {
					temp1 = 0.0
					temp2 = 0.0
					for l = 1; l <= *k; l++ {
						temp1 += complex64(cmplx.Conj(complex128((*a)[l-1][i-1]))) * (*b)[l-1][j-1]
						temp2 += complex64(cmplx.Conj(complex128((*b)[l-1][i-1]))) * (*a)[l-1][j-1]
					}
					if i == j {
						if *beta == 0.0 {
							(*c)[j-1][j-1] = complex(real((*alpha)*temp1+complex64(cmplx.Conj(complex128(*alpha)))*temp2), 0.0)
						} else {
							(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1])+real((*alpha)*temp1+complex64(cmplx.Conj(complex128(*alpha)))*temp2), 0.0)
						}
					} else {
						if *beta == 0.0 {
							(*c)[i-1][j-1] = (*alpha)*temp1 + complex64(cmplx.Conj(complex128(*alpha)))*temp2
						} else {
							(*c)[i-1][j-1] = complex(*beta, 0.0)*(*c)[i-1][j-1] + (*alpha)*temp1 + complex64(cmplx.Conj(complex128(*alpha)))*temp2
						}
					}
				}
			}
		}
	}
}
