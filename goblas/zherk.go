package goblas

import (
	"math/cmplx"
)

// Zherk ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Zherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION alpha,beta
//       INTEGER k,lda,ldc,n
//       CHARACTER trans,uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 a(lda,*),c(ldc,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Zherk  performs one of the hermitian rank k operations
//
//    c := alpha*a*a**H + beta*c,
//
// or
//
//    c := alpha*a**H*a + beta*c,
//
// where  alpha and beta  are  real scalars,  c is an  n by n  hermitian
// matrix and  a  is an  n by k  matrix in the  first case and a  k by n
// matrix in the second case.
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
//              uplo = 'L' or 'L'   Only the  lower triangular part of  c
//                                  is to be referenced.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry,  trans  specifies the operation to be performed as
//           follows:
//
//              trans = 'N' or 'N'   c := alpha*a*a**H + beta*c.
//
//              trans = 'C' or 'C'   c := alpha*a**H*a + beta*c.
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
//           On entry with  trans = 'N' or 'N',  k  specifies  the number
//           of  columns   of  the   matrix   a,   and  on   entry   with
//           trans = 'C' or 'C',  k  specifies  the number of rows of the
//           matrix a.  k must be at least zero.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is DOUBLE PRECISION .
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX*16 array, dimension ( lda, ka ), where ka is
//           k  when  trans = 'N' or 'N',  and is  n  otherwise.
//           Before entry with  trans = 'N' or 'N',  the  leading  n by k
//           part of the array  a  must contain the matrix  a,  otherwise
//           the leading  k by n  part of the array  a  must contain  the
//           matrix a.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in  the  calling  (sub)  program.   When  trans = 'N' or 'N'
//           then  lda must be at least  max( 1, n ), otherwise  lda must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is DOUBLE PRECISION.
//           On entry, beta specifies the scalar beta.
// \endverbatim
//
// \param[in,out] c
// \verbatim
//          c is COMPLEX*16 array, dimension ( ldc, n )
//           Before entry  with  uplo = 'U' or 'u',  the leading  n by n
//           upper triangular part of the array c must contain the upper
//           triangular part  of the  hermitian matrix  and the strictly
//           lower triangular part of c is not referenced.  On exit, the
//           upper triangular part of the array  c is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry  with  uplo = 'L' or 'L',  the leading  n by n
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
// \author NAG Ltd.
//
// \date December 2016
//
// \ingroup complex16_blas_level3
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
//
//  -- Modified 8-Nov-93 to set c(j,j) to DBLE( c(j,j) ) when beta = 1.
//     Ed Anderson, Cray Research Inc.
// \endverbatim
//
//  =====================================================================
func Zherk(major, uplo, trans *byte, n, k *int, alpha *float64, a *[][]complex128, lda *int, beta *float64, c *[][]complex128, ldc *int) {
	var temp complex128
	var rtemp float64
	var i, info, j, l, nrowa int
	var upper bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
	} else if *ldc < max(1, *n) {
		info = 11
	}
	if info != 0 {
		Xerbla(func() *string { y := "Zherk"; return &y }(), &info)
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
					(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1]), 0.0)
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
					(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1]), 0.0)
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
		//        Form  c := alpha*a*a**H + beta*c.
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
					(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1]), 0.0)
				} else {
					(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 {
						temp = complex(*alpha, 0.0) * cmplx.Conj((*a)[j-1][l-1])
						for i = 1; i <= j-1; i++ {
							(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
						}
						(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1])+real(temp*(*a)[i-1][l-1]), 0.0)
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
					(*c)[j-1][j-1] = complex((*beta)*real((*c)[j-1][j-1]), 0.0)
					for i = j + 1; i <= *n; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
				} else {
					(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1]), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 {
						temp = complex(*alpha, 0.0) * cmplx.Conj((*a)[j-1][l-1])
						(*c)[j-1][j-1] = complex(real((*c)[j-1][j-1])+real(temp*(*a)[j-1][l-1]), 0.0)
						for i = j + 1; i <= *n; i++ {
							(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
						}
					}
				}
			}
		}
	} else {
		//
		//        Form  c := alpha*a**H*a + beta*c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= j-1; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += cmplx.Conj((*a)[l-1][i-1]) * (*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = complex(*alpha, 0.0) * temp
					} else {
						(*c)[i-1][j-1] = complex(*alpha, 0.0)*temp + complex(*beta, 0.0)*(*c)[i-1][j-1]
					}
				}
				rtemp = 0.0
				for l = 1; l <= *k; l++ {
					rtemp += real(cmplx.Conj((*a)[l-1][j-1]) * (*a)[l-1][j-1])
				}
				if *beta == 0.0 {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp, 0.0)
				} else {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp+(*beta)*real((*c)[j-1][j-1]), 0.0)
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				rtemp = 0.0
				for l = 1; l <= *k; l++ {
					rtemp += real(cmplx.Conj((*a)[l-1][j-1]) * (*a)[l-1][j-1])
				}
				if *beta == 0.0 {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp, 0.0)
				} else {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp+(*beta)*real((*c)[j-1][j-1]), 0.0)
				}
				for i = j + 1; i <= *n; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp += cmplx.Conj((*a)[l-1][i-1]) * (*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = complex(*alpha, 0.0) * temp
					} else {
						(*c)[i-1][j-1] = complex(*alpha, 0.0)*temp + complex(*beta, 0.0)*(*c)[i-1][j-1]
					}
				}
			}
		}
	}
}
