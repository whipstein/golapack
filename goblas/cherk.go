package goblas

import "math/cmplx"

// Cherk ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cherk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER K,LDA,LDC,N
//       CHARACTER TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX A(LDA,*),C(LDC,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Cherk  performs one of the hermitian rank k operations
//
//    C := alpha*A*A**H + beta*C,
//
// or
//
//    C := alpha*A**H*A + beta*C,
//
// where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
// matrix and  A  is an  n by k  matrix in the  first case and a  k by n
// matrix in the second case.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On  entry,   UPLO  specifies  whether  the  upper  or  lower
//           triangular  part  of the  array  C  is to be  referenced  as
//           follows:
//
//              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
//                                  is to be referenced.
//
//              UPLO = 'L' or 'L'   Only the  lower triangular part of  C
//                                  is to be referenced.
// \endverbatim
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER*1
//           On entry,  TRANS  specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.
//
//              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry,  N specifies the order of the matrix C.  N must be
//           at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
//           of  columns   of  the   matrix   A,   and  on   entry   with
//           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
//           matrix A.  K must be at least zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is REAL
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX array, dimension ( LDA, ka ), where ka is
//           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
//           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
//           part of the array  A  must contain the matrix  A,  otherwise
//           the leading  k by n  part of the array  A  must contain  the
//           matrix A.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
//           then  LDA must be at least  max( 1, n ), otherwise  LDA must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] BETA
// \verbatim
//          BETA is REAL
//           On entry, BETA specifies the scalar beta.
// \endverbatim
//
// \param[in,out] C
// \verbatim
//          C is COMPLEX array, dimension ( LDC, N )
//           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
//           upper triangular part of the array C must contain the upper
//           triangular part  of the  hermitian matrix  and the strictly
//           lower triangular part of C is not referenced.  On exit, the
//           upper triangular part of the array  C is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry  with  UPLO = 'L' or 'L',  the leading  n by n
//           lower triangular part of the array C must contain the lower
//           triangular part  of the  hermitian matrix  and the strictly
//           upper triangular part of C is not referenced.  On exit, the
//           lower triangular part of the array  C is overwritten by the
//           lower triangular part of the updated matrix.
//           Note that the imaginary parts of the diagonal elements need
//           not be set,  they are assumed to be zero,  and on exit they
//           are set to zero.
// \endverbatim
//
// \param[in] LDC
// \verbatim
//          LDC is INTEGER
//           On entry, LDC specifies the first dimension of C as declared
//           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
//
//  -- Modified 8-Nov-93 to set C(j,j) to REAL( C(j,j) ) when BETA = 1.
//     Ed Anderson, Cray Research Inc.
// \endverbatim
//
//  =====================================================================
func Cherk(major, uplo, trans *byte, n, k *int, alpha *float32, a *[][]complex64, lda *int, beta *float32, c *[][]complex64, ldc *int) {
	var temp complex64
	var rtemp float32
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
		Xerbla(func() *string { y := "Cherk"; return &y }(), &info)
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
					(*c)[j-1][j-1] = complex((*beta)*real(((*c)[j-1][j-1])), 0.0)
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
					(*c)[j-1][j-1] = complex((*beta)*real(((*c)[j-1][j-1])), 0.0)
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
		//        Form  C := alpha*A*A**H + beta*C.
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
					(*c)[j-1][j-1] = complex((*beta)*real(((*c)[j-1][j-1])), 0.0)
				} else {
					(*c)[j-1][j-1] = complex(real(((*c)[j-1][j-1])), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 {
						temp = complex(*alpha, 0.0) * complex64(cmplx.Conj(complex128((*a)[j-1][l-1])))
						for i = 1; i <= j-1; i++ {
							(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
						}
						(*c)[j-1][j-1] = complex(real(((*c)[j-1][j-1]))+real(temp*(*a)[i-1][l-1]), 0.0)
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
					(*c)[j-1][j-1] = complex((*beta)*real(((*c)[j-1][j-1])), 0.0)
					for i = j + 1; i <= *n; i++ {
						(*c)[i-1][j-1] = complex(*beta, 0.0) * (*c)[i-1][j-1]
					}
				} else {
					(*c)[j-1][j-1] = complex(real(((*c)[j-1][j-1])), 0.0)
				}
				for l = 1; l <= *k; l++ {
					if (*a)[j-1][l-1] != 0.0 {
						temp = complex(*alpha, 0.0) * complex64(cmplx.Conj(complex128((*a)[j-1][l-1])))
						(*c)[j-1][j-1] = complex(real(((*c)[j-1][j-1]))+real(temp*(*a)[j-1][l-1]), 0.0)
						for i = j + 1; i <= *n; i++ {
							(*c)[i-1][j-1] += temp * (*a)[i-1][l-1]
						}
					}
				}
			}
		}
	} else {
		//
		//        Form  C := alpha*A**H*A + beta*C.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= j-1; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp = temp + complex64(cmplx.Conj(complex128((*a)[l-1][i-1])))*(*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = complex((*alpha), 0.0) * temp
					} else {
						(*c)[i-1][j-1] = complex(*alpha, 0.0)*temp + complex(*beta, 0.0)*(*c)[i-1][j-1]
					}
				}
				rtemp = 0.0
				for l = 1; l <= *k; l++ {
					rtemp = rtemp + real(complex64(cmplx.Conj(complex128((*a)[l-1][j-1])))*(*a)[l-1][j-1])
				}
				if *beta == 0.0 {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp, 0.0)
				} else {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp+(*beta)*real(((*c)[j-1][j-1])), 0.0)
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				rtemp = 0.0
				for l = 1; l <= *k; l++ {
					rtemp = rtemp + real(complex64(cmplx.Conj(complex128((*a)[l-1][j-1])))*(*a)[l-1][j-1])
				}
				if *beta == 0.0 {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp, 0.0)
				} else {
					(*c)[j-1][j-1] = complex((*alpha)*rtemp+(*beta)*real(((*c)[j-1][j-1])), 0.0)
				}
				for i = j + 1; i <= *n; i++ {
					temp = 0.0
					for l = 1; l <= *k; l++ {
						temp = temp + complex64(cmplx.Conj(complex128((*a)[l-1][i-1])))*(*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = complex((*alpha), 0.0) * temp
					} else {
						(*c)[i-1][j-1] = complex(*alpha, 0.0)*temp + complex(*beta, 0.0)*(*c)[i-1][j-1]
					}
				}
			}
		}
	}
}
