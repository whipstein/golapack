package golapack

// Ssyrk ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER K,LDA,LDC,N
//       CHARACTER TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,*),C(LDC,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// SSYRK  performs one of the symmetric rank k operations
//
//    C := alpha*A*A**T + beta*C,
//
// or
//
//    C := alpha*A**T*A + beta*C,
//
// where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
// and  A  is an  n by k  matrix in the first case and a  k by n  matrix
// in the second case.
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
//              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
//                                  is to be referenced.
// \endverbatim
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER*1
//           On entry,  TRANS  specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
//
//              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
//
//              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
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
//           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
//           of rows of the matrix  A.  K must be at least zero.
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
//          A is REAL array, dimension ( LDA, ka ), where ka is
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
//          C is REAL array, dimension ( LDC, N )
//           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
//           upper triangular part of the array C must contain the upper
//           triangular part  of the  symmetric matrix  and the strictly
//           lower triangular part of C is not referenced.  On exit, the
//           upper triangular part of the array  C is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
//           lower triangular part of the array C must contain the lower
//           triangular part  of the  symmetric matrix  and the strictly
//           upper triangular part of C is not referenced.  On exit, the
//           lower triangular part of the array  C is overwritten by the
//           lower triangular part of the updated matrix.
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
// \ingroup single_blas_level3
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
func Ssyrk(uplo *byte, trans *byte, n *int, k *int, alpha *float32, a *[]float32, aoff, lda *int, beta *float32, c *[]float32, coff, ldc *int) {
	var upper bool
	var temp float32
	var i, info, j, l, nrowa int
	var one float32 = 1.0
	var zero float32 = 0.0

	//
	//     Test the input parameters.
	//
	if Lsame(trans, func() *byte { y := byte('N'); return &y }()) {
		nrowa = (*n)
	} else {
		nrowa = (*k)
	}
	upper = Lsame(uplo, func() *byte { y := byte('U'); return &y }())

	info = 0
	if (!upper) && (!Lsame(uplo, func() *byte { y := byte('L'); return &y }())) {
		info = 1
	} else if (!Lsame(trans, func() *byte { y := byte('N'); return &y }())) && (!Lsame(trans, func() *byte { y := byte('T'); return &y }())) && (!Lsame(trans, func() *byte { y := byte('C'); return &y }())) {
		info = 2
	} else if (*n) < 0 {
		info = 3
	} else if (*k) < 0 {
		info = 4
	} else if (*lda) < maxint(1, nrowa) {
		info = 7
	} else if (*ldc) < maxint(1, *n) {
		info = 10
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("SSYRK "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*n) == 0) || ((((*alpha) == zero) || ((*k) == 0)) && ((*beta) == one)) {
		return
	}
	//
	//     And when  alpha.eq.zero.
	//
	if (*alpha) == zero {
		if upper {
			if (*beta) == zero {
				for j = 1; j <= (*n); j++ {
					for i = 1; i <= j; i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				}
			} else {
				for j = 1; j <= (*n); j++ {
					for i = 1; i <= j; i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		} else {
			if (*beta) == zero {
				for j = 1; j <= (*n); j++ {
					for i = j; i <= (*n); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				}
			} else {
				for j = 1; j <= (*n); j++ {
					for i = j; i <= (*n); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
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
		//        Form  C := alpha*A*A**T + beta*C.
		//
		if upper {
			for j = 1; j <= (*n); j++ {
				if (*beta) == zero {
					for i = 1; i <= j; i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				} else if (*beta) != one {
					for i = 1; i <= j; i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
				for l = 1; l <= (*k); l++ {
					if (*a)[j-1+(l-1)*(*lda)+(*aoff)] != zero {
						temp = (*alpha) * (*a)[j-1+(l-1)*(*lda)+(*aoff)]
						for i = 1; i <= j; i++ {
							(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*c)[i-1+(j-1)*(*ldc)+(*coff)] + temp*(*a)[i-1+(l-1)*(*lda)+(*aoff)]
						}
					}
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				if (*beta) == zero {
					for i = j; i <= (*n); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				} else if (*beta) != one {
					for i = j; i <= (*n); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
				for l = 1; l <= (*k); l++ {
					if (*a)[j-1+(l-1)*(*lda)+(*aoff)] != zero {
						temp = (*alpha) * (*a)[j-1+(l-1)*(*lda)+(*aoff)]
						for i = j; i <= (*n); i++ {
							(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*c)[i-1+(j-1)*(*ldc)+(*coff)] + temp*(*a)[i-1+(l-1)*(*lda)+(*aoff)]
						}
					}
				}
			}
		}
	} else {
		//
		//        Form  C := alpha*A**T*A + beta*C.
		//
		if upper {
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= j; i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + (*a)[l-1+(i-1)*(*lda)+(*aoff)]*(*a)[l-1+(j-1)*(*lda)+(*aoff)]
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				for i = j; i <= (*n); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + (*a)[l-1+(i-1)*(*lda)+(*aoff)]*(*a)[l-1+(j-1)*(*lda)+(*aoff)]
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		}
	}
}
