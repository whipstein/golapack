package golapack

// Zgemm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 ALPHA,BETA
//       INTEGER K,LDA,LDB,LDC,M,N
//       CHARACTER TRANSA,TRANSB
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// ZGEMM  performs one of the matrix-matrix operations
//
//    C := alpha*op( A )*op( B ) + beta*C,
//
// where  op( X ) is one of
//
//    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
//
// alpha and beta are scalars, and A, B and C are matrices, with op( A )
// an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] TRANSA
// \verbatim
//          TRANSA is CHARACTER*1
//           On entry, TRANSA specifies the form of op( A ) to be used in
//           the matrix multiplication as follows:
//
//              TRANSA = 'N' or 'n',  op( A ) = A.
//
//              TRANSA = 'T' or 't',  op( A ) = A**T.
//
//              TRANSA = 'C' or 'c',  op( A ) = A**H.
// \endverbatim
//
// \param[in] TRANSB
// \verbatim
//          TRANSB is CHARACTER*1
//           On entry, TRANSB specifies the form of op( B ) to be used in
//           the matrix multiplication as follows:
//
//              TRANSB = 'N' or 'n',  op( B ) = B.
//
//              TRANSB = 'T' or 't',  op( B ) = B**T.
//
//              TRANSB = 'C' or 'c',  op( B ) = B**H.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is INTEGER
//           On entry,  M  specifies  the number  of rows  of the  matrix
//           op( A )  and of the  matrix  C.  M  must  be at least  zero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry,  N  specifies the number  of columns of the matrix
//           op( B ) and the number of columns of the matrix C. N must be
//           at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry,  K  specifies  the number of columns of the matrix
//           op( A ) and the number of rows of the matrix op( B ). K must
//           be at least  zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is COMPLEX*16
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
//           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
//           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
//           part of the array  A  must contain the matrix  A,  otherwise
//           the leading  k by m  part of the array  A  must contain  the
//           matrix A.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
//           LDA must be at least  max( 1, m ), otherwise  LDA must be at
//           least  max( 1, k ).
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
//           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
//           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
//           part of the array  B  must contain the matrix  B,  otherwise
//           the leading  n by k  part of the array  B  must contain  the
//           matrix B.
// \endverbatim
//
// \param[in] LDB
// \verbatim
//          LDB is INTEGER
//           On entry, LDB specifies the first dimension of B as declared
//           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
//           LDB must be at least  max( 1, k ), otherwise  LDB must be at
//           least  max( 1, n ).
// \endverbatim
//
// \param[in] BETA
// \verbatim
//          BETA is COMPLEX*16
//           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
//           supplied as zero then C need not be set on input.
// \endverbatim
//
// \param[in,out] C
// \verbatim
//          C is COMPLEX*16 array, dimension ( LDC, N )
//           Before entry, the leading  m by n  part of the array  C must
//           contain the matrix  C,  except when  beta  is zero, in which
//           case C need not be set on entry.
//           On exit, the array  C  is overwritten by the  m by n  matrix
//           ( alpha*op( A )*op( B ) + beta*C ).
// \endverbatim
//
// \param[in] LDC
// \verbatim
//          LDC is INTEGER
//           On entry, LDC specifies the first dimension of C as declared
//           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
// \endverbatim
//
//  =====================================================================
func Zgemm(transa *byte, transb *byte, m *int, n *int, k *int, alpha *complex128, a *[]complex128, aoff, lda *int, b *[]complex128, boff, ldb *int, beta *complex128, c *[]complex128, coff, ldc *int) {
	var conja, conjb, nota, notb bool
	var temp complex128
	var i, info, j, l, nrowa, nrowb int
	var one complex128 = (1.0 + 0.0*1i)
	var zero complex128 = (0.0 + 0.0*1i)

	//
	//     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
	//     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
	//     B  respectively are to be  transposed but  not conjugated  and set
	//     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
	//     and the number of rows of  B  respectively.
	//
	nota = Lsame(transa, func() *byte { y := byte('N'); return &y }())
	notb = Lsame(transb, func() *byte { y := byte('N'); return &y }())
	conja = Lsame(transa, func() *byte { y := byte('C'); return &y }())
	conjb = Lsame(transb, func() *byte { y := byte('C'); return &y }())
	if nota {
		nrowa = (*m)
	} else {
		nrowa = (*k)
	}
	if notb {
		nrowb = (*k)
	} else {
		nrowb = (*n)
	}
	//
	//     Test the input parameters.
	//
	info = 0
	if (!nota) && (!conja) && (!Lsame(transa, func() *byte { y := byte('T'); return &y }())) {
		info = 1
	} else if (!notb) && (!conjb) && (!Lsame(transb, func() *byte { y := byte('T'); return &y }())) {
		info = 2
	} else if (*m) < 0 {
		info = 3
	} else if (*n) < 0 {
		info = 4
	} else if (*k) < 0 {
		info = 5
	} else if (*lda) < maxint(1, nrowa) {
		info = 8
	} else if (*ldb) < maxint(1, nrowb) {
		info = 10
	} else if (*ldc) < maxint(1, *m) {
		info = 13
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("ZGEMM "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if ((*m) == 0) || ((*n) == 0) || ((((*alpha) == zero) || ((*k) == 0)) && ((*beta) == one)) {
		return
	}
	//
	//     And when  alpha.eq.zero.
	//
	if (*alpha) == zero {
		if (*beta) == zero {
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
				}
			}
		} else {
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
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
			//           Form  C := alpha*A*B + beta*C.
			//
			for j = 1; j <= (*n); j++ {
				if (*beta) == zero {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				} else if (*beta) != one {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
				for l = 1; l <= (*k); l++ {
					temp = (*alpha) * (*b)[l-1+(j-1)*(*ldb)+(*boff)]
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*c)[i-1+(j-1)*(*ldc)+(*coff)] + temp*(*a)[i-1+(l-1)*(*lda)+(*aoff)]
					}
				}
			}
		} else if conja {
			//
			//           Form  C := alpha*A**H*B + beta*C.
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + conjc128((*a)[l-1+(i-1)*(*lda)+(*aoff)])*(*b)[l-1+(j-1)*(*ldb)+(*boff)]
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		} else {
			//
			//           Form  C := alpha*A**T*B + beta*C
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + (*a)[l-1+(i-1)*(*lda)+(*aoff)]*(*b)[l-1+(j-1)*(*ldb)+(*boff)]
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		}
	} else if nota {
		if conjb {
			//
			//           Form  C := alpha*A*B**H + beta*C.
			//
			for j = 1; j <= (*n); j++ {
				if (*beta) == zero {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				} else if (*beta) != one {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
				for l = 1; l <= (*k); l++ {
					temp = (*alpha) * conjc128((*b)[j-1+(l-1)*(*ldb)+(*boff)])
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*c)[i-1+(j-1)*(*ldc)+(*coff)] + temp*(*a)[i-1+(l-1)*(*lda)+(*aoff)]
					}
				}
			}
		} else {
			//
			//           Form  C := alpha*A*B**T + beta*C
			//
			for j = 1; j <= (*n); j++ {
				if (*beta) == zero {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = zero
					}
				} else if (*beta) != one {
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*beta) * (*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
				for l = 1; l <= (*k); l++ {
					temp = (*alpha) * (*b)[j-1+(l-1)*(*ldb)+(*boff)]
					for i = 1; i <= (*m); i++ {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*c)[i-1+(j-1)*(*ldc)+(*coff)] + temp*(*a)[i-1+(l-1)*(*lda)+(*aoff)]
					}
				}
			}
		}
	} else if conja {
		if conjb {
			//
			//           Form  C := alpha*A**H*B**H + beta*C.
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + conjc128((*a)[l-1+(i-1)*(*lda)+(*aoff)])*conjc128((*b)[j-1+(l-1)*(*ldb)+(*boff)])
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		} else {
			//
			//           Form  C := alpha*A**H*B**T + beta*C
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + conjc128((*a)[l-1+(i-1)*(*lda)+(*aoff)])*(*b)[j-1+(l-1)*(*ldb)+(*boff)]
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		}
	} else {
		if conjb {
			//
			//           Form  C := alpha*A**T*B**H + beta*C
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + (*a)[l-1+(i-1)*(*lda)+(*aoff)]*conjc128((*b)[j-1+(l-1)*(*ldb)+(*boff)])
					}
					if (*beta) == zero {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha) * temp
					} else {
						(*c)[i-1+(j-1)*(*ldc)+(*coff)] = (*alpha)*temp + (*beta)*(*c)[i-1+(j-1)*(*ldc)+(*coff)]
					}
				}
			}
		} else {
			//
			//           Form  C := alpha*A**T*B**T + beta*C
			//
			for j = 1; j <= (*n); j++ {
				for i = 1; i <= (*m); i++ {
					temp = zero
					for l = 1; l <= (*k); l++ {
						temp = temp + (*a)[l-1+(i-1)*(*lda)+(*aoff)]*(*b)[j-1+(l-1)*(*ldb)+(*boff)]
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
