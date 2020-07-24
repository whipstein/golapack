package golapack

// Ztrsm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
//
//       .. Scalar Arguments ..
//       COMPLEX*16 ALPHA
//       INTEGER LDA,LDB,M,N
//       CHARACTER DIAG,SIDE,TRANSA,UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX*16 A(LDA,*),B(LDB,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// ZTRSM  solves one of the matrix equations
//
//    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
//
// where alpha is a scalar, X and B are m by n matrices, A is a unit, or
// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
//
//    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
//
// The matrix X is overwritten on B.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SIDE
// \verbatim
//          SIDE is CHARACTER*1
//           On entry, SIDE specifies whether op( A ) appears on the left
//           or right of X as follows:
//
//              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
//
//              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
// \endverbatim
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the matrix A is an upper or
//           lower triangular matrix as follows:
//
//              UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//              UPLO = 'L' or 'l'   A is a lower triangular matrix.
// \endverbatim
//
// \param[in] TRANSA
// \verbatim
//          TRANSA is CHARACTER*1
//           On entry, TRANSA specifies the form of op( A ) to be used in
//           the matrix multiplication as follows:
//
//              TRANSA = 'N' or 'n'   op( A ) = A.
//
//              TRANSA = 'T' or 't'   op( A ) = A**T.
//
//              TRANSA = 'C' or 'c'   op( A ) = A**H.
// \endverbatim
//
// \param[in] DIAG
// \verbatim
//          DIAG is CHARACTER*1
//           On entry, DIAG specifies whether or not A is unit triangular
//           as follows:
//
//              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
//
//              DIAG = 'N' or 'n'   A is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is INTEGER
//           On entry, M specifies the number of rows of B. M must be at
//           least zero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the number of columns of B.  N must be
//           at least zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is COMPLEX*16
//           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
//           zero then  A is not referenced and  B need not be set before
//           entry.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX*16 array, dimension ( LDA, k ),
//           where k is m when SIDE = 'L' or 'l'
//             and k is n when SIDE = 'R' or 'r'.
//           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
//           upper triangular part of the array  A must contain the upper
//           triangular matrix  and the strictly lower triangular part of
//           A is not referenced.
//           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
//           lower triangular part of the array  A must contain the lower
//           triangular matrix  and the strictly upper triangular part of
//           A is not referenced.
//           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
//           A  are not referenced either,  but are assumed to be  unity.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
//           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
//           then LDA must be at least max( 1, n ).
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is COMPLEX*16 array, dimension ( LDB, N )
//           Before entry,  the leading  m by n part of the array  B must
//           contain  the  right-hand  side  matrix  B,  and  on exit  is
//           overwritten by the solution matrix  X.
// \endverbatim
//
// \param[in] LDB
// \verbatim
//          LDB is INTEGER
//           On entry, LDB specifies the first dimension of B as declared
//           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
func Ztrsm(side *byte, uplo *byte, transa *byte, diag *byte, m *int, n *int, alpha *complex128, a *[]complex128, aoff, lda *int, b *[]complex128, boff, ldb *int) {
	var lside, noconj, nounit, upper bool
	var temp complex128
	var i, info, j, k, nrowa int
	var one complex128 = (1.0 + 0.0*1i)
	var zero complex128 = (0.0 + 0.0*1i)

	//
	//     Test the input parameters.
	//
	lside = Lsame(side, func() *byte { y := byte('L'); return &y }())
	if lside {
		nrowa = (*m)
	} else {
		nrowa = (*n)
	}
	noconj = Lsame(transa, func() *byte { y := byte('T'); return &y }())
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	upper = Lsame(uplo, func() *byte { y := byte('U'); return &y }())
	//
	info = 0
	if (!lside) && (!Lsame(side, func() *byte { y := byte('R'); return &y }())) {
		info = 1
	} else if (!upper) && (!Lsame(uplo, func() *byte { y := byte('L'); return &y }())) {
		info = 2
	} else if (!Lsame(transa, func() *byte { y := byte('N'); return &y }())) && (!Lsame(transa, func() *byte { y := byte('T'); return &y }())) && (!Lsame(transa, func() *byte { y := byte('C'); return &y }())) {
		info = 3
	} else if (!Lsame(diag, func() *byte { y := byte('U'); return &y }())) && (!Lsame(diag, func() *byte { y := byte('N'); return &y }())) {
		info = 4
	} else if (*m) < 0 {
		info = 5
	} else if (*n) < 0 {
		info = 6
	} else if (*lda) < maxint(1, nrowa) {
		info = 9
	} else if (*ldb) < maxint(1, *m) {
		info = 11
	}
	if info != 0 {
		Xerbla(func() *[]byte { y := []byte("ZTRSM "); return &y }(), &info)
		return
	}
	//
	//     Quick return if possible.
	//
	if (*m) == 0 || (*n) == 0 {
		return
	}
	//
	//     And when  alpha.eq.zero.
	//
	if (*alpha) == zero {
		for j = 1; j <= (*n); j++ {
			for i = 1; i <= (*m); i++ {
				(*b)[i-1+(j-1)*(*ldb)+(*boff)] = zero
			}
		}
		return
	}
	//
	//     Start the operations.
	//
	if lside {
		if Lsame(transa, func() *byte { y := byte('N'); return &y }()) {
			//
			//           Form  B := alpha*inv( A )*B.
			//
			if upper {
				for j = 1; j <= (*n); j++ {
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
					for k = (*m); k >= 1; k-- {
						if (*b)[k-1+(j-1)*(*ldb)+(*boff)] != zero {
							if nounit {
								(*b)[k-1+(j-1)*(*ldb)+(*boff)] = (*b)[k-1+(j-1)*(*ldb)+(*boff)] / (*a)[k-1+(k-1)*(*lda)+(*aoff)]
							}
							for i = 1; i <= k-1; i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - (*b)[k-1+(j-1)*(*ldb)+(*boff)]*(*a)[i-1+(k-1)*(*lda)+(*aoff)]
							}
						}
					}
				}
			} else {
				for j = 1; j <= (*n); j++ {
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
					for k = 1; k <= (*m); k++ {
						if (*b)[k-1+(j-1)*(*ldb)+(*boff)] != zero {
							if nounit {
								(*b)[k-1+(j-1)*(*ldb)+(*boff)] = (*b)[k-1+(j-1)*(*ldb)+(*boff)] / (*a)[k-1+(k-1)*(*lda)+(*aoff)]
							}
							for i = k + 1; i <= (*m); i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - (*b)[k-1+(j-1)*(*ldb)+(*boff)]*(*a)[i-1+(k-1)*(*lda)+(*aoff)]
							}
						}
					}
				}
			}
		} else {
			//
			//           Form  B := alpha*inv( A**T )*B
			//           or    B := alpha*inv( A**H )*B.
			//
			if upper {
				for j = 1; j <= (*n); j++ {
					for i = 1; i <= (*m); i++ {
						temp = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						if noconj {
							for k = 1; k <= i-1; k++ {
								temp = temp - (*a)[k-1+(i-1)*(*lda)+(*aoff)]*(*b)[k-1+(j-1)*(*ldb)+(*boff)]
							}
							if nounit {
								temp = temp / (*a)[i-1+(i-1)*(*lda)+(*aoff)]
							}
						} else {
							for k = 1; k <= i-1; k++ {
								temp = temp - conjc128((*a)[k-1+(i-1)*(*lda)+(*aoff)])*(*b)[k-1+(j-1)*(*ldb)+(*boff)]
							}
							if nounit {
								temp = temp / conjc128((*a)[i-1+(i-1)*(*lda)+(*aoff)])
							}
						}
						(*b)[i-1+(j-1)*(*ldb)+(*boff)] = temp
					}
				}
			} else {
				for j = 1; j <= (*n); j++ {
					for i = (*m); i >= 1; i-- {
						temp = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						if noconj {
							for k = i + 1; k <= (*m); k++ {
								temp = temp - (*a)[k-1+(i-1)*(*lda)+(*aoff)]*(*b)[k-1+(j-1)*(*ldb)+(*boff)]
							}
							if nounit {
								temp = temp / (*a)[i-1+(i-1)*(*lda)+(*aoff)]
							}
						} else {
							for k = i + 1; k <= (*m); k++ {
								temp = temp - conjc128((*a)[k-1+(i-1)*(*lda)+(*aoff)])*(*b)[k-1+(j-1)*(*ldb)+(*boff)]
							}
							if nounit {
								temp = temp / conjc128((*a)[i-1+(i-1)*(*lda)+(*aoff)])
							}
						}
						(*b)[i-1+(j-1)*(*ldb)+(*boff)] = temp
					}
				}
			}
		}
	} else {
		if Lsame(transa, func() *byte { y := byte('N'); return &y }()) {
			//
			//           Form  B := alpha*B*inv( A ).
			//
			if upper {
				for j = 1; j <= (*n); j++ {
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
					for k = 1; k <= j-1; k++ {
						if (*a)[k-1+(j-1)*(*lda)+(*aoff)] != zero {
							for i = 1; i <= (*m); i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - (*a)[k-1+(j-1)*(*lda)+(*aoff)]*(*b)[i-1+(k-1)*(*ldb)+(*boff)]
							}
						}
					}
					if nounit {
						temp = one / (*a)[j-1+(j-1)*(*lda)+(*aoff)]
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = temp * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
				}
			} else {
				for j = (*n); j >= 1; j-- {
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
					for k = j + 1; k <= (*n); k++ {
						if (*a)[k-1+(j-1)*(*lda)+(*aoff)] != zero {
							for i = 1; i <= (*m); i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - (*a)[k-1+(j-1)*(*lda)+(*aoff)]*(*b)[i-1+(k-1)*(*ldb)+(*boff)]
							}
						}
					}
					if nounit {
						temp = one / (*a)[j-1+(j-1)*(*lda)+(*aoff)]
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(j-1)*(*ldb)+(*boff)] = temp * (*b)[i-1+(j-1)*(*ldb)+(*boff)]
						}
					}
				}
			}
		} else {
			//
			//           Form  B := alpha*B*inv( A**T )
			//           or    B := alpha*B*inv( A**H ).
			//
			if upper {
				for k = (*n); k >= 1; k-- {
					if nounit {
						if noconj {
							temp = one / (*a)[k-1+(k-1)*(*lda)+(*aoff)]
						} else {
							temp = one / conjc128((*a)[k-1+(k-1)*(*lda)+(*aoff)])
						}
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(k-1)*(*ldb)+(*boff)] = temp * (*b)[i-1+(k-1)*(*ldb)+(*boff)]
						}
					}
					for j = 1; j <= k-1; j++ {
						if (*a)[j-1+(k-1)*(*lda)+(*aoff)] != zero {
							if noconj {
								temp = (*a)[j-1+(k-1)*(*lda)+(*aoff)]
							} else {
								temp = conjc128((*a)[j-1+(k-1)*(*lda)+(*aoff)])
							}
							for i = 1; i <= (*m); i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - temp*(*b)[i-1+(k-1)*(*ldb)+(*boff)]
							}
						}
					}
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(k-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(k-1)*(*ldb)+(*boff)]
						}
					}
				}
			} else {
				for k = 1; k <= (*n); k++ {
					if nounit {
						if noconj {
							temp = one / (*a)[k-1+(k-1)*(*lda)+(*aoff)]
						} else {
							temp = one / conjc128((*a)[k-1+(k-1)*(*lda)+(*aoff)])
						}
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(k-1)*(*ldb)+(*boff)] = temp * (*b)[i-1+(k-1)*(*ldb)+(*boff)]
						}
					}
					for j = k + 1; j <= (*n); j++ {
						if (*a)[j-1+(k-1)*(*lda)+(*aoff)] != zero {
							if noconj {
								temp = (*a)[j-1+(k-1)*(*lda)+(*aoff)]
							} else {
								temp = conjc128((*a)[j-1+(k-1)*(*lda)+(*aoff)])
							}
							for i = 1; i <= (*m); i++ {
								(*b)[i-1+(j-1)*(*ldb)+(*boff)] = (*b)[i-1+(j-1)*(*ldb)+(*boff)] - temp*(*b)[i-1+(k-1)*(*ldb)+(*boff)]
							}
						}
					}
					if (*alpha) != one {
						for i = 1; i <= (*m); i++ {
							(*b)[i-1+(k-1)*(*ldb)+(*boff)] = (*alpha) * (*b)[i-1+(k-1)*(*ldb)+(*boff)]
						}
					}
				}
			}
		}
	}
}
