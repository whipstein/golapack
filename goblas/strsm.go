package goblas

// Strsm ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE STRSM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
//
//       .. Scalar Arguments ..
//       REAL alpha
//       INTEGER lda,ldb,m,n
//       CHARACTER diag,side,transa,uplo
//       ..
//       .. Array Arguments ..
//       REAL a(lda,*),b(ldb,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// STRSM  solves 1.0 of the matrix equations
//
//    op( a )*X = alpha*b,   or   X*op( a ) = alpha*b,
//
// where alpha is a scalar, X and b are m by n matrices, a is a unit, or
// non-unit,  upper or lower triangular matrix  and  op( a )  is 1.0  of
//
//    op( a ) = a   or   op( a ) = a**T.
//
// The matrix X is overwritten on b.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] side
// \verbatim
//          side is CHARACTER*1
//           On entry, side specifies whether op( a ) appears on the left
//           or right of X as follows:
//
//              side = 'L' or 'l'   op( a )*X = alpha*b.
//
//              side = 'R' or 'r'   X*op( a ) = alpha*b.
// \endverbatim
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On entry, uplo specifies whether the matrix a is an upper or
//           lower triangular matrix as follows:
//
//              uplo = 'U' or 'u'   a is an upper triangular matrix.
//
//              uplo = 'L' or 'l'   a is a lower triangular matrix.
// \endverbatim
//
// \param[in] transa
// \verbatim
//          transa is CHARACTER*1
//           On entry, transa specifies the form of op( a ) to be used in
//           the matrix multiplication as follows:
//
//              transa = 'N' or 'N'   op( a ) = a.
//
//              transa = 'T' or 't'   op( a ) = a**T.
//
//              transa = 'C' or 'c'   op( a ) = a**T.
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//           On entry, diag specifies whether or not a is unit triangular
//           as follows:
//
//              diag = 'U' or 'u'   a is assumed to be unit triangular.
//
//              diag = 'N' or 'N'   a is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] m
// \verbatim
//          m is INTEGER
//           On entry, m specifies the number of rows of b. m must be at
//           least 0.0.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry, n specifies the number of columns of b.  n must be
//           at least 0.0.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is REAL
//           On entry,  alpha specifies the scalar  alpha. When  alpha is
//           0.0 then  a is not referenced and  b need not be set before
//           entry.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is REAL array, dimension ( lda, k ),
//           where k is m when side = 'L' or 'l'
//             and k is n when side = 'R' or 'r'.
//           Before entry  with  uplo = 'U' or 'u',  the  leading  k by k
//           upper triangular part of the array  a must contain the upper
//           triangular matrix  and the strictly lower triangular part of
//           a is not referenced.
//           Before entry  with  uplo = 'L' or 'l',  the  leading  k by k
//           lower triangular part of the array  a must contain the lower
//           triangular matrix  and the strictly upper triangular part of
//           a is not referenced.
//           Note that when  diag = 'U' or 'u',  the diagonal elements of
//           a  are not referenced either,  but are assumed to be  unity.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in the calling (sub) program.  When  side = 'L' or 'l'  then
//           lda  must be at least  max( 1, m ),  when  side = 'R' or 'r'
//           then lda must be at least max( 1, n ).
// \endverbatim
//
// \param[in,out] b
// \verbatim
//          b is REAL array, dimension ( ldb, n )
//           Before entry,  the leading  m by n part of the array  b must
//           contain  the  right-hand  side  matrix  b,  and  on exit  is
//           overwritten by the solution matrix  X.
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
//
//  -- Written on 8-February-1989.
//     Jack Dongarra, Argonne National Laboratory.
//     Iain Duff, AERE Harwell.
//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
//     Sven Hammarling, Numerical Algorithms Group Ltd.
// \endverbatim
//
//  =====================================================================
func Strsm(major, side, uplo, transa, diag *byte, m, n *int, alpha *float32, a *[][]float32, lda *int, b *[][]float32, ldb *int) {
	var temp float32
	var i, info, j, k, nrowa int
	var lside, nounit, upper bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     Test the input parameters.
	//
	lside = Lsame(side, func() *byte { y := byte('L'); return &y }())
	if lside {
		nrowa = *m
	} else {
		nrowa = *n
	}
	nounit = Lsame(diag, func() *byte { y := byte('N'); return &y }())
	upper = Lsame(uplo, func() *byte { y := byte('U'); return &y }())

	if !Lsame(major, func() *byte { y := byte('C'); return &y }()) && !Lsame(major, func() *byte { y := byte('R'); return &y }()) {
		info = 1
	} else if !lside && !Lsame(side, func() *byte { y := byte('R'); return &y }()) {
		info = 2
	} else if !upper && !Lsame(uplo, func() *byte { y := byte('L'); return &y }()) {
		info = 3
	} else if !Lsame(transa, func() *byte { y := byte('N'); return &y }()) && !Lsame(transa, func() *byte { y := byte('T'); return &y }()) && !Lsame(transa, func() *byte { y := byte('C'); return &y }()) {
		info = 4
	} else if !Lsame(diag, func() *byte { y := byte('U'); return &y }()) && !Lsame(diag, func() *byte { y := byte('N'); return &y }()) {
		info = 5
	} else if *m < 0 {
		info = 6
	} else if *n < 0 {
		info = 7
	} else if *lda < max(1, nrowa) {
		info = 10
	} else if *ldb < max(1, *m) {
		info = 12
	}
	if info != 0 {
		name := "Strsm"
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
	if *m == 0 || *n == 0 {
		return
	}
	//
	//     And when  alpha.eq.0.0.
	//
	if *alpha == 0.0 {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*b)[i-1][j-1] = 0.0
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
			//           Form  b := alpha*inv( a )*b.
			//
			if upper {
				for j = 1; j <= *n; j++ {
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= *alpha
						}
					}
					for k = *m; k >= 1; k-- {
						if (*b)[k-1][j-1] != 0.0 {
							if nounit {
								(*b)[k-1][j-1] /= (*a)[k-1][k-1]
							}
							for i = 1; i <= k-1; i++ {
								(*b)[i-1][j-1] -= (*b)[k-1][j-1] * (*a)[i-1][k-1]
							}
						}
					}
				}
			} else {
				for j = 1; j <= *n; j++ {
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= *alpha
						}
					}
					for k = 1; k <= *m; k++ {
						if (*b)[k-1][j-1] != 0.0 {
							if nounit {
								(*b)[k-1][j-1] /= (*a)[k-1][k-1]
							}
							for i = k + 1; i <= *m; i++ {
								(*b)[i-1][j-1] -= (*b)[k-1][j-1] * (*a)[i-1][k-1]
							}
						}
					}
				}
			}
		} else {
			//
			//           Form  b := alpha*inv( a**T )*b.
			//
			if upper {
				for j = 1; j <= *n; j++ {
					for i = 1; i <= *m; i++ {
						temp = (*alpha) * (*b)[i-1][j-1]
						for k = 1; k <= i-1; k++ {
							temp -= (*a)[k-1][i-1] * (*b)[k-1][j-1]
						}
						if nounit {
							temp /= (*a)[i-1][i-1]
						}
						(*b)[i-1][j-1] = temp
					}
				}
			} else {
				for j = 1; j <= *n; j++ {
					for i = *m; i >= 1; i-- {
						temp = (*alpha) * (*b)[i-1][j-1]
						for k = i + 1; k <= *m; k++ {
							temp -= (*a)[k-1][i-1] * (*b)[k-1][j-1]
						}
						if nounit {
							temp /= (*a)[i-1][i-1]
						}
						(*b)[i-1][j-1] = temp
					}
				}
			}
		}
	} else {
		if Lsame(transa, func() *byte { y := byte('N'); return &y }()) {
			//
			//           Form  b := alpha*b*inv( a ).
			//
			if upper {
				for j = 1; j <= *n; j++ {
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= *alpha
						}
					}
					for k = 1; k <= j-1; k++ {
						if (*a)[k-1][j-1] != 0.0 {
							for i = 1; i <= *m; i++ {
								(*b)[i-1][j-1] -= (*a)[k-1][j-1] * (*b)[i-1][k-1]
							}
						}
					}
					if nounit {
						temp = 1.0 / (*a)[j-1][j-1]
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= temp
						}
					}
				}
			} else {
				for j = *n; j >= 1; j-- {
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= *alpha
						}
					}
					for k = j + 1; k <= *n; k++ {
						if (*a)[k-1][j-1] != 0.0 {
							for i = 1; i <= *m; i++ {
								(*b)[i-1][j-1] -= (*a)[k-1][j-1] * (*b)[i-1][k-1]
							}
						}
					}
					if nounit {
						temp = 1.0 / (*a)[j-1][j-1]
						for i = 1; i <= *m; i++ {
							(*b)[i-1][j-1] *= temp
						}
					}
				}
			}
		} else {
			//
			//           Form  b := alpha*b*inv( a**T ).
			//
			if upper {
				for k = *n; k >= 1; k-- {
					if nounit {
						temp = 1.0 / (*a)[k-1][k-1]
						for i = 1; i <= *m; i++ {
							(*b)[i-1][k-1] *= temp
						}
					}
					for j = 1; j <= k-1; j++ {
						if (*a)[j-1][k-1] != 0.0 {
							temp = (*a)[j-1][k-1]
							for i = 1; i <= *m; i++ {
								(*b)[i-1][j-1] -= temp * (*b)[i-1][k-1]
							}
						}
					}
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][k-1] *= *alpha
						}
					}
				}
			} else {
				for k = 1; k <= *n; k++ {
					if nounit {
						temp = 1.0 / (*a)[k-1][k-1]
						for i = 1; i <= *m; i++ {
							(*b)[i-1][k-1] *= temp
						}
					}
					for j = k + 1; j <= *n; j++ {
						if (*a)[j-1][k-1] != 0.0 {
							temp = (*a)[j-1][k-1]
							for i = 1; i <= *m; i++ {
								(*b)[i-1][j-1] -= temp * (*b)[i-1][k-1]
							}
						}
					}
					if *alpha != 1.0 {
						for i = 1; i <= *m; i++ {
							(*b)[i-1][k-1] *= *alpha
						}
					}
				}
			}
		}
	}

	return
}
