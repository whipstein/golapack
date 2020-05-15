package goblas
// \brief \b Strsm
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Strsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
//
//       .. Scalar Arguments ..
//       REAL ALPHA
//       INTEGER LDA,LDB,M,N
//       CHARACTER DIAG,SIDE,TRANSA,UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,//),B(LDB,//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Strsm  solves one of the matrix equations
//
//    op( A)//X = alpha//B,   or   X//op( A) = alpha//B,
//
// where alpha is a scalar, X and B are m by n matrices, A is a unit, or
// non-unit,  upper or lower triangular matrix  and  op( A)  is one  of
//
//    op( A) = A   or   op( A) = A////T.
//
// The matrix X is overwritten on B.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SIDE
// \verbatim
//          SIDE is CHARACTER//1
//           On entry, SIDE specifies whether op( A) appears on the left
//           or right of X as follows:
//
//              SIDE = 'L' or 'l'   op( A)//X = alpha//B.
//
//              SIDE = 'R' or 'r'   X//op( A) = alpha//B.
// \endverbatim
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
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
//          TRANSA is CHARACTER//1
//           On entry, TRANSA specifies the form of op( A) to be used in
//           the matrix multiplication as follows:
//
//              TRANSA = 'N' or 'n'   op( A) = A.
//
//              TRANSA = 'T' or 't'   op( A) = A////T.
//
//              TRANSA = 'C' or 'c'   op( A) = A////T.
// \endverbatim
//
// \param[in] DIAG
// \verbatim
//          DIAG is CHARACTER//1
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
//          ALPHA is REAL
//           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
//           zero then  A is not referenced and  B need not be set before
//           entry.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, k),
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
//           LDA  must be at least  max( 1, m),  when  SIDE = 'R' or 'r'
//           then LDA must be at least max( 1, n).
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is REAL array, dimension ( LDB, N)
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
//           max( 1, m).
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
func Strsm(side *byte, uplo *byte, transa *byte, diag *byte, m *int, n *int, alpha *float64, a *[][]float64, lda *int, b *[][]float64, ldb *int) {
	temp := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	k := new(int)
	nrowa := new(int)
	lside := new(bool)
	nounit := new(bool)
	upper := new(bool)
	one := new(float64)
	zero := new(float64)
	//*
	//*  -- Reference BLAS level3 routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     December 2016
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//*  =====================================================================
	//*
	//*     .. External Functions ..
	//*     ..
	//*     .. External Subroutines ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Parameters ..
	(*one) = 1.0e+0
	(*zero) = 0.0e+0
	//*     ..
	//*
	//*     Test the input parameters.
	//*
	(*lside) = (*Lsame(side, func() *byte{y := byte('l'); return &y}()))
	if *lside {
		(*nrowa) = (*m)
	} else {
		(*nrowa) = (*n)
	}
	(*nounit) = (*Lsame(diag, func() *byte{y := byte('n'); return &y}()))
	(*upper) = (*Lsame(uplo, func() *byte{y := byte('u'); return &y}()))
	//*
	(*info) = 0
	if ( . !(*lside)) && ( . !Lsame((*side), "r")) {
		(*info) = 1
	} else if ( . !(*upper)) && ( . !Lsame((*uplo), "l")) {
		(*info) = 2
	} else if ( . !Lsame((*transa), "n")) && ( . !Lsame((*transa), "t")) && ( . !Lsame((*transa), "c")) {
		(*info) = 3
	} else if ( . !Lsame((*diag), "u")) && ( . !Lsame((*diag), "n")) {
		(*info) = 4
	} else if (*m) < 0 {
		(*info) = 5
	} else if (*n) < 0 {
		(*info) = 6
	} else if (*lda) < max(func() *int{y := 1; return &y}(), nrowa) {
		(*info) = 9
	} else if (*ldb) < max(func() *int{y := 1; return &y}(), m) {
		(*info) = 11
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte{y := []byte("strsm "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if (*m) == 0 || (*n) == 0 {
		return
	}
	//*
	//*     And when  alpha.eq.zero.
	//*
	if (*alpha) == (*zero) {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*b)[(*i)-1][(*j)-1] = (*zero)
			//Label10:
			}
		//Label20:
		}
		return
	}
	//*
	//*     Start the operations.
	//*
	if *lside {
		if Lsame(transa, func() *byte{y := byte('n'); return &y}()) {
			//*
			//*           Form  B := alpha*inv( A)*B.
			//*
			if *upper {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*alpha) * (*b)[(*i)-1][(*j)-1]
						//Label30:
						}
					}
					for (*k) = (*m); (*k) <= 1; (*k) += -1 {
						if (*b)[(*k)-1][(*j)-1] != (*zero) {
							if *nounit {
								(*b)[(*k)-1][(*j)-1] = (*b)[(*k)-1][(*j)-1] / (*a)[(*k)-1][(*k)-1]
							}
							for (*i) = 1; (*i) <= (*k)-1; (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*b)[(*k)-1][(*j)-1]*(*a)[(*i)-1][(*k)-1]
							//Label40:
							}
						}
					//Label50:
					}
				//Label60:
				}
			} else {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*alpha) * (*b)[(*i)-1][(*j)-1]
						//Label70:
						}
					}
					for (*k) = 1; (*k) <= (*m); (*k)++ {
						if (*b)[(*k)-1][(*j)-1] != (*zero) {
							if *nounit {
								(*b)[(*k)-1][(*j)-1] = (*b)[(*k)-1][(*j)-1] / (*a)[(*k)-1][(*k)-1]
							}
							for (*i) = (*k) + 1; (*i) <= (*m); (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*b)[(*k)-1][(*j)-1]*(*a)[(*i)-1][(*k)-1]
							//Label80:
							}
						}
					//Label90:
					}
				//Label100:
				}
			}
		} else {
			//*
			//*           Form  B := alpha*inv( A**T)*B.
			//*
			if *upper {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*temp) = (*alpha) * (*b)[(*i)-1][(*j)-1]
						for (*k) = 1; (*k) <= (*i)-1; (*k)++ {
							(*temp) = (*temp) - (*a)[(*k)-1][(*i)-1]*(*b)[(*k)-1][(*j)-1]
						//Label110:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[(*i)-1][(*i)-1]
						}
						(*b)[(*i)-1][(*j)-1] = (*temp)
					//Label120:
					}
				//Label130:
				}
			} else {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					for (*i) = (*m); (*i) <= 1; (*i) += -1 {
						(*temp) = (*alpha) * (*b)[(*i)-1][(*j)-1]
						for (*k) = (*i) + 1; (*k) <= (*m); (*k)++ {
							(*temp) = (*temp) - (*a)[(*k)-1][(*i)-1]*(*b)[(*k)-1][(*j)-1]
						//Label140:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[(*i)-1][(*i)-1]
						}
						(*b)[(*i)-1][(*j)-1] = (*temp)
					//Label150:
					}
				//Label160:
				}
			}
		}
	} else {
		if Lsame(transa, func() *byte{y := byte('n'); return &y}()) {
			//*
			//*           Form  B := alpha*B*inv( A).
			//*
			if *upper {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*alpha) * (*b)[(*i)-1][(*j)-1]
						//Label170:
						}
					}
					for (*k) = 1; (*k) <= (*j)-1; (*k)++ {
						if (*a)[(*k)-1][(*j)-1] != (*zero) {
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*a)[(*k)-1][(*j)-1]*(*b)[(*i)-1][(*k)-1]
							//Label180:
							}
						}
					//Label190:
					}
					if *nounit {
						(*temp) = (*one) / (*a)[(*j)-1][(*j)-1]
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*temp) * (*b)[(*i)-1][(*j)-1]
						//Label200:
						}
					}
				//Label210:
				}
			} else {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*alpha) * (*b)[(*i)-1][(*j)-1]
						//Label220:
						}
					}
					for (*k) = (*j) + 1; (*k) <= (*n); (*k)++ {
						if (*a)[(*k)-1][(*j)-1] != (*zero) {
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*a)[(*k)-1][(*j)-1]*(*b)[(*i)-1][(*k)-1]
							//Label230:
							}
						}
					//Label240:
					}
					if *nounit {
						(*temp) = (*one) / (*a)[(*j)-1][(*j)-1]
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*j)-1] = (*temp) * (*b)[(*i)-1][(*j)-1]
						//Label250:
						}
					}
				//Label260:
				}
			}
		} else {
			//*
			//*           Form  B := alpha*B*inv( A**T).
			//*
			if *upper {
				for (*k) = (*n); (*k) <= 1; (*k) += -1 {
					if *nounit {
						(*temp) = (*one) / (*a)[(*k)-1][(*k)-1]
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*k)-1] = (*temp) * (*b)[(*i)-1][(*k)-1]
						//Label270:
						}
					}
					for (*j) = 1; (*j) <= (*k)-1; (*j)++ {
						if (*a)[(*j)-1][(*k)-1] != (*zero) {
							(*temp) = (*a)[(*j)-1][(*k)-1]
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*temp)*(*b)[(*i)-1][(*k)-1]
							//Label280:
							}
						}
					//Label290:
					}
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*k)-1] = (*alpha) * (*b)[(*i)-1][(*k)-1]
						//Label300:
						}
					}
				//Label310:
				}
			} else {
				for (*k) = 1; (*k) <= (*n); (*k)++ {
					if *nounit {
						(*temp) = (*one) / (*a)[(*k)-1][(*k)-1]
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*k)-1] = (*temp) * (*b)[(*i)-1][(*k)-1]
						//Label320:
						}
					}
					for (*j) = (*k) + 1; (*j) <= (*n); (*j)++ {
						if (*a)[(*j)-1][(*k)-1] != (*zero) {
							(*temp) = (*a)[(*j)-1][(*k)-1]
							for (*i) = 1; (*i) <= (*m); (*i)++ {
								(*b)[(*i)-1][(*j)-1] = (*b)[(*i)-1][(*j)-1] - (*temp)*(*b)[(*i)-1][(*k)-1]
							//Label330:
							}
						}
					//Label340:
					}
					if (*alpha) != (*one) {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*b)[(*i)-1][(*k)-1] = (*alpha) * (*b)[(*i)-1][(*k)-1]
						//Label350:
						}
					}
				//Label360:
				}
			}
		}
	}
	//*
	return
	//*
	//*     End of Strsm .
	//*
}