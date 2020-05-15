package goblas
// \brief \b Ssymm
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ssymm(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
//
//       .. Scalar Arguments ..
//       REAL ALPHA,BETA
//       INTEGER LDA,LDB,LDC,M,N
//       CHARACTER SIDE,UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,//),B(LDB,//),C(LDC,//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ssymm  performs one of the matrix-matrix operations
//
//    C := alpha//A//B + beta//C,
//
// or
//
//    C := alpha//B//A + beta//C,
//
// where alpha and beta are scalars,  A is a symmetric matrix and  B and
// C are  m by n matrices.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SIDE
// \verbatim
//          SIDE is CHARACTER//1
//           On entry,  SIDE  specifies whether  the  symmetric matrix  A
//           appears on the  left or right  in the  operation as follows:
//
//              SIDE = 'L' or 'l'   C := alpha//A//B + beta//C,
//
//              SIDE = 'R' or 'r'   C := alpha//B//A + beta//C,
// \endverbatim
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
//           On  entry,   UPLO  specifies  whether  the  upper  or  lower
//           triangular  part  of  the  symmetric  matrix   A  is  to  be
//           referenced as follows:
//
//              UPLO = 'U' or 'u'   Only the upper triangular part of the
//                                  symmetric matrix is to be referenced.
//
//              UPLO = 'L' or 'l'   Only the lower triangular part of the
//                                  symmetric matrix is to be referenced.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is INTEGER
//           On entry,  M  specifies the number of rows of the matrix  C.
//           M  must be at least zero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the number of columns of the matrix C.
//           N  must be at least zero.
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
//          A is REAL array, dimension ( LDA, ka), where ka is
//           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
//           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
//           the array  A  must contain the  symmetric matrix,  such that
//           when  UPLO = 'U' or 'u', the leading m by m upper triangular
//           part of the array  A  must contain the upper triangular part
//           of the  symmetric matrix and the  strictly  lower triangular
//           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
//           the leading  m by m  lower triangular part  of the  array  A
//           must  contain  the  lower triangular part  of the  symmetric
//           matrix and the  strictly upper triangular part of  A  is not
//           referenced.
//           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
//           the array  A  must contain the  symmetric matrix,  such that
//           when  UPLO = 'U' or 'u', the leading n by n upper triangular
//           part of the array  A  must contain the upper triangular part
//           of the  symmetric matrix and the  strictly  lower triangular
//           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
//           the leading  n by n  lower triangular part  of the  array  A
//           must  contain  the  lower triangular part  of the  symmetric
//           matrix and the  strictly upper triangular part of  A  is not
//           referenced.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
//           LDA must be at least  max( 1, m), otherwise  LDA must be at
//           least  max( 1, n).
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is REAL array, dimension ( LDB, N)
//           Before entry, the leading  m by n part of the array  B  must
//           contain the matrix B.
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
// \param[in] BETA
// \verbatim
//          BETA is REAL
//           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
//           supplied as zero then C need not be set on input.
// \endverbatim
//
// \param[in,out] C
// \verbatim
//          C is REAL array, dimension ( LDC, N)
//           Before entry, the leading  m by n  part of the array  C must
//           contain the matrix  C,  except when  beta  is zero, in which
//           case C need not be set on entry.
//           On exit, the array  C  is overwritten by the  m by n updated
//           matrix.
// \endverbatim
//
// \param[in] LDC
// \verbatim
//          LDC is INTEGER
//           On entry, LDC specifies the first dimension of C as declared
//           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
//  -- Written on 8-February-1989.
//     Jack Dongarra, Argonne National Laboratory.
//     Iain Duff, AERE Harwell.
//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
//     Sven Hammarling, Numerical Algorithms Group Ltd.
// \endverbatim
//
//  =====================================================================
func Ssymm(side *byte, uplo *byte, m *int, n *int, alpha *float64, a *[][]float64, lda *int, b *[][]float64, ldb *int, beta *float64, c *[][]float64, ldc *int) {
	temp1 := new(float64)
	temp2 := new(float64)
	i := new(int)
	info := new(int)
	j := new(int)
	k := new(int)
	nrowa := new(int)
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
	//*     Set NROWA as the number of rows of A.
	//*
	if Lsame(side, func() *byte{y := byte('l'); return &y }()) {
		(*nrowa) = (*m)
	} else {
		(*nrowa) = (*n)
	}
	(*upper) = (*Lsame(uplo, func() *byte{y := byte('u'); return &y }()))
	//*
	//*     Test the input parameters.
	//*
	(*info) = 0
	if ( . !Lsame((*side), "l")) && ( . !Lsame((*side), "r")) {
		(*info) = 1
	} else if ( . !(*upper)) && ( . !Lsame((*uplo), "l")) {
		(*info) = 2
	} else if (*m) < 0 {
		(*info) = 3
	} else if (*n) < 0 {
		(*info) = 4
	} else if (*lda) < max(func() *int{y := 1; return &y }(), nrowa) {
		(*info) = 7
	} else if (*ldb) < max(func() *int{y := 1; return &y }(), m) {
		(*info) = 9
	} else if (*ldc) < max(func() *int{y := 1; return &y }(), m) {
		(*info) = 12
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte{y := []byte("ssymm "); return &y }(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*m) == 0) || ((*n) == 0) || ( ((*alpha) == (*zero)) && ((*beta) == (*one))) {
		return
	}
	//*
	//*     And when  alpha.eq.zero.
	//*
	if (*alpha) == (*zero) {
		if (*beta) == (*zero) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*zero)
				//Label10:
				}
			//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*beta) * (*c)[(*i)-(1)][(*j)-(1)]
				//Label30:
				}
			//Label40:
			}
		}
		return
	}
	//*
	//*     Start the operations.
	//*
	if Lsame(side, func() *byte{y := byte('l'); return &y }()) {
		//*
		//*        Form  C := alpha*A*B + beta*C.
		//*
		if *upper {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp1) = (*alpha) * (*b)[(*i)-(1)][(*j)-(1)]
					(*temp2) = (*zero)
					for (*k) = 1; (*k) <= (*i)-1; (*k)++ {
						(*c)[(*k)-(1)][(*j)-(1)] = (*c)[(*k)-(1)][(*j)-(1)] + (*temp1)*(*a)[(*k)-(1)][(*i)-(1)]
						(*temp2) = (*temp2) + (*b)[(*k)-(1)][(*j)-(1)]*(*a)[(*k)-(1)][(*i)-(1)]
					//Label50:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-(1)][(*j)-(1)] = (*temp1)*(*a)[(*i)-(1)][(*i)-(1)] + (*alpha)*(*temp2)
					} else {
						(*c)[(*i)-(1)][(*j)-(1)] = (*beta)*(*c)[(*i)-(1)][(*j)-(1)] + (*temp1)*(*a)[(*i)-(1)][(*i)-(1)] + (*alpha)*(*temp2)
					}
				//Label60:
				}
			//Label70:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = (*m); (*i) <= 1; (*i) += -1 {
					(*temp1) = (*alpha) * (*b)[(*i)-(1)][(*j)-(1)]
					(*temp2) = (*zero)
					for (*k) = (*i) + 1; (*k) <= (*m); (*k)++ {
						(*c)[(*k)-(1)][(*j)-(1)] = (*c)[(*k)-(1)][(*j)-(1)] + (*temp1)*(*a)[(*k)-(1)][(*i)-(1)]
						(*temp2) = (*temp2) + (*b)[(*k)-(1)][(*j)-(1)]*(*a)[(*k)-(1)][(*i)-(1)]
					//Label80:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-(1)][(*j)-(1)] = (*temp1)*(*a)[(*i)-(1)][(*i)-(1)] + (*alpha)*(*temp2)
					} else {
						(*c)[(*i)-(1)][(*j)-(1)] = (*beta)*(*c)[(*i)-(1)][(*j)-(1)] + (*temp1)*(*a)[(*i)-(1)][(*i)-(1)] + (*alpha)*(*temp2)
					}
				//Label90:
				}
			//Label100:
			}
		}
	} else {
		//*
		//*        Form  C := alpha*B*A + beta*C.
		//*
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			(*temp1) = (*alpha) * (*a)[(*j)-(1)][(*j)-(1)]
			if (*beta) == (*zero) {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*temp1) * (*b)[(*i)-(1)][(*j)-(1)]
				//Label110:
				}
			} else {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*beta)*(*c)[(*i)-(1)][(*j)-(1)] + (*temp1)*(*b)[(*i)-(1)][(*j)-(1)]
				//Label120:
				}
			}
			for (*k) = 1; (*k) <= (*j)-1; (*k)++ {
				if *upper {
					(*temp1) = (*alpha) * (*a)[(*k)-(1)][(*j)-(1)]
				} else {
					(*temp1) = (*alpha) * (*a)[(*j)-(1)][(*k)-(1)]
				}
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*c)[(*i)-(1)][(*j)-(1)] + (*temp1)*(*b)[(*i)-(1)][(*k)-(1)]
				//Label130:
				}
			//Label140:
			}
			for (*k) = (*j) + 1; (*k) <= (*n); (*k)++ {
				if *upper {
					(*temp1) = (*alpha) * (*a)[(*j)-(1)][(*k)-(1)]
				} else {
					(*temp1) = (*alpha) * (*a)[(*k)-(1)][(*j)-(1)]
				}
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-(1)][(*j)-(1)] = (*c)[(*i)-(1)][(*j)-(1)] + (*temp1)*(*b)[(*i)-(1)][(*k)-(1)]
				//Label150:
				}
			//Label160:
			}
		//Label170:
		}
	}
	//*
	return
	//*
	//*     End of Ssymm .
	//*
}
