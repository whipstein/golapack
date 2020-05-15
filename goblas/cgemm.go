package goblas

import 
// \brief \b Cgemm
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Cgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
//
//       .. Scalar Arguments ..
//       COMPLEX ALPHA,BETA
//       INTEGER K,LDA,LDB,LDC,M,N
//       CHARACTER TRANSA,TRANSB
//       ..
//       .. Array Arguments ..
//       COMPLEX A(LDA,//),B(LDB,//),C(LDC,//)
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
//    C := alpha//op( A)//op( B) + beta//C,
//
// where  op( X) is one of
//
//    op( X) = X   or   op( X) = X////T   or   op( X) = X////H,
//
// alpha and beta are scalars, and A, B and C are matrices, with op( A)
// an m by k matrix,  op( B)  a  k by n matrix and  C an m by n matrix.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] TRANSA
// \verbatim
//          TRANSA is CHARACTER//1
//           On entry, TRANSA specifies the form of op( A) to be used in
//           the matrix multiplication as follows:
//
//              TRANSA = 'N' or 'n',  op( A) = A.
//
//              TRANSA = 'T' or 't',  op( A) = A////T.
//
//              TRANSA = 'C' or 'c',  op( A) = A////H.
// \endverbatim
//
// \param[in] TRANSB
// \verbatim
//          TRANSB is CHARACTER//1
//           On entry, TRANSB specifies the form of op( B) to be used in
//           the matrix multiplication as follows:
//
//              TRANSB = 'N' or 'n',  op( B) = B.
//
//              TRANSB = 'T' or 't',  op( B) = B////T.
//
//              TRANSB = 'C' or 'c',  op( B) = B////H.
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is INTEGER
//           On entry,  M  specifies  the number  of rows  of the  matrix
//           op( A)  and of the  matrix  C.  M  must  be at least  zero.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry,  N  specifies the number  of columns of the matrix
//           op( B) and the number of columns of the matrix C. N must be
//           at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry,  K  specifies  the number of columns of the matrix
//           op( A) and the number of rows of the matrix op( B). K must
//           be at least  zero.
// \endverbatim
//
// \param[in] ALPHA
// \verbatim
//          ALPHA is COMPLEX
//           On entry, ALPHA specifies the scalar alpha.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX array, dimension ( LDA, ka), where ka is
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
//           LDA must be at least  max( 1, m), otherwise  LDA must be at
//           least  max( 1, k).
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is COMPLEX array, dimension ( LDB, kb), where kb is
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
//           LDB must be at least  max( 1, k), otherwise  LDB must be at
//           least  max( 1, n).
// \endverbatim
//
// \param[in] BETA
// \verbatim
//          BETA is COMPLEX
//           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
//           supplied as zero then C need not be set on input.
// \endverbatim
//
// \param[in,out] C
// \verbatim
//          C is COMPLEX array, dimension ( LDC, N)
//           Before entry, the leading  m by n  part of the array  C must
//           contain the matrix  C,  except when  beta  is zero, in which
//           case C need not be set on entry.
//           On exit, the array  C  is overwritten by the  m by n  matrix
//           ( alpha//op( A)//op( B) + beta//C).
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
func Cgemm(transa *byte, transb *byte, m *int, n *int, k *int, alpha *complex128, a *[][]complex128, lda *int, b *[][]complex128, ldb *int, beta *complex128, c *[][]complex128, ldc *int) {
	temp := new(complex128)
	i := new(int)
	info := new(int)
	j := new(int)
	l := new(int)
	ncola := new(int)
	nrowa := new(int)
	nrowb := new(int)
	conja := new(bool)
	conjb := new(bool)
	nota := new(bool)
	notb := new(bool)
	one := new(complex128)
	zero := new(complex128)
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
	(*one) = (1.0e+0 + (0.0e+0)*1i)
	(*zero) = (0.0e+0 + (0.0e+0)*1i)
	//*     ..
	//*
	//*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
	//*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
	//*     B  respectively are to be  transposed but  not conjugated  and set
	//*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
	//*     and the number of rows of  B  respectively.
	//*
	(*nota) = (*Lsame(transa, func() *byte{y := byte('n'); return &y}()))
	(*notb) = (*Lsame(transb, func() *byte{y := byte('n'); return &y}()))
	(*conja) = (*Lsame(transa, func() *byte{y := byte('c'); return &y}()))
	(*conjb) = (*Lsame(transb, func() *byte{y := byte('c'); return &y}()))
	if *nota {
		(*nrowa) = (*m)
		(*ncola) = (*k)
	} else {
		(*nrowa) = (*k)
		(*ncola) = (*m)
	}
	if *notb {
		(*nrowb) = (*k)
	} else {
		(*nrowb) = (*n)
	}
	//*
	//*     Test the input parameters.
	//*
	(*info) = 0
	if ( . !(*nota)) && ( . !(*conja)) && ( . !Lsame((*transa), "t")) {
		(*info) = 1
	} else if ( . !(*notb)) && ( . !(*conjb)) && ( . !Lsame((*transb), "t")) {
		(*info) = 2
	} else if (*m) < 0 {
		(*info) = 3
	} else if (*n) < 0 {
		(*info) = 4
	} else if (*k) < 0 {
		(*info) = 5
	} else if (*lda) < max(func() *int{y := 1; return &y}(), nrowa) {
		(*info) = 8
	} else if (*ldb) < max(func() *int{y := 1; return &y}(), nrowb) {
		(*info) = 10
	} else if (*ldc) < max(func() *int{y := 1; return &y}(), m) {
		(*info) = 13
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte{y := []byte("cgemm "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if ((*m) == 0) || ((*n) == 0) || ( ( ((*alpha) == (*zero)) || ((*k) == 0)) && ((*beta) == (*one))) {
		return
	}
	//*
	//*     And when  alpha.eq.zero.
	//*
	if (*alpha) == (*zero) {
		if (*beta) == (*zero) {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-1][(*j)-1] = (*zero)
				//Label10:
				}
			//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*c)[(*i)-1][(*j)-1] = (*beta) * (*c)[(*i)-1][(*j)-1]
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
	if *notb {
		if *nota {
			//*
			//*           Form  C := alpha*A*B + beta*C.
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*beta) == (*zero) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*zero)
					//Label50:
					}
				} else if (*beta) != (*one) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*beta) * (*c)[(*i)-1][(*j)-1]
					//Label60:
					}
				}
				for (*l) = 1; (*l) <= (*k); (*l)++ {
					(*temp) = (*alpha) * (*b)[(*l)-1][(*j)-1]
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*c)[(*i)-1][(*j)-1] + (*temp)*(*a)[(*i)-1][(*l)-1]
					//Label70:
					}
				//Label80:
				}
			//Label90:
			}
		} else if *conja {
			//*
			//*           Form  C := alpha*A**H*B + beta*C.
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + CONJG(((*a)[(*l)-1][(*i)-1]))*(*b)[(*l)-1][(*j)-1]
					//Label100:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label110:
				}
			//Label120:
			}
		} else {
			//*
			//*           Form  C := alpha*A**T*B + beta*C
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + (*a)[(*l)-1][(*i)-1]*(*b)[(*l)-1][(*j)-1]
					//Label130:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label140:
				}
			//Label150:
			}
		}
	} else if *nota {
		if *conjb {
			//*
			//*           Form  C := alpha*A*B**H + beta*C.
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*beta) == (*zero) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*zero)
					//Label160:
					}
				} else if (*beta) != (*one) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*beta) * (*c)[(*i)-1][(*j)-1]
					//Label170:
					}
				}
				for (*l) = 1; (*l) <= (*k); (*l)++ {
					(*temp) = (*alpha) * CONJG(((*b)[(*j)-1][(*l)-1]))
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*c)[(*i)-1][(*j)-1] + (*temp)*(*a)[(*i)-1][(*l)-1]
					//Label180:
					}
				//Label190:
				}
			//Label200:
			}
		} else {
			//*
			//*           Form  C := alpha*A*B**T + beta*C
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				if (*beta) == (*zero) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*zero)
					//Label210:
					}
				} else if (*beta) != (*one) {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*beta) * (*c)[(*i)-1][(*j)-1]
					//Label220:
					}
				}
				for (*l) = 1; (*l) <= (*k); (*l)++ {
					(*temp) = (*alpha) * (*b)[(*j)-1][(*l)-1]
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*c)[(*i)-1][(*j)-1] = (*c)[(*i)-1][(*j)-1] + (*temp)*(*a)[(*i)-1][(*l)-1]
					//Label230:
					}
				//Label240:
				}
			//Label250:
			}
		}
	} else if *conja {
		if *conjb {
			//*
			//*           Form  C := alpha*A**H*B**H + beta*C.
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + CONJG(((*a)[(*l)-1][(*i)-1]))*CONJG(((*b)[(*j)-1][(*l)-1]))
					//Label260:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label270:
				}
			//Label280:
			}
		} else {
			//*
			//*           Form  C := alpha*A**H*B**T + beta*C
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + CONJG(((*a)[(*l)-1][(*i)-1]))*(*b)[(*j)-1][(*l)-1]
					//Label290:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label300:
				}
			//Label310:
			}
		}
	} else {
		if *conjb {
			//*
			//*           Form  C := alpha*A**T*B**H + beta*C
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + (*a)[(*l)-1][(*i)-1]*CONJG(((*b)[(*j)-1][(*l)-1]))
					//Label320:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label330:
				}
			//Label340:
			}
		} else {
			//*
			//*           Form  C := alpha*A**T*B**T + beta*C
			//*
			for (*j) = 1; (*j) <= (*n); (*j)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*temp) = (*zero)
					for (*l) = 1; (*l) <= (*k); (*l)++ {
						(*temp) = (*temp) + (*a)[(*l)-1][(*i)-1]*(*b)[(*j)-1][(*l)-1]
					//Label350:
					}
					if (*beta) == (*zero) {
						(*c)[(*i)-1][(*j)-1] = (*alpha) * (*temp)
					} else {
						(*c)[(*i)-1][(*j)-1] = (*alpha)*(*temp) + (*beta)*(*c)[(*i)-1][(*j)-1]
					}
				//Label360:
				}
			//Label370:
			}
		}
	}
	//*
	return
	//*
	//*     End of Cgemm .
	//*
}
