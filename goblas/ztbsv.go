package goblas

import 

// \brief \b Ztbsv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Ztbsv(UPLO,TRANS,DIAG,N,K,A,LDA,X,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,K,LDA,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       COMPLEX//16 A(LDA,//),X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Ztbsv  solves one of the systems of equations
//
//    A//x = b,   or   A////T//x = b,   or   A////H//x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular band matrix, with ( k + 1)
// diagonals.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] UPLO
// \verbatim
//          UPLO is CHARACTER//1
//           On entry, UPLO specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//              UPLO = 'L' or 'l'   A is a lower triangular matrix.
// \endverbatim
//
// \param[in] TRANS
// \verbatim
//          TRANS is CHARACTER//1
//           On entry, TRANS specifies the equations to be solved as
//           follows:
//
//              TRANS = 'N' or 'n'   A//x = b.
//
//              TRANS = 'T' or 't'   A////T//x = b.
//
//              TRANS = 'C' or 'c'   A////H//x = b.
// \endverbatim
//
// \param[in] DIAG
// \verbatim
//          DIAG is CHARACTER//1
//           On entry, DIAG specifies whether or not A is unit
//           triangular as follows:
//
//              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
//
//              DIAG = 'N' or 'n'   A is not assumed to be unit
//                                  triangular.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is INTEGER
//           On entry, N specifies the order of the matrix A.
//           N must be at least zero.
// \endverbatim
//
// \param[in] K
// \verbatim
//          K is INTEGER
//           On entry with UPLO = 'U' or 'u', K specifies the number of
//           super-diagonals of the matrix A.
//           On entry with UPLO = 'L' or 'l', K specifies the number of
//           sub-diagonals of the matrix A.
//           K must satisfy  0 .le. K.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is COMPLEX//16 array, dimension ( LDA, N)
//           Before entry with UPLO = 'U' or 'u', the leading ( k + 1)
//           by n part of the array A must contain the upper triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row
//           ( k + 1) of the array, the first super-diagonal starting at
//           position 2 in row k, and so on. The top left k by k triangle
//           of the array A is not referenced.
//           The following program segment will transfer an upper
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, J = 1, N
//                    M = K + 1 - J
//                    DO 10, I = MAX( 1, J - K), J
//                       A( M + I, J) = matrix( I, J)
//              10    CONTINUE
//              20 CONTINUE
//
//           Before entry with UPLO = 'L' or 'l', the leading ( k + 1)
//           by n part of the array A must contain the lower triangular
//           band part of the matrix of coefficients, supplied column by
//           column, with the leading diagonal of the matrix in row 1 of
//           the array, the first sub-diagonal starting at position 1 in
//           row 2, and so on. The bottom right k by k triangle of the
//           array A is not referenced.
//           The following program segment will transfer a lower
//           triangular band matrix from conventional full matrix storage
//           to band storage:
//
//                 DO 20, J = 1, N
//                    M = 1 - J
//                    DO 10, I = J, MIN( N, J + K)
//                       A( M + I, J) = matrix( I, J)
//              10    CONTINUE
//              20 CONTINUE
//
//           Note that when DIAG = 'U' or 'u' the elements of the array A
//           corresponding to the diagonal elements of the matrix are not
//           referenced, but are assumed to be unity.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           ( k + 1).
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is COMPLEX//16 array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)).
//           Before entry, the incremented array X must contain the n
//           element right-hand side vector b. On exit, X is overwritten
//           with the solution vector x.
// \endverbatim
//
// \param[in] incx
// \verbatim
//          incx is INTEGER
//           On entry, incx specifies the increment for the elements of
//           X. incx must not be zero.
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
// \ingroup complex16_blas_level2
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  Level 2 Blas routine.
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Ztbsv(uplo *byte, trans *byte, diag *byte, n *int, k *int, a *[][]complex128, lda *int, x *[]complex128, incx *int) {
	zero := new(complex128)
	temp := new(complex128)
	i := new(int)
	info := new(int)
	ix := new(int)
	j := new(int)
	jx := new(int)
	kplus1 := new(int)
	kx := new(int)
	l := new(int)
	noconj := new(bool)
	nounit := new(bool)
	//*
	//*  -- Reference BLAS level2 routine (version 3.7.0) --
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
	//*     .. Parameters ..
	(*zero) = (0.0e+0 + (0.0e+0)*1i)
	//*     ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. External Functions ..
	//*     ..
	//*     .. External Subroutines ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	//*
	//*     Test the input parameters.
	//*
	(*info) = 0
	if !Lsame((*uplo), "u") && !Lsame((*uplo), "l") {
		(*info) = 1
	} else if !Lsame((*trans), "n") && !Lsame((*trans), "t") && !Lsame((*trans), "c") {
		(*info) = 2
	} else if !Lsame((*diag), "u") && !Lsame((*diag), "n") {
		(*info) = 3
	} else if (*n) < 0 {
		(*info) = 4
	} else if (*k) < 0 {
		(*info) = 5
	} else if (*lda) < ((*k) + 1) {
		(*info) = 7
	} else if (*incx) == 0 {
		(*info) = 9
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte {y :=[]byte("ztbsv "); return &y}(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if (*n) == 0 {
		return
	}
	//*
	(*noconj) = (*Lsame(trans, func() *byte {y := byte('t'); return &y}()))
	(*nounit) = (*Lsame(diag, func() *byte {y := byte('n'); return &y}()))
	//*
	//*     Set up the start point in X if the increment is not unity. This
	//*     will be  ( N - 1)*incx  too small for descending loops.
	//*
	if (*incx) <= 0 {
		(*kx) = 1 - ((*n)-1)*(*incx)
	} else if (*incx) != 1 {
		(*kx) = 1
	}
	//*
	//*     Start the operations. In this version the elements of A are
	//*     accessed by sequentially with one pass through A.
	//*
	if Lsame(trans, func() *byte {y := byte('n'); return &y}()) {
		//*
		//*        Form  x := inv( A)*x.
		//*
		if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
			(*kplus1) = (*k) + 1
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*x)[(*j)-1] != (*zero) {
						(*l) = (*kplus1) - (*j)
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] / (*a)[(*kplus1)-1][(*j)-1]
						}
						(*temp) = (*x)[(*j)-1]
						for (*i) = (*j) - 1; (*i) <= (MAX(1, (*j)-(*k))); (*i) += -1 {
							(*x)[(*i)-1] = (*x)[(*i)-1] - (*temp)*(*a)[(*l)+(*i)-1][(*j)-1]
							//Label10:
						}
					}
					//Label20:
				}
			} else {
				(*kx) = (*kx) + ((*n)-1)*(*incx)
				(*jx) = (*kx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*kx) = (*kx) - (*incx)
					if (*x)[(*jx)-1] != (*zero) {
						(*ix) = (*kx)
						(*l) = (*kplus1) - (*j)
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] / (*a)[(*kplus1)-1][(*j)-1]
						}
						(*temp) = (*x)[(*jx)-1]
						for (*i) = (*j) - 1; (*i) <= (MAX(1, (*j)-(*k))); (*i) += -1 {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] - (*temp)*(*a)[(*l)+(*i)-1][(*j)-1]
							(*ix) = (*ix) - (*incx)
							//Label30:
						}
					}
					(*jx) = (*jx) - (*incx)
					//Label40:
				}
			}
		} else {
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*x)[(*j)-1] != (*zero) {
						(*l) = 1 - (*j)
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] / (*a)[0][(*j)-1]
						}
						(*temp) = (*x)[(*j)-1]
						for (*i) = (*j) + 1; (*i) <= (MIN((*n), (*j)+(*k))); (*i)++ {
							(*x)[(*i)-1] = (*x)[(*i)-1] - (*temp)*(*a)[(*l)+(*i)-1][(*j)-1]
							//Label50:
						}
					}
					//Label60:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*kx) = (*kx) + (*incx)
					if (*x)[(*jx)-1] != (*zero) {
						(*ix) = (*kx)
						(*l) = 1 - (*j)
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] / (*a)[0][(*j)-1]
						}
						(*temp) = (*x)[(*jx)-1]
						for (*i) = (*j) + 1; (*i) <= (MIN((*n), (*j)+(*k))); (*i)++ {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] - (*temp)*(*a)[(*l)+(*i)-1][(*j)-1]
							(*ix) = (*ix) + (*incx)
							//Label70:
						}
					}
					(*jx) = (*jx) + (*incx)
					//Label80:
				}
			}
		}
	} else {
		//*
		//*        Form  x := inv( A**T)*x  or  x := inv( A**H)*x.
		//*
		if Lsame(uplo, func() *byte {y := byte('u'); return &y}()) {
			(*kplus1) = (*k) + 1
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*j)-1]
					(*l) = (*kplus1) - (*j)
					if *noconj {
						for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
							(*temp) = (*temp) - (*a)[(*l)+(*i)-1][(*j)-1]*(*x)[(*i)-1]
							//Label90:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[(*kplus1)-1][(*j)-1]
						}
					} else {
						for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
							(*temp) = (*temp) - DCONJG(((*a)[(*l)+(*i)-1][(*j)-1]))*(*x)[(*i)-1]
							//Label100:
						}
						if *nounit {
							(*temp) = (*temp) / DCONJG(((*a)[(*kplus1)-1][(*j)-1]))
						}
					}
					(*x)[(*j)-1] = (*temp)
					//Label110:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*kx)
					(*l) = (*kplus1) - (*j)
					if *noconj {
						for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
							(*temp) = (*temp) - (*a)[(*l)+(*i)-1][(*j)-1]*(*x)[(*ix)-1]
							(*ix) = (*ix) + (*incx)
							//Label120:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[(*kplus1)-1][(*j)-1]
						}
					} else {
						for (*i) = MAX(1, (*j)-(*k)); (*i) <= (*j)-1; (*i)++ {
							(*temp) = (*temp) - DCONJG(((*a)[(*l)+(*i)-1][(*j)-1]))*(*x)[(*ix)-1]
							(*ix) = (*ix) + (*incx)
							//Label130:
						}
						if *nounit {
							(*temp) = (*temp) / DCONJG(((*a)[(*kplus1)-1][(*j)-1]))
						}
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) + (*incx)
					if (*j) > (*k) {
						(*kx) = (*kx) + (*incx)
					}
					//Label140:
				}
			}
		} else {
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*j)-1]
					(*l) = 1 - (*j)
					if *noconj {
						for (*i) = MIN((*n), (*j)+(*k)); (*i) <= (*j)+1; (*i) += -1 {
							(*temp) = (*temp) - (*a)[(*l)+(*i)-1][(*j)-1]*(*x)[(*i)-1]
							//Label150:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[0][(*j)-1]
						}
					} else {
						for (*i) = MIN((*n), (*j)+(*k)); (*i) <= (*j)+1; (*i) += -1 {
							(*temp) = (*temp) - DCONJG(((*a)[(*l)+(*i)-1][(*j)-1]))*(*x)[(*i)-1]
							//Label160:
						}
						if *nounit {
							(*temp) = (*temp) / DCONJG(((*a)[0][(*j)-1]))
						}
					}
					(*x)[(*j)-1] = (*temp)
					//Label170:
				}
			} else {
				(*kx) = (*kx) + ((*n)-1)*(*incx)
				(*jx) = (*kx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*kx)
					(*l) = 1 - (*j)
					if *noconj {
						for (*i) = MIN((*n), (*j)+(*k)); (*i) <= (*j)+1; (*i) += -1 {
							(*temp) = (*temp) - (*a)[(*l)+(*i)-1][(*j)-1]*(*x)[(*ix)-1]
							(*ix) = (*ix) - (*incx)
							//Label180:
						}
						if *nounit {
							(*temp) = (*temp) / (*a)[0][(*j)-1]
						}
					} else {
						for (*i) = MIN((*n), (*j)+(*k)); (*i) <= (*j)+1; (*i) += -1 {
							(*temp) = (*temp) - DCONJG(((*a)[(*l)+(*i)-1][(*j)-1]))*(*x)[(*ix)-1]
							(*ix) = (*ix) - (*incx)
							//Label190:
						}
						if *nounit {
							(*temp) = (*temp) / DCONJG(((*a)[0][(*j)-1]))
						}
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) - (*incx)
					if ((*n) - (*j)) >= (*k) {
						(*kx) = (*kx) - (*incx)
					}
					//Label200:
				}
			}
		}
	}
	//*
	return
	//*
	//*     End of Ztbsv .
	//*
}
