package goblas

// \brief \b Strmv
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Strmv(UPLO,TRANS,DIAG,N,A,LDA,X,incx)
//
//       .. Scalar Arguments ..
//       INTEGER incx,LDA,N
//       CHARACTER DIAG,TRANS,UPLO
//       ..
//       .. Array Arguments ..
//       REAL A(LDA,//),X(//)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Strmv  performs one of the matrix-vector operations
//
//    x := A//x,   or   x := A////T//x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix.
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
//           On entry, TRANS specifies the operation to be performed as
//           follows:
//
//              TRANS = 'N' or 'n'   x := A//x.
//
//              TRANS = 'T' or 't'   x := A////T//x.
//
//              TRANS = 'C' or 'c'   x := A////T//x.
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
// \param[in] A
// \verbatim
//          A is REAL array, dimension ( LDA, N)
//           Before entry with  UPLO = 'U' or 'u', the leading n by n
//           upper triangular part of the array A must contain the upper
//           triangular matrix and the strictly lower triangular part of
//           A is not referenced.
//           Before entry with UPLO = 'L' or 'l', the leading n by n
//           lower triangular part of the array A must contain the lower
//           triangular matrix and the strictly upper triangular part of
//           A is not referenced.
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced either, but are assumed to be unity.
// \endverbatim
//
// \param[in] LDA
// \verbatim
//          LDA is INTEGER
//           On entry, LDA specifies the first dimension of A as declared
//           in the calling (sub) program. LDA must be at least
//           max( 1, n).
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is REAL array, dimension at least
//           ( 1 + ( n - 1)//abs( incx)).
//           Before entry, the incremented array X must contain the n
//           element vector x. On exit, X is overwritten with the
//           transformed vector x.
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
// \ingroup single_blas_level2
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  Level 2 Blas routine.
//  The vector and matrix arguments are not referenced when n = 0, or M = 0
//
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.
// \endverbatim
//
//  =====================================================================
func Strmv(uplo *byte, trans *byte, diag *byte, n *int, a *[][]float64, lda *int, x *[]float64, incx *int) {
	zero := new(float64)
	temp := new(float64)
	i := new(int)
	info := new(int)
	ix := new(int)
	j := new(int)
	jx := new(int)
	kx := new(int)
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
	(*zero) = 0.0e+0
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
	} else if (*lda) < max(func() *int { y := 1; return &y }(), n) {
		(*info) = 6
	} else if (*incx) == 0 {
		(*info) = 8
	}
	if (*info) != 0 {
		Xerbla(func() *[]byte { y := []byte("strmv "); return &y }(), info)
		return
	}
	//*
	//*     Quick return if possible.
	//*
	if (*n) == 0 {
		return
	}
	//*
	(*nounit) = (*Lsame(diag, func() *byte { y := byte('n'); return &y }()))
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
	//*     accessed sequentially with one pass through A.
	//*
	if Lsame(trans, func() *byte { y := byte('n'); return &y }()) {
		//*
		//*        Form  x := A*x.
		//*
		if Lsame(uplo, func() *byte { y := byte('u'); return &y }()) {
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*x)[(*j)-1] != (*zero) {
						(*temp) = (*x)[(*j)-1]
						for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
							(*x)[(*i)-1] = (*x)[(*i)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
							//Label10:
						}
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] * (*a)[(*j)-1][(*j)-1]
						}
					}
					//Label20:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					if (*x)[(*jx)-1] != (*zero) {
						(*temp) = (*x)[(*jx)-1]
						(*ix) = (*kx)
						for (*i) = 1; (*i) <= (*j)-1; (*i)++ {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
							(*ix) = (*ix) + (*incx)
							//Label30:
						}
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] * (*a)[(*j)-1][(*j)-1]
						}
					}
					(*jx) = (*jx) + (*incx)
					//Label40:
				}
			}
		} else {
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*x)[(*j)-1] != (*zero) {
						(*temp) = (*x)[(*j)-1]
						for (*i) = (*n); (*i) <= (*j)+1; (*i) += -1 {
							(*x)[(*i)-1] = (*x)[(*i)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
							//Label50:
						}
						if *nounit {
							(*x)[(*j)-1] = (*x)[(*j)-1] * (*a)[(*j)-1][(*j)-1]
						}
					}
					//Label60:
				}
			} else {
				(*kx) = (*kx) + ((*n)-1)*(*incx)
				(*jx) = (*kx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					if (*x)[(*jx)-1] != (*zero) {
						(*temp) = (*x)[(*jx)-1]
						(*ix) = (*kx)
						for (*i) = (*n); (*i) <= (*j)+1; (*i) += -1 {
							(*x)[(*ix)-1] = (*x)[(*ix)-1] + (*temp)*(*a)[(*i)-1][(*j)-1]
							(*ix) = (*ix) - (*incx)
							//Label70:
						}
						if *nounit {
							(*x)[(*jx)-1] = (*x)[(*jx)-1] * (*a)[(*j)-1][(*j)-1]
						}
					}
					(*jx) = (*jx) - (*incx)
					//Label80:
				}
			}
		}
	} else {
		//*
		//*        Form  x := A**T*x.
		//*
		if Lsame(uplo, func() *byte { y := byte('u'); return &y }()) {
			if (*incx) == 1 {
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*j)-1]
					if *nounit {
						(*temp) = (*temp) * (*a)[(*j)-1][(*j)-1]
					}
					for (*i) = (*j) - 1; (*i) <= 1; (*i) += -1 {
						(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*i)-1]
						//Label90:
					}
					(*x)[(*j)-1] = (*temp)
					//Label100:
				}
			} else {
				(*jx) = (*kx) + ((*n)-1)*(*incx)
				for (*j) = (*n); (*j) <= 1; (*j) += -1 {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*jx)
					if *nounit {
						(*temp) = (*temp) * (*a)[(*j)-1][(*j)-1]
					}
					for (*i) = (*j) - 1; (*i) <= 1; (*i) += -1 {
						(*ix) = (*ix) - (*incx)
						(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*ix)-1]
						//Label110:
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) - (*incx)
					//Label120:
				}
			}
		} else {
			if (*incx) == 1 {
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*j)-1]
					if *nounit {
						(*temp) = (*temp) * (*a)[(*j)-1][(*j)-1]
					}
					for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
						(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*i)-1]
						//Label130:
					}
					(*x)[(*j)-1] = (*temp)
					//Label140:
				}
			} else {
				(*jx) = (*kx)
				for (*j) = 1; (*j) <= (*n); (*j)++ {
					(*temp) = (*x)[(*jx)-1]
					(*ix) = (*jx)
					if *nounit {
						(*temp) = (*temp) * (*a)[(*j)-1][(*j)-1]
					}
					for (*i) = (*j) + 1; (*i) <= (*n); (*i)++ {
						(*ix) = (*ix) + (*incx)
						(*temp) = (*temp) + (*a)[(*i)-1][(*j)-1]*(*x)[(*ix)-1]
						//Label150:
					}
					(*x)[(*jx)-1] = (*temp)
					(*jx) = (*jx) + (*incx)
					//Label160:
				}
			}
		}
	}
	//*
	return
	//*
	//*     End of Strmv .
	//*
}
