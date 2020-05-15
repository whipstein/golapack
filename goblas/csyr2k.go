package goblas

// Csyr2k ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE Csyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
//
//       .. Scalar Arguments ..
//       COMPLEX alpha,beta
//       INTEGER k,lda,ldb,ldc,n
//       CHARACTER trans,uplo
//       ..
//       .. Array Arguments ..
//       COMPLEX a(lda,*),b(ldb,*),c(ldc,*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Csyr2k  performs one of the symmetric rank 2k operations
//
//    c := alpha*a*b**T + alpha*b*a**T + beta*c,
//
// or
//
//    c := alpha*a**T*b + alpha*b**T*a + beta*c,
//
// where  alpha and beta  are scalars,  c is an  n by n symmetric matrix
// and  a and b  are  n by k  matrices  in the  first  case  and  k by n
// matrices in the second case.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//           On  entry,   uplo  specifies  whether  the  upper  or  lower
//           triangular  part  of the  array  c  is to be  referenced  as
//           follows:
//
//              uplo = 'U' or 'u'   Only the  upper triangular part of  c
//                                  is to be referenced.
//
//              uplo = 'L' or 'l'   Only the  lower triangular part of  c
//                                  is to be referenced.
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//           On entry,  trans  specifies the operation to be performed as
//           follows:
//
//              trans = 'N' or 'n'    c := alpha*a*b**T + alpha*b*a**T +
//                                         beta*c.
//
//              trans = 'T' or 't'    c := alpha*a**T*b + alpha*b**T*a +
//                                         beta*c.
// \endverbatim
//
// \param[in] n
// \verbatim
//          n is INTEGER
//           On entry,  n specifies the order of the matrix c.  n must be
//           at least zero.
// \endverbatim
//
// \param[in] k
// \verbatim
//          k is INTEGER
//           On entry with  trans = 'N' or 'n',  k  specifies  the number
//           of  columns  of the  matrices  a and b,  and on  entry  with
//           trans = 'T' or 't',  k  specifies  the number of rows of the
//           matrices  a and b.  k must be at least zero.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is COMPLEX
//           On entry, alpha specifies the scalar alpha.
// \endverbatim
//
// \param[in] a
// \verbatim
//          a is COMPLEX array, dimension ( lda, ka ), where ka is
//           k  when  trans = 'N' or 'n',  and is  n  otherwise.
//           Before entry with  trans = 'N' or 'n',  the  leading  n by k
//           part of the array  a  must contain the matrix  a,  otherwise
//           the leading  k by n  part of the array  a  must contain  the
//           matrix a.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is INTEGER
//           On entry, lda specifies the first dimension of a as declared
//           in  the  calling  (sub)  program.   When  trans = 'N' or 'n'
//           then  lda must be at least  max( 1, n ), otherwise  lda must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] b
// \verbatim
//          b is COMPLEX array, dimension ( ldb, kb ), where kb is
//           k  when  trans = 'N' or 'n',  and is  n  otherwise.
//           Before entry with  trans = 'N' or 'n',  the  leading  n by k
//           part of the array  b  must contain the matrix  b,  otherwise
//           the leading  k by n  part of the array  b  must contain  the
//           matrix b.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is INTEGER
//           On entry, ldb specifies the first dimension of b as declared
//           in  the  calling  (sub)  program.   When  trans = 'N' or 'n'
//           then  ldb must be at least  max( 1, n ), otherwise  ldb must
//           be at least  max( 1, k ).
// \endverbatim
//
// \param[in] beta
// \verbatim
//          beta is COMPLEX
//           On entry, beta specifies the scalar beta.
// \endverbatim
//
// \param[in,out] c
// \verbatim
//          c is COMPLEX array, dimension ( ldc, n )
//           Before entry  with  uplo = 'U' or 'u',  the leading  n by n
//           upper triangular part of the array c must contain the upper
//           triangular part  of the  symmetric matrix  and the strictly
//           lower triangular part of c is not referenced.  On exit, the
//           upper triangular part of the array  c is overwritten by the
//           upper triangular part of the updated matrix.
//           Before entry  with  uplo = 'L' or 'l',  the leading  n by n
//           lower triangular part of the array c must contain the lower
//           triangular part  of the  symmetric matrix  and the strictly
//           upper triangular part of c is not referenced.  On exit, the
//           lower triangular part of the array  c is overwritten by the
//           lower triangular part of the updated matrix.
// \endverbatim
//
// \param[in] ldc
// \verbatim
//          ldc is INTEGER
//           On entry, ldc specifies the first dimension of c as declared
//           in  the  calling  (sub)  program.   ldc  must  be  at  least
//           max( 1, n ).
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG ld.
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
//  lvel 3 Blas routine.
//
//  -- Written on 8-February-1989.
//     Jack Dongarra, Argonne National lboratory.
//     Iain Duff, AERE Harwell.
//     Jeremy Du Croz, Numerical Algorithms Group ld.
//     Sven Hammarling, Numerical Algorithms Group ld.
// \endverbatim
//
//  =====================================================================
func Csyr2k(major, uplo, trans *byte, n, k *int, alpha *complex64, a *[][]complex64, lda *int, b *[][]complex64, ldb *int, beta *complex64, c *[][]complex64, ldc *int) {
	var temp1, temp2 complex64
	var i, info, j, l, nrowa int
	var upper bool
	//
	//  -- Reference BLAS level3 routine (version 3.7.0) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG ld..--
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
	} else if !Lsame(trans, func() *byte { y := byte('N'); return &y }()) && !Lsame(trans, func() *byte { y := byte('T'); return &y }()) {
		info = 3
	} else if *n < 0 {
		info = 4
	} else if *k < 0 {
		info = 5
	} else if *lda < max(1, nrowa) {
		info = 8
	} else if *ldb < max(1, nrowa) {
		info = 10
	} else if *ldc < max(1, *n) {
		info = 13
	}
	if info != 0 {
		Xerbla(func() *string { y := "Csyr2k"; return &y }(), &info)
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
					for i = 1; i <= j; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
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
					for i = j; i <= *n; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
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
		//        Form  c := alpha*a*b**T + alpha*b*a**T + c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				if *beta == 0.0 {
					for i = 1; i <= j; i++ {
						(*c)[i-1][j-1] = 0.0
					}
				} else if *beta != 1.0 {
					for i = 1; i <= j; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
				}
				for l = 1; l <= *k; l++ {
					if ((*a)[j-1][l-1] != 0.0) || ((*b)[j-1][l-1] != 0.0) {
						temp1 = (*alpha) * (*b)[j-1][l-1]
						temp2 = (*alpha) * (*a)[j-1][l-1]
						for i = 1; i <= j; i++ {
							(*c)[i-1][j-1] += (*a)[i-1][l-1]*temp1 + (*b)[i-1][l-1]*temp2
						}
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
					for i = j; i <= *n; i++ {
						(*c)[i-1][j-1] = (*beta) * (*c)[i-1][j-1]
					}
				}
				for l = 1; l <= *k; l++ {
					if ((*a)[j-1][l-1] != 0.0) || ((*b)[j-1][l-1] != 0.0) {
						temp1 = (*alpha) * (*b)[j-1][l-1]
						temp2 = (*alpha) * (*a)[j-1][l-1]
						for i = j; i <= *n; i++ {
							(*c)[i-1][j-1] += (*a)[i-1][l-1]*temp1 + (*b)[i-1][l-1]*temp2
						}
					}
				}
			}
		}
	} else {
		//
		//        Form  c := alpha*a**T*b + alpha*b**T*a + c.
		//
		if upper {
			for j = 1; j <= *n; j++ {
				for i = 1; i <= j; i++ {
					temp1 = 0.0
					temp2 = 0.0
					for l = 1; l <= *k; l++ {
						temp1 += (*a)[l-1][i-1] * (*b)[l-1][j-1]
						temp2 += (*b)[l-1][i-1] * (*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha)*temp1 + (*alpha)*temp2
					} else {
						(*c)[i-1][j-1] = (*beta)*(*c)[i-1][j-1] + (*alpha)*temp1 + (*alpha)*temp2
					}
				}
			}
		} else {
			for j = 1; j <= *n; j++ {
				for i = j; i <= *n; i++ {
					temp1 = 0.0
					temp2 = 0.0
					for l = 1; l <= *k; l++ {
						temp1 += (*a)[l-1][i-1] * (*b)[l-1][j-1]
						temp2 += (*b)[l-1][i-1] * (*a)[l-1][j-1]
					}
					if *beta == 0.0 {
						(*c)[i-1][j-1] = (*alpha)*temp1 + (*alpha)*temp2
					} else {
						(*c)[i-1][j-1] = (*beta)*(*c)[i-1][j-1] + (*alpha)*temp1 + (*alpha)*temp2
					}
				}
			}
		}
	}
}
