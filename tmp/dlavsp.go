package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// dlavsp  performs one of the matrix-vector operations
//    x := A*x  or  x := A'*x,
// where x is an N element vector and  A is one of the factors
// from the block U*D*U' or L*D*l' factorization computed by Dsptrf.
//
// If trans = 'N', multiplies by U  or U * D  (or L  or L * D)
// If trans = 'T', multiplies by U' or d * U' (or L' or d * L')
// If trans = 'C', multiplies by U' or d * U' (or L' or d * L')
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE dlavsp( uplo, trans, diag, n, nrhs, a, ipiv, b, ldb,
//                          info)
//
//       .. Scalar Arguments ..
//       CHARACTER          diag, trans, uplo
//       intEGER            info, ldb, n, nrhs
//       ..
//       .. Array Arguments ..
//       intEGER            ipiv(*)
//       DOUBLE PRECISION   a(*), B( ldb, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// dlavsp  performs one of the matrix-vector operations
//    x := A*x  or  x := A'*x,
// where x is an N element vector and  A is one of the factors
// from the block U*D*U' or L*D*l' factorization computed by Dsptrf.
//
// If trans = 'N', multiplies by U  or U * D  (or L  or L * D)
// If trans = 'T', multiplies by U' or d * U' (or L' or d * L')
// If trans = 'C', multiplies by U' or d * U' (or L' or d * L')
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] uplo
// \verbatim
//          uplo is CHARACTER*1
//          Specifies whether the factor stored in A is upper or lower
//          triangular.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular
// \endverbatim
//
// \param[in] trans
// \verbatim
//          trans is CHARACTER*1
//          Specifies the operation to be performed:
//          = 'N':  x := A*x
//          = 'T':  x := A'*x
//          = 'C':  x := A'*x
// \endverbatim
//
// \param[in] diag
// \verbatim
//          diag is CHARACTER*1
//          Specifies whether or not the diagonal blocks are unit
//          matrices.  If the diagonal blocks are assumed to be unit,
//          then A = U or A = l, otherwise A = U*D or A = L*D.
//          = 'U':  diagonal blocks are assumed to be unit matrices.
//          = 'N':  diagonal blocks are assumed to be non-unit matrices.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of rows and columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] nrhs
// \verbatim
//          nrhs is intEGER
//          The number of right hand sides, i.e., the number of vectors
//          x to be multiplied by A.  nrhs >= 0.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
//          The block diagonal matrix D and the multipliers used to
//          obtain the factor U or l, stored as a packed triangular
//          matrix as computed by Dsptrf.
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is intEGER array, dimension (n)
//          The pivot indices from Dsptrf.
// \endverbatim
//
// \param[in,out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb,nrhs)
//          On entry, B contains nrhs vectors of length N.
//          On exit, B is overwritten with the product A * B.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of the array B.  ldb >= max(1,N).
// \endverbatim
//
// \param[out] info
// \verbatim
//          info is intEGER
//          = 0: successful exit
//          < 0: if info = -k, the k-th argument had an illegal value
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
// \ingroup double_lin
//
//  =====================================================================
func dlavsp(uplo *byte, trans *byte, diag *byte, n *int, nrhs *int, a *[]float64, ipiv *[]int, b *[][]float64, ldb *int, info *int) {
	one := new(float64)
	nounit := new(bool)
	j := new(int)
	k := new(int)
	kc := new(int)
	kcnext := new(int)
	kp := new(int)
	d11 := new(float64)
	d12 := new(float64)
	d21 := new(float64)
	d22 := new(float64)
	t1 := new(float64)
	t2 := new(float64)
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*one) = 1.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Test the input parameters.
	//
	(*(info)) = 0
	if !blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) && !blas.Lsame((uplo), func() *byte {y := byte('L'); return &y }()) {
		(*(info)) = -1
	} else if !blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) && !blas.Lsame((trans), func() *byte {y := byte('T'); return &y }()) && !blas.Lsame((trans), func() *byte {y := byte('C'); return &y }()) {
		(*(info)) = -2
	} else if !blas.Lsame((diag), func() *byte {y := byte('U'); return &y }()) && !blas.Lsame((diag), func() *byte {y := byte('N'); return &y }()) {
		(*(info)) = -3
	} else if (*(n)) < 0 {
		(*(info)) = -4
	} else if (*(ldb)) < (MAX(1, (*(n)))) {
		(*(info)) = -8
	}
	if (*(info)) != 0 {
		Xerbla(func() *[]byte {y := []byte("dlavsp "); return &y }(), -(*(info)))
		return
	}
	//
	//     Quick return if possible.
	//
	if (*(n)) == 0 {
		return
	}
	//
	(*nounit) = (*blas.Lsame((diag), func() *byte {y := byte('N'); return &y }()))
	//------------------------------------------
	//
	//     Compute  b := A * B  (No transpose)
	//
	//------------------------------------------
	if blas.Lsame((trans), func() *byte {y := byte('N'); return &y }()) {
		//
		//        Compute  b := U*B
		//        where U = P(m)*inv(U(m))* ... *P1*inv(U1)
		//
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			//
			//        Loop forward applying the transformations.
			//
			(*k) = 1
			(*kc) = 1
		Label10:
			;
			if (*k) > (*(n)) {
				goto Label30
			}
			//
			//           1 x 1 pivot block
			//
			if (*(ipiv))[(*k)-1] > 0 {
				//
				//              Multiply by the diagonal element if forming U * D.
				//
				if *nounit {
					Dscal((nrhs), &((*(a))[(*kc)+(*k)-0]), &((*(b))[(*k)-1][0]), (ldb))
				}
				//
				//              Multiply by P(k) * inv(U(k))  if K > 1.
				//
				if (*k) > 1 {
					//
					//                 Apply the transformation.
					//
					Dger((*k)-1, (nrhs), one, &((*(a))[(*kc)-1]), func() *int {y := 1; return &y }(), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[0][0]), (ldb))
					//
					//                 Interchange if P(k) != I.
					//
					(*kp) = (*(ipiv))[(*k)-1]
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
				}
				(*kc) = (*kc) + (*k)
				(*k) = (*k) + 1
			} else {
				//
				//              2 x 2 pivot block
				//
				(*kcnext) = (*kc) + (*k)
				//
				//              Multiply by the diagonal block if forming U * D.
				//
				if *nounit {
					(*d11) = (*(a))[(*kcnext)-0]
					(*d22) = (*(a))[(*kcnext)+(*k)-1]
					(*d12) = (*(a))[(*kcnext)+(*k)-0]
					(*d21) = (*d12)
					for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
						(*t1) = (*(b))[(*k)-1][(*j)-1]
						(*t2) = (*(b))[(*k)+0][(*j)-1]
						(*(b))[(*k)-1][(*j)-1] = (*d11)*(*t1) + (*d12)*(*t2)
						(*(b))[(*k)+0][(*j)-1] = (*d21)*(*t1) + (*d22)*(*t2)
						//Label20:
					}
				}
				//
				//              Multiply by  P(k) * inv(U(k))  if K > 1.
				//
				if (*k) > 1 {
					//
					//                 Apply the transformations.
					//
					Dger((*k)-1, (nrhs), one, &((*(a))[(*kc)-1]), func() *int {y := 1; return &y }(), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[0][0]), (ldb))
					Dger((*k)-1, (nrhs), one, &((*(a))[(*kcnext)-1]), func() *int {y := 1; return &y }(), &((*(b))[(*k)+0][0]), (ldb), &((*(b))[0][0]), (ldb))
					//
					//                 Interchange if P(k) != I.
					//
					(*kp) = (ABS(((*(ipiv))[(*k)-1])))
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
				}
				(*kc) = (*kcnext) + (*k) + 1
				(*k) = (*k) + 2
			}
			goto Label10
		Label30:

			//
			//        Compute  b := L*B
			//        where L = P1*inv(L1)* ... *P(m)*inv(L(m)) .
			//
		} else {
			//
			//           Loop backward applying the transformations to B.
			//
			(*k) = (*(n))
			(*kc) = (*(n))*((*(n))+1)/2 + 1
		Label40:
			;
			if (*k) < 1 {
				goto Label60
			}
			(*kc) = (*kc) - ((*(n)) - (*k) + 1)
			//
			//           Test the pivot index.  If greater than zero, a 1 x 1
			//           pivot was used, otherwise a 2 x 2 pivot was used.
			//
			if (*(ipiv))[(*k)-1] > 0 {
				//
				//              1 x 1 pivot block:
				//
				//              Multiply by the diagonal element if forming L * D.
				//
				if *nounit {
					Dscal((nrhs), &((*(a))[(*kc)-1]), &((*(b))[(*k)-1][0]), (ldb))
				}
				//
				//              Multiply by  P(k) * inv(L(k))  if K < N.
				//
				if (*k) != (*(n)) {
					(*kp) = (*(ipiv))[(*k)-1]
					//
					//                 Apply the transformation.
					//
					Dger((*(n))-(*k), (nrhs), one, &((*(a))[(*kc)+0]), func() *int {y := 1; return &y }(), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*k)+0][0]), (ldb))
					//
					//                 Interchange if a permutation was applied at the
					//                 K-th step of the factorization.
					//
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
				}
				(*k) = (*k) - 1
				//
			} else {
				//
				//              2 x 2 pivot block:
				//
				(*kcnext) = (*kc) - ((*(n)) - (*k) + 2)
				//
				//              Multiply by the diagonal block if forming L * D.
				//
				if *nounit {
					(*d11) = (*(a))[(*kcnext)-1]
					(*d22) = (*(a))[(*kc)-1]
					(*d21) = (*(a))[(*kcnext)+0]
					(*d12) = (*d21)
					for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
						(*t1) = (*(b))[(*k)-0][(*j)-1]
						(*t2) = (*(b))[(*k)-1][(*j)-1]
						(*(b))[(*k)-0][(*j)-1] = (*d11)*(*t1) + (*d12)*(*t2)
						(*(b))[(*k)-1][(*j)-1] = (*d21)*(*t1) + (*d22)*(*t2)
						//Label50:
					}
				}
				//
				//              Multiply by  P(k) * inv(L(k))  if K < N.
				//
				if (*k) != (*(n)) {
					//
					//                 Apply the transformation.
					//
					Dger((*(n))-(*k), (nrhs), one, &((*(a))[(*kc)+0]), func() *int {y := 1; return &y }(), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*k)+0][0]), (ldb))
					Dger((*(n))-(*k), (nrhs), one, &((*(a))[(*kcnext)+1]), func() *int {y := 1; return &y }(), &((*(b))[(*k)-0][0]), (ldb), &((*(b))[(*k)+0][0]), (ldb))
					//
					//                 Interchange if a permutation was applied at the
					//                 K-th step of the factorization.
					//
					(*kp) = (ABS(((*(ipiv))[(*k)-1])))
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
				}
				(*kc) = (*kcnext)
				(*k) = (*k) - 2
			}
			goto Label40
		Label60:
		}
		//----------------------------------------
		//
		//     Compute  b := A' * B  (transpose)
		//
		//----------------------------------------
	} else {
		//
		//        Form  b := U'*B
		//        where U  = P(m)*inv(U(m))* ... *P1*inv(U1)
		//        and   U' = inv(U'1)*P1* ... *inv(U'(m))*P(m)
		//
		if blas.Lsame((uplo), func() *byte {y := byte('U'); return &y }()) {
			//
			//           Loop backward applying the transformations.
			//
			(*k) = (*(n))
			(*kc) = (*(n))*((*(n))+1)/2 + 1
		Label70:
			;
			if (*k) < 1 {
				goto Label90
			}
			(*kc) = (*kc) - (*k)
			//
			//           1 x 1 pivot block.
			//
			if (*(ipiv))[(*k)-1] > 0 {
				if (*k) > 1 {
					//
					//                 Interchange if P(k) != I.
					//
					(*kp) = (*(ipiv))[(*k)-1]
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
					//
					//                 Apply the transformation
					//
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*k)-1, (nrhs), one, (b), (ldb), &((*(a))[(*kc)-1]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)-1][0]), (ldb))
				}
				if *nounit {
					Dscal((nrhs), &((*(a))[(*kc)+(*k)-0]), &((*(b))[(*k)-1][0]), (ldb))
				}
				(*k) = (*k) - 1
				//
				//           2 x 2 pivot block.
				//
			} else {
				(*kcnext) = (*kc) - ((*k) - 1)
				if (*k) > 2 {
					//
					//                 Interchange if P(k) != I.
					//
					(*kp) = (ABS(((*(ipiv))[(*k)-1])))
					if (*kp) != (*k)-1 {
						Dswap((nrhs), &((*(b))[(*k)-0][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
					//
					//                 Apply the transformations
					//
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*k)-2, (nrhs), one, (b), (ldb), &((*(a))[(*kc)-1]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)-1][0]), (ldb))
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*k)-2, (nrhs), one, (b), (ldb), &((*(a))[(*kcnext)-1]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)-0][0]), (ldb))
				}
				//
				//              Multiply by the diagonal block if non-unit.
				//
				if *nounit {
					(*d11) = (*(a))[(*kc)-0]
					(*d22) = (*(a))[(*kc)+(*k)-0]
					(*d12) = (*(a))[(*kc)+(*k)-1]
					(*d21) = (*d12)
					for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
						(*t1) = (*(b))[(*k)-0][(*j)-1]
						(*t2) = (*(b))[(*k)-1][(*j)-1]
						(*(b))[(*k)-0][(*j)-1] = (*d11)*(*t1) + (*d12)*(*t2)
						(*(b))[(*k)-1][(*j)-1] = (*d21)*(*t1) + (*d22)*(*t2)
						//Label80:
					}
				}
				(*kc) = (*kcnext)
				(*k) = (*k) - 2
			}
			goto Label70
		Label90:

			//
			//        Form  b := L'*B
			//        where L  = P1*inv(L1)* ... *P(m)*inv(L(m))
			//        and   L' = inv(L(m))*P(m)* ... *inv(L1)*P1
			//
		} else {
			//
			//           Loop forward applying the L-transformations.
			//
			(*k) = 1
			(*kc) = 1
		Label100:
			;
			if (*k) > (*(n)) {
				goto Label120
			}
			//
			//           1 x 1 pivot block
			//
			if (*(ipiv))[(*k)-1] > 0 {
				if (*k) < (*(n)) {
					//
					//                 Interchange if P(k) != I.
					//
					(*kp) = (*(ipiv))[(*k)-1]
					if (*kp) != (*k) {
						Dswap((nrhs), &((*(b))[(*k)-1][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
					//
					//                 Apply the transformation
					//
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*(n))-(*k), (nrhs), one, &((*(b))[(*k)+0][0]), (ldb), &((*(a))[(*kc)+0]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)-1][0]), (ldb))
				}
				if *nounit {
					Dscal((nrhs), &((*(a))[(*kc)-1]), &((*(b))[(*k)-1][0]), (ldb))
				}
				(*kc) = (*kc) + (*(n)) - (*k) + 1
				(*k) = (*k) + 1
				//
				//           2 x 2 pivot block.
				//
			} else {
				(*kcnext) = (*kc) + (*(n)) - (*k) + 1
				if (*k) < (*(n))-1 {
					//
					//              Interchange if P(k) != I.
					//
					(*kp) = (ABS(((*(ipiv))[(*k)-1])))
					if (*kp) != (*k)+1 {
						Dswap((nrhs), &((*(b))[(*k)+0][0]), (ldb), &((*(b))[(*kp)-1][0]), (ldb))
					}
					//
					//                 Apply the transformation
					//
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*(n))-(*k)-1, (nrhs), one, &((*(b))[(*k)+1][0]), (ldb), &((*(a))[(*kcnext)+0]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)+0][0]), (ldb))
					Dgemv(func() *[]byte {y := []byte("Transpose"); return &y }(), (*(n))-(*k)-1, (nrhs), one, &((*(b))[(*k)+1][0]), (ldb), &((*(a))[(*kc)+1]), func() *int {y := 1; return &y }(), one, &((*(b))[(*k)-1][0]), (ldb))
				}
				//
				//              Multiply by the diagonal block if non-unit.
				//
				if *nounit {
					(*d11) = (*(a))[(*kc)-1]
					(*d22) = (*(a))[(*kcnext)-1]
					(*d21) = (*(a))[(*kc)+0]
					(*d12) = (*d21)
					for (*j) = 1; (*j) <= (*(nrhs)); (*j)++ {
						(*t1) = (*(b))[(*k)-1][(*j)-1]
						(*t2) = (*(b))[(*k)+0][(*j)-1]
						(*(b))[(*k)-1][(*j)-1] = (*d11)*(*t1) + (*d12)*(*t2)
						(*(b))[(*k)+0][(*j)-1] = (*d21)*(*t1) + (*d22)*(*t2)
						//Label110:
					}
				}
				(*kc) = (*kcnext) + ((*(n)) - (*k))
				(*k) = (*k) + 2
			}
			goto Label100
		Label120:
		}
		//
	}
	return
	//
	//     End of dlavsp
	//
}
