package goblas

import 

// Dgbt01 re_constructs a band matrix  A  from its L*U factorization and
// computes the residual:
//    norm(L*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dgbt01( m, n, kl, ku, a, lda, afac, ldafac, ipiv, work,
//                          resid)
//
//       .. Scalar Arguments ..
//       intEGER            kl, ku, lda, ldafac, m, N
//       DOUBLE PRECISION   resid
//       ..
//       .. Array Arguments ..
//       intEGER            ipiv(*)
//       DOUBLE PRECISION   a( lda, *), afac( ldafac, *), work(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dgbt01 re_constructs a band matrix  A  from its L*U factorization and
// computes the residual:
//    norm(L*U - A) / ( N * norm(a) * eps),
// where eps is the machine epsilon.
//
// The expression L*U - A is computed one column at a time, so A and
// afac are not modified.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          The number of rows of the matrix A.  M >= 0.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          The number of columns of the matrix A.  N >= 0.
// \endverbatim
//
// \param[in] kl
// \verbatim
//          kl is intEGER
//          The number of subdiagonals within the band of A.  kl >= 0.
// \endverbatim
//
// \param[in] ku
// \verbatim
//          ku is intEGER
//          The number of superdiagonals within the band of A.  ku >= 0.
// \endverbatim
//
// \param[in,out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda,N)
//          The original matrix A in band storage, stored in rows 1 to
//          kl+ku+1.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER.
//          The leading dimension of the array A.  lda >= max(1,kl+ku+1).
// \endverbatim
//
// \param[in] afac
// \verbatim
//          afac is DOUBLE PRECISION array, dimension (ldafac,N)
//          The factored form of the matrix A.  afac contains the banded
//          factors L and U from the L*U factorization, as computed by
//          DGBTRF.  U is stored as an upper triangular band matrix with
//          kl+ku superdiagonals in rows 1 to kl+ku+1, and the
//          multipliers used during the factorization are stored in rows
//          kl+ku+2 to 2*kl+ku+1.  See DGBTRF for further details.
// \endverbatim
//
// \param[in] ldafac
// \verbatim
//          ldafac is intEGER
//          The leading dimension of the array afac.
//          ldafac >= max(1,2*kl*ku+1).
// \endverbatim
//
// \param[in] ipiv
// \verbatim
//          ipiv is intEGER array, dimension (min(m,N))
//          The pivot indices from DGBTRF.
// \endverbatim
//
// \param[out] work
// \verbatim
//          work is DOUBLE PRECISION array, dimension (2*kl+ku+1)
// \endverbatim
//
// \param[out] resid
// \verbatim
//          resid is DOUBLE PRECISION
//          norm(L*U - A) / ( N * norm(a) * eps)
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
func Dgbt01(m *int, n *int, kl *int, ku *int, a *[][]float64, lda *int, afac *[][]float64, ldafac *int, ipiv *[]int, work *[]float64, resid *float64) {
	zero := new(float64)
	one := new(float64)
	i := new(int)
	i1 := new(int)
	i2 := new(int)
	il := new(int)
	ip := new(int)
	iw := new(int)
	j := new(int)
	jl := new(int)
	ju := new(int)
	jua := new(int)
	kd := new(int)
	lenj := new(int)
	anorm := new(float64)
	eps := new(float64)
	t := new(float64)
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
	(*zero) = 0.0e+0
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
	//     Quick exit if M = 0 or N = 0.
	//
	(*(resid)) = (*zero)
	if (*(m)) <= 0 || (*(n)) <= 0 {
		return
	}
	//
	//     Determine eps and the norm of A.
	//
	(*eps) = (*Dlamch(func() *[]byte {y :=[]byte("epsilon"); return &y }()))
	(*kd) = (*(ku)) + 1
	(*anorm) = (*zero)
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		(*i1) = (MAX((*kd)+1-(*j), 1))
		(*i2) = (Min((*kd)+(*(m))-(*j), (*(kl))+(*kd)))
		if (*i2) >= (*i1) {
			(*anorm) = (MAX((*anorm), Dasum((*i2)-(*i1)+1, &((*(a))[(*i1)-1][(*j)-1]), func() *int {y := 1; return &y }())))
		}
		//Label10:
	}
	//
	//     Compute one column at a time of L*U - A.
	//
	(*kd) = (*(kl)) + (*(ku)) + 1
	for (*j) = 1; (*j) <= (*(n)); (*j)++ {
		//
		//        Copy the J-th column of U to work.
		//
		(*ju) = (Min((*(kl))+(*(ku)), (*j)-1))
		(*jl) = (Min((*(kl)), (*(m))-(*j)))
		(*lenj) = Min((*(m)), (*j)) - (*j) + (*ju) + 1
		if (*lenj) > 0 {
			Dcopy(lenj, &((*(afac))[(*kd)-(*ju)-1][(*j)-1]), func() *int {y := 1; return &y }(), (work), func() *int {y := 1; return &y }())
			for (*i) = (*lenj) + 1; (*i) <= (*ju)+(*jl)+1; (*i)++ {
				(*(work))[(*i)-1] = (*zero)
				//Label20:
			}
			//
			//           Multiply by the unit lower triangular matrix L.  Note that L
			//           is stored as a product of transformations and permutations.
			//
			for (*i) = Min((*(m))-1, (*j)); (*i) <= (*j)-(*ju); (*i) += -1 {
				(*il) = (Min((*(kl)), (*(m))-(*i)))
				if (*il) > 0 {
					(*iw) = (*i) - (*j) + (*ju) + 1
					(*t) = (*(work))[(*iw)-1]
					Daxpy(il, t, &((*(afac))[(*kd)+0][(*i)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*iw)+0]), func() *int {y := 1; return &y }())
					(*iP) = (*(ipiv))[(*i)-1]
					if (*i) != (*iP) {
						(*iP) = (*iP) - (*j) + (*ju) + 1
						(*(work))[(*iw)-1] = (*(work))[(*iP)-1]
						(*(work))[(*iP)-1] = (*t)
					}
				}
				//Label30:
			}
			//
			//           Subtract the corresponding column of A.
			//
			(*jua) = (Min((*ju), (*(ku))))
			if (*jua)+(*jl)+1 > 0 {
				Daxpy((*jua)+(*jl)+1, -(*one), &((*(a))[(*(ku))+1-(*jua)-1][(*j)-1]), func() *int {y := 1; return &y }(), &((*(work))[(*ju)+1-(*jua)-1]), func() *int {y := 1; return &y }())
			}
			//
			//           Compute the 1-norm of the column.
			//
			(*(resid)) = (MAX((*(resid)), Dasum((*ju)+(*jl)+1, (work), func() *int {y := 1; return &y }())))
		}
		//Label40:
	}
	//
	//     Compute norm( L*U - A) / ( N * norm(a) * eps)
	//
	if (*anorm) <= (*zero) {
		if (*(resid)) != (*zero) {
			(*(resid)) = (*one) / (*eps)
		}
	} else {
		(*(resid)) = (((*(resid)) / DBLE((*(n)))) / (*anorm)) / (*eps)
	}
	//
	return
	//
	//     End of Dgbt01
	//
}
