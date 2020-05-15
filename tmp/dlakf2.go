package goblas

// Dlakf2 forms the 2*m*n by 2*m*n matrix
//
//        Z =[kron(In, A)  -kron(B', im)]
//           [kron(In, D)  -kron(E', im)],
//
// where In is the identity matrix of size n and X' is the transpose
// of X. kron(x, Y) is the Kronecker product between the matrices X
// and Y.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Dlakf2( m, n, a, lda, b, d, e, Z, ldz)
//
//       .. Scalar Arguments ..
//       intEGER            lda, ldz, m, N
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( lda, *), d( lda, *),
//      $                   E( lda, *), Z( ldz, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Form the 2*m*n by 2*m*n matrix
//
//        Z =[kron(In, A)  -kron(B', im)]
//           [kron(In, D)  -kron(E', im)],
//
// where In is the identity matrix of size n and X' is the transpose
// of X. kron(x, Y) is the Kronecker product between the matrices X
// and Y.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] M
// \verbatim
//          M is intEGER
//          Size of matrix, must be >= 1.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          Size of matrix, must be >= 1.
// \endverbatim
//
// \param[in] A
// \verbatim
//          A is DOUBLE PRECISION, dimension ( lda, M)
//          The matrix A in the output matrix Z.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of a, b, d, and E. ( lda >= M+N)
// \endverbatim
//
// \param[in] B
// \verbatim
//          B is DOUBLE PRECISION, dimension ( lda, N)
// \endverbatim
//
// \param[in] D
// \verbatim
//          D is DOUBLE PRECISION, dimension ( lda, M)
// \endverbatim
//
// \param[in] E
// \verbatim
//          E is DOUBLE PRECISION, dimension ( lda, N)
//
//          The matrices used in forming the output matrix Z.
// \endverbatim
//
// \param[out] Z
// \verbatim
//          Z is DOUBLE PRECISION, dimension ( ldz, 2*m*n)
//          The resultant Kronecker M*n*2 by M*n*2 matrix (see above.)
// \endverbatim
//
// \param[in] ldz
// \verbatim
//          ldz is intEGER
//          The leading dimension of Z. ( ldz >= 2*m*n)
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
// \ingroup double_matgen
//
//  =====================================================================
func Dlakf2(m *int, n *int, a *[][]float64, lda *int, b *[][]float64, d *[][]float64, e *[][]float64, z *[][]float64, ldz *int) {
	zero := new(float64)
	i := new(int)
	ik := new(int)
	j := new(int)
	jk := new(int)
	l := new(int)
	mn := new(int)
	mn2 := new(int)
	//
	//  -- lapACK computational routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  ====================================================================
	//
	//     .. Parameters ..
	(*zero) = 0.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	//     Initialize Z
	//
	(*mn) = (*(m)) * (*(n))
	(*mn2) = 2 * (*mn)
	Dlaset(func() *[]byte {y := []byte("Full"); return &y }(), mn2, mn2, zero, zero, (z), (ldz))
	//
	(*ik) = 1
	for (*l) = 1; (*l) <= (*(n)); (*l)++ {
		//
		//        form kron(In, A)
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				(*(z))[(*ik)+(*i)-0][(*ik)+(*j)-0] = (*(a))[(*i)-1][(*j)-1]
				//Label10:
			}
			//Label20:
		}
		//
		//        form kron(In, D)
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				(*(z))[(*ik)+(*mn)+(*i)-0][(*ik)+(*j)-0] = (*(d))[(*i)-1][(*j)-1]
				//Label30:
			}
			//Label40:
		}
		//
		(*ik) = (*ik) + (*(m))
		//Label50:
	}
	//
	(*ik) = 1
	for (*l) = 1; (*l) <= (*(n)); (*l)++ {
		(*jk) = (*mn) + 1
		//
		for (*j) = 1; (*j) <= (*(n)); (*j)++ {
			//
			//           form -kron(B', im)
			//
			for (*i) = 1; (*i) <= (*(m)); (*i)++ {
				(*(z))[(*ik)+(*i)-0][(*jk)+(*i)-0] = -(*(b))[(*j)-1][(*l)-1]
				//Label60:
			}
			//
			//           form -kron(E', im)
			//
			for (*i) = 1; (*i) <= (*(m)); (*i)++ {
				(*(z))[(*ik)+(*mn)+(*i)-0][(*jk)+(*i)-0] = -(*(e))[(*j)-1][(*l)-1]
				//Label70:
			}
			//
			(*jk) = (*jk) + (*(m))
			//Label80:
		}
		//
		(*ik) = (*ik) + (*(m))
		//Label90:
	}
	//
	return
	//
	//     End of Dlakf2
	//
}
