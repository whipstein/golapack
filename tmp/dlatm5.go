package goblas

import 

// dlatm5 generates matrices involved in the Generalized Sylvester
// equation:
//
//     A * R - L * B = C
//     d * R - L * E = F
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE dlatm5( prtype, m, n, a, lda, b, ldb, c, ldc, d, LDD,
//                          e, lde, f, ldf, r, LDR, l, ldl, alpha, qblcka,
//                          qblckb)
//
//       .. Scalar Arguments ..
//       intEGER            lda, ldb, ldc, LDD, lde, ldf, ldl, LDR, m, n,
//      $                   prtype, qblcka, qblckb
//       DOUBLE PRECISION   alpha
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   a( lda, *), B( ldb, *), c( ldc, *),
//      $                   d( LDD, *), E( lde, *), F( ldf, *),
//      $                   L( ldl, *), r( LDR, *)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// dlatm5 generates matrices involved in the Generalized Sylvester
// equation:
//
//     A * R - L * B = C
//     d * R - L * E = F
//
// They also satisfy (the diagonalization condition)
//
// [I -L] ([A  -C], [D -F])[I  R] = ([A], [D])
// [I] ([B] [E])[I]   ([B] [E])
//
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] prtype
// \verbatim
//          prtype is intEGER
//          "Points" to a certain type of the matrices to generate
//          (see further details).
// \endverbatim
//
// \param[in] M
// \verbatim
//          M is intEGER
//          Specifies the order of A and D and the number of rows in
//          c, f,  R and L.
// \endverbatim
//
// \param[in] N
// \verbatim
//          N is intEGER
//          Specifies the order of B and E and the number of columns in
//          c, f, R and L.
// \endverbatim
//
// \param[out] A
// \verbatim
//          A is DOUBLE PRECISION array, dimension (lda, M).
//          On exit A M-by-M is initialized according to prtype.
// \endverbatim
//
// \param[in] lda
// \verbatim
//          lda is intEGER
//          The leading dimension of A.
// \endverbatim
//
// \param[out] B
// \verbatim
//          B is DOUBLE PRECISION array, dimension (ldb, N).
//          On exit B N-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] ldb
// \verbatim
//          ldb is intEGER
//          The leading dimension of B.
// \endverbatim
//
// \param[out] C
// \verbatim
//          C is DOUBLE PRECISION array, dimension (ldc, N).
//          On exit C M-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] ldc
// \verbatim
//          ldc is intEGER
//          The leading dimension of C.
// \endverbatim
//
// \param[out] D
// \verbatim
//          D is DOUBLE PRECISION array, dimension (LDD, M).
//          On exit D M-by-M is initialized according to prtype.
// \endverbatim
//
// \param[in] LDD
// \verbatim
//          LDD is intEGER
//          The leading dimension of D.
// \endverbatim
//
// \param[out] E
// \verbatim
//          E is DOUBLE PRECISION array, dimension (lde, N).
//          On exit E N-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] lde
// \verbatim
//          lde is intEGER
//          The leading dimension of E.
// \endverbatim
//
// \param[out] F
// \verbatim
//          F is DOUBLE PRECISION array, dimension (ldf, N).
//          On exit F M-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] ldf
// \verbatim
//          ldf is intEGER
//          The leading dimension of F.
// \endverbatim
//
// \param[out] R
// \verbatim
//          R is DOUBLE PRECISION array, dimension (LDR, N).
//          On exit R M-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] LDR
// \verbatim
//          LDR is intEGER
//          The leading dimension of R.
// \endverbatim
//
// \param[out] L
// \verbatim
//          L is DOUBLE PRECISION array, dimension (ldl, N).
//          On exit L M-by-N is initialized according to prtype.
// \endverbatim
//
// \param[in] ldl
// \verbatim
//          ldl is intEGER
//          The leading dimension of L.
// \endverbatim
//
// \param[in] alpha
// \verbatim
//          alpha is DOUBLE PRECISION
//          Parameter used in generating prtype = 1 and 5 matrices.
// \endverbatim
//
// \param[in] qblcka
// \verbatim
//          qblcka is intEGER
//          When prtype = 3, specifies the distance between 2-by-2
//          blocks on the diagonal in A. Otherwise, qblcka is not
//          referenced. qblcka > 1.
// \endverbatim
//
// \param[in] qblckb
// \verbatim
//          qblckb is intEGER
//          When prtype = 3, specifies the distance between 2-by-2
//          blocks on the diagonal in B. Otherwise, qblckb is not
//          referenced. qblckb > 1.
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
// \date June 2016
//
// \ingroup double_matgen
//
// \par Further Details:
//  =====================
//
// \verbatim
//
//  prtype = 1: A and B are Jordan blocks, D and E are identity matrices
//
//             a : if (i == j) then a(i, j) = 1.0
//                 if (j == i + 1) then a(i, j) = -1.0
//                 else a(i, j) = 0.0,            i, j = 1...M
//
//             b : if (i == j) then B(i, j) = 1.0 - alpha
//                 if (j == i + 1) then B(i, j) = 1.0
//                 else B(i, j) = 0.0,            i, j = 1...N
//
//             d : if (i == j) then d(i, j) = 1.0
//                 else d(i, j) = 0.0,            i, j = 1...M
//
//             e : if (i == j) then E(i, j) = 1.0
//                 else E(i, j) = 0.0,            i, j = 1...N
//
//             L =  R are chosen from[-10...10],
//                  which specifies the right hand sides (C, F).
//
//  prtype = 2 or 3: Triangular and/or quasi- triangular.
//
//             a : if (i <= j) then a(i, j) =[-1...1]
//                 else a(i, j) = 0.0,             i, j = 1...M
//
//                 if (prtype = 3) then
//                    a(k + 1, k + 1) = a(k, k)
//                    a(k + 1, k) =[-1...1]
//                    sign(a(k, k + 1) = -(sin(a(k + 1, k))
//                        k = 1, M - 1, qblcka
//
//             b : if (i <= j) then B(i, j) =[-1...1]
//                 else B(i, j) = 0.0,            i, j = 1...N
//
//                 if (prtype = 3) then
//                    B(k + 1, k + 1) = B(k, k)
//                    B(k + 1, k) =[-1...1]
//                    sign(B(k, k + 1) = -(sign(B(k + 1, k))
//                        k = 1, N - 1, qblckb
//
//             d : if (i <= j) then d(i, j) =[-1...1].
//                 else d(i, j) = 0.0,            i, j = 1...M
//
//
//             e : if (i <= j) then d(i, j) =[-1...1]
//                 else E(i, j) = 0.0,            i, j = 1...N
//
//                 l, R are chosen from[-10...10],
//                 which specifies the right hand sides (C, F).
//
//  prtype = 4 Full
//             a(i, j) =[-10...10]
//             d(i, j) =[-1...1]    i,j = 1...M
//             B(i, j) =[-10...10]
//             E(i, j) =[-1...1]    i,j = 1...N
//             r(i, j) =[-10...10]
//             L(i, j) =[-1...1]    i = 1..M ,j = 1...N
//
//             l, R specifies the right hand sides (C, F).
//
//  prtype = 5 special case common and/or close eigs.
// \endverbatim
//
//  =====================================================================
func dlatm5(prtype *int, m *int, n *int, a *[][]float64, lda *int, b *[][]float64, ldb *int, c *[][]float64, ldc *int, d *[][]float64, LDd *int, e *[][]float64, lde *int, f *[][]float64, ldf *int, r *[][]float64, ldr *int, l *[][]float64, ldl *int, alpha *float64, qblcka *int, qblckb *int) {
	one := new(float64)
	zero := new(float64)
	twenty := new(float64)
	half := new(float64)
	two := new(float64)
	i := new(int)
	j := new(int)
	k := new(int)
	imeps := new(float64)
	reeps := new(float64)
	//
	//  -- lapACK computational routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2016
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
	(*zero) = 0.0e+0
	(*twenty) = 2.0e+1
	(*half) = 0.5e+0
	(*two) = 2.0e+0
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*(prtype)) == 1 {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				if (*i) == (*j) {
					(*(a))[(*i)-(1)][(*j)-(1)] = (*one)
					(*(d))[(*i)-(1)][(*j)-(1)] = (*one)
				} else if (*i) == (*j)-1 {
					(*(a))[(*i)-(1)][(*j)-(1)] = -(*one)
					(*(d))[(*i)-(1)][(*j)-(1)] = (*zero)
				} else {
					(*(a))[(*i)-(1)][(*j)-(1)] = (*zero)
					(*(d))[(*i)-(1)][(*j)-(1)] = (*zero)
				}
				//Label10:
			}
			//Label20:
		}
		//
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				if (*i) == (*j) {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*one) - (*(alpha))
					(*(e))[(*i)-(1)][(*j)-(1)] = (*one)
				} else if (*i) == (*j)-1 {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*one)
					(*(e))[(*i)-(1)][(*j)-(1)] = (*zero)
				} else {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*zero)
					(*(e))[(*i)-(1)][(*j)-(1)] = (*zero)
				}
				//Label30:
			}
			//Label40:
		}
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(r))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)/(*j)))) * (*twenty)
				(*(l))[(*i)-(1)][(*j)-(1)] = (*(r))[(*i)-(1)][(*j)-(1)]
				//Label50:
			}
			//Label60:
		}
		//
	} else if (*(prtype)) == 2 || (*(prtype)) == 3 {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				if (*i) <= (*j) {
					(*(a))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)))) * (*two)
					(*(d))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*two)
				} else {
					(*(a))[(*i)-(1)][(*j)-(1)] = (*zero)
					(*(d))[(*i)-(1)][(*j)-(1)] = (*zero)
				}
				//Label70:
			}
			//Label80:
		}
		//
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				if (*i) <= (*j) {
					(*(b))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)+(*j)))) * (*two)
					(*(e))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*j)))) * (*two)
				} else {
					(*(b))[(*i)-(1)][(*j)-(1)] = (*zero)
					(*(e))[(*i)-(1)][(*j)-(1)] = (*zero)
				}
				//Label90:
			}
			//Label100:
		}
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(r))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*twenty)
				(*(l))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)+(*j)))) * (*twenty)
				//Label110:
			}
			//Label120:
		}
		//
		if (*(prtype)) == 3 {
			if (*(qblcka)) <= 1 {
				(*(qblcka)) = 2
			}
			for (*k) = 1; (*k) <= (*(m))-1; (*k) += (*(qblcka)) {
				(*(a))[(*k)+0][(*k)+0] = (*(a))[(*k)-(1)][(*k)-(1)]
				(*(a))[(*k)+0][(*k)-(1)] = -Sin(&((*(a))[(*k)-(1)][(*k)+0]))
				//Label130:
			}
			//
			if (*(qblckb)) <= 1 {
				(*(qblckb)) = 2
			}
			for (*k) = 1; (*k) <= (*(n))-1; (*k) += (*(qblckb)) {
				(*(b))[(*k)+0][(*k)+0] = (*(b))[(*k)-(1)][(*k)-(1)]
				(*(b))[(*k)+0][(*k)-(1)] = -Sin(&((*(b))[(*k)-(1)][(*k)+0]))
				//Label140:
			}
		}
		//
	} else if (*(prtype)) == 4 {
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(m)); (*j)++ {
				(*(a))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*twenty)
				(*(d))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)+(*j)))) * (*two)
				//Label150:
			}
			//Label160:
		}
		//
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(b))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)+(*j)))) * (*twenty)
				(*(e))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*two)
				//Label170:
			}
			//Label180:
		}
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(r))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*j)/(*i)))) * (*twenty)
				(*(l))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*two)
				//Label190:
			}
			//Label200:
		}
		//
	} else if (*(prtype)) >= 5 {
		(*reeps) = (*half) * (*two) * (*twenty) / (*(alpha))
		(*imeps) = ((*half) - (*two)) / (*(alpha))
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			for (*j) = 1; (*j) <= (*(n)); (*j)++ {
				(*(r))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)*(*j)))) * (*(alpha)) / (*twenty)
				(*(l))[(*i)-(1)][(*j)-(1)] = ((*half) - Sin(DBLE((*i)+(*j)))) * (*(alpha)) / (*twenty)
				//Label210:
			}
			//Label220:
		}
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			(*(d))[(*i)-(1)][(*i)-(1)] = (*one)
			//Label230:
		}
		//
		for (*i) = 1; (*i) <= (*(m)); (*i)++ {
			if (*i) <= 4 {
				(*(a))[(*i)-(1)][(*i)-(1)] = (*one)
				if (*i) > 2 {
					(*(a))[(*i)-(1)][(*i)-(1)] = (*one) + (*reeps)
				}
				if MOD((*i), int(2)) != 0 && (*i) < (*(m)) {
					(*(a))[(*i)-(1)][(*i)+0] = (*imeps)
				} else if (*i) > 1 {
					(*(a))[(*i)-(1)][(*i)-0] = -(*imeps)
				}
			} else if (*i) <= 8 {
				if (*i) <= 6 {
					(*(a))[(*i)-(1)][(*i)-(1)] = (*reeps)
				} else {
					(*(a))[(*i)-(1)][(*i)-(1)] = -(*reeps)
				}
				if MOD((*i), int(2)) != 0 && (*i) < (*(m)) {
					(*(a))[(*i)-(1)][(*i)+0] = (*one)
				} else if (*i) > 1 {
					(*(a))[(*i)-(1)][(*i)-0] = -(*one)
				}
			} else {
				(*(a))[(*i)-(1)][(*i)-(1)] = (*one)
				if MOD((*i), int(2)) != 0 && (*i) < (*(m)) {
					(*(a))[(*i)-(1)][(*i)+0] = (*imeps) * 2
				} else if (*i) > 1 {
					(*(a))[(*i)-(1)][(*i)-0] = -(*imeps) * 2
				}
			}
			//Label240:
		}
		//
		for (*i) = 1; (*i) <= (*(n)); (*i)++ {
			(*(e))[(*i)-(1)][(*i)-(1)] = (*one)
			if (*i) <= 4 {
				(*(b))[(*i)-(1)][(*i)-(1)] = -(*one)
				if (*i) > 2 {
					(*(b))[(*i)-(1)][(*i)-(1)] = (*one) - (*reeps)
				}
				if MOD((*i), int(2)) != 0 && (*i) < (*(n)) {
					(*(b))[(*i)-(1)][(*i)+0] = (*imeps)
				} else if (*i) > 1 {
					(*(b))[(*i)-(1)][(*i)-0] = -(*imeps)
				}
			} else if (*i) <= 8 {
				if (*i) <= 6 {
					(*(b))[(*i)-(1)][(*i)-(1)] = (*reeps)
				} else {
					(*(b))[(*i)-(1)][(*i)-(1)] = -(*reeps)
				}
				if MOD((*i), int(2)) != 0 && (*i) < (*(n)) {
					(*(b))[(*i)-(1)][(*i)+0] = (*one) + (*imeps)
				} else if (*i) > 1 {
					(*(b))[(*i)-(1)][(*i)-0] = -(*one) - (*imeps)
				}
			} else {
				(*(b))[(*i)-(1)][(*i)-(1)] = (*one) - (*reeps)
				if MOD((*i), int(2)) != 0 && (*i) < (*(n)) {
					(*(b))[(*i)-(1)][(*i)+0] = (*imeps) * 2
				} else if (*i) > 1 {
					(*(b))[(*i)-(1)][(*i)-0] = -(*imeps) * 2
				}
			}
			//Label250:
		}
	}
	//
	//     Compute rhs (C, F)
	//
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), one, (a), (lda), (r), (LDR), zero, (c), (ldc))
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (n), -(*one), (l), (ldl), (b), (ldb), one, (c), (ldc))
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (m), one, (d), (LDD), (r), (LDR), zero, (F), (ldf))
	Dgemm(func() *byte {y := byte('N'); return &y }(), func() *byte {y := byte('N'); return &y }(), (m), (n), (n), -(*one), (l), (ldl), (e), (lde), one, (F), (ldf))
	//
	//     End of dlatm5
	//
}
