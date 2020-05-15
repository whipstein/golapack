package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Alahd prints header information for the different test paths.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Alahd( iounit, path)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        path
//       intEGER            iounit
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Alahd prints header information for the different test paths.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] iounit
// \verbatim
//          iounit is intEGER
//          The unit number to which the header information should be
//          printed.
// \endverbatim
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          The name of the path for which the header information is to
//          be printed.  Current paths are
//             _GE:  General matrices
//             _GB:  General band
//             _GT:  General Tridiagonal
//             _PO:  Symmetric or Hermitian positive definite
//             _PS:  Symmetric or Hermitian positive semi-definite
//             _PP:  Symmetric or Hermitian positive definite packed
//             _PB:  Symmetric or Hermitian positive definite band
//             _PT:  Symmetric or Hermitian positive definite tridiagonal
//             _SY:  Symmetric indefinite,
//                     with partial (Bunch-Kaufman) pivoting
//             _SR:  Symmetric indefinite,
//                     with rook (bounded Bunch-Kaufman) pivoting
//             _SK:  Symmetric indefinite,
//                     with rook (bounded Bunch-Kaufman) pivoting
//                     ( new storage format for factors:
//                       L and diagonal of D is stored in a,
//                       subdiagonal of D is stored in E)
//             _SP:  Symmetric indefinite packed,
//                     with partial (Bunch-Kaufman) pivoting
//             _HA:  (complex) Hermitian ,
//                     with Aasen Algorithm
//             _HE:  (complex) Hermitian indefinite,
//                     with partial (Bunch-Kaufman) pivoting
//             _HR:  (complex) Hermitian indefinite,
//                     with rook (bounded Bunch-Kaufman) pivoting
//             _HK:  (complex) Hermitian indefinite,
//                     with rook (bounded Bunch-Kaufman) pivoting
//                     ( new storage format for factors:
//                       L and diagonal of D is stored in a,
//                       subdiagonal of D is stored in E)
//             _HP:  (complex) Hermitian indefinite packed,
//                     with partial (Bunch-Kaufman) pivoting
//             _TR:  Triangular
//             _TP:  Triangular packed
//             _TB:  Triangular band
//             _QR:  QR (general matrices)
//             _LQ:  LQ (general matrices)
//             _QL:  QL (general matrices)
//             _RQ:  RQ (general matrices)
//             _QP:  QR with column pivoting
//             _TZ:  Trapezoidal
//             _LS:  Least Squares driver routines
//             _LU:  LU variants
//             _CH:  Cholesky variants
//             _QS:  QR variants
//             _QT:  QRT (general matrices)
//             _QX:  QRT (triangular-pentagonal matrices)
//             _ts:  QR routines for tall-skinny and short-wide matrices
//             _HH:  Householder re_construction for tall-skinny matrices
//          The first character must be one of s, d, c, or Z (C or Z only
//          if complex).
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
// \date June 2019
//
// \ingroup aux_lin
//
//  =====================================================================
func Alahd(iounit *int, path *[]byte) {
	corz := new(bool)
	sord := new(bool)
	c1 := new(byte)
	c3 := new(byte)
	p2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	eigcnm := func() *[]byte {
		arr := make([]byte, 4)
		return &arr
	}()
	subnam := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	sym := func() *[]byte {
		arr := make([]byte, 9)
		return &arr
	}()
	//
	//  -- lapACK test routine (version 3.9.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     June 2019
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	if (*(iounit)) <= 0 {
		return
	}
	(*c1) = (*(path))[0]
	(*c3) = (*(path))[2]
	(*p2)[0] = (*(path))[1]
	(*sord) = blas.Lsame(c1, func() *byte {y := byte('S'); return &y }()) || blas.Lsame(c1, func() *byte {y := byte('D'); return &y }())
	(*corz) = blas.Lsame(c1, func() *byte {y := byte('C'); return &y }()) || blas.Lsame(c1, func() *byte {y := byte('Z'); return &y }())
	if !((*sord) || (*corz)) {
		return
	}
	//
	if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GE"); return &y }()) {
		//
		//        GE: General dense
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  General dense matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        7. Last n/2 columns zero\n    2. Upper triangular                8. Random, cndnum = sqrt(0.1/eps)\n    3. Lower triangular                9. Random, cndnum = 0.1/eps\n    4. Random, cndnum = 2             10. Scaled near underflow\n    5. First column zero              11. Scaled near overflow\n    6. Last column zero\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GB"); return &y }()) {
		//
		//        GB: General band
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  General band matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. Random, cndnum = 2              5. Random, cndnum = sqrt(0.1/eps)\n    2. First column zero               6. Random, cndnum = .01/eps\n    3. Last column zero                7. Scaled near underflow\n    4. Last n/2 columns zero           8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GT"); return &y }()) {
		//
		//        GT: General tridiagonal
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  General tridiagonal\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (1-6 have specified condition numbers):\n    1. diagonal                        7. Random, unspecified cndnum\n    2. Random, cndnum = 2              8. First column zero\n    3. Random, cndnum = sqrt(0.1/eps)  9. Last column zero\n    4. Random, cndnum = 0.1/eps       10. Last n/2 columns zero\n    5. Scaled near underflow          11. Scaled near underflow\n    6. Scaled near overflow           12. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PO"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PP"); return &y }())) {
		//
		//        PO: Positive definite full
		//        PP: Positive definite packed
		//
		if *sord {
			(*sym) = *func() *[]byte {y := []byte("Symmetric"); return &y }()
		} else {
			(*sym) = *func() *[]byte {y := []byte("Hermitian"); return &y }()
		}
		if blas.Lsame(c3, func() *byte {y := byte('O'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite matrices\n"); return &y }(), (*(path)), (*sym))
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite packed matrices\n"); return &y }(), (*(path)), (*sym))
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Random, cndnum = sqrt(0.1/eps)\n    2. Random, cndnum = 2              7. Random, cndnum = 0.1/eps\n   *3. First row and column zero       8. Scaled near underflow\n   *4. Last row and column zero        9. Scaled near overflow\n   *5. Middle row and column zero\n   (* - tests error exits from %3sTRF, no test ratios are computed)\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U' * U - A) / ( N * norm(a) * eps), or\n       norm( L * L' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PS"); return &y }()) {
		//
		//        PS: Positive semi-definite full
		//
		if *sord {
			(*sym) = *func() *[]byte {y := []byte("Symmetric"); return &y }()
		} else {
			(*sym) = *func() *[]byte {y := []byte("Hermitian"); return &y }()
		}
		if (blas.Lsame(c1, func() *byte {y := byte('S'); return &y }())) || (blas.Lsame(c1, func() *byte {y := byte('C'); return &y }())) {
			(*eigcnm) = *func() *[]byte {y := []byte("1E04"); return &y }()
		} else {
			(*eigcnm) = *func() *[]byte {y := []byte("1d12"); return &y }()
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite packed matrices\n"); return &y }(), (*(path)), (*sym))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal\n    2. Random, cndnum = 2              \n   *3. Nonzero eigenvalues of: d(1:rank-1)=1 and d(rank) = 1.0/%4s\n   *4. Nonzero eigenvalues of: d1=1 and  d(2:rank) = 1.0/%4s\n   *5. Nonzero eigenvalues of: d(i) = %4s**(-(I-1)/(rank-1))  I=1:rank\n    6. Random, cndnum = sqrt(0.1/eps)\n    7. Random, cndnum = 0.1/eps\n    8. Scaled near underflow\n    9. Scaled near overflow\n   (* - Semi-definite tests)\n")
			return &y
		}(), (*eigcnm), (*eigcnm), (*eigcnm))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Difference:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   rank minus computed rank, returned by %sPSTRF\n"); return &y }(), (*c1))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratio:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   norm( P * U' * U * P' - A) / ( N * norm(a) * eps), or\n   norm( P * L * L' * P' - A) / ( N * norm(a) * eps)\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PB"); return &y }()) {
		//
		//        PB: Positive definite band
		//
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite band matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite band matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. Random, cndnum = 2              5. Random, cndnum = sqrt(0.1/eps)\n   *2. First row and column zero       6. Random, cndnum = 0.1/eps\n   *3. Last row and column zero        7. Scaled near underflow\n   *4. Middle row and column zero      8. Scaled near overflow\n   (* - tests error exits from %3sTRF, no test ratios are computed)\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U' * U - A) / ( N * norm(a) * eps), or\n       norm( L * L' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PT"); return &y }()) {
		//
		//        PT: Positive definite tridiagonal
		//
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite tridiagonal\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %9s positive definite tridiagonal\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		}
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (1-6 have specified condition numbers):\n    1. diagonal                        7. Random, unspecified cndnum\n    2. Random, cndnum = 2              8. First row and column zero\n    3. Random, cndnum = sqrt(0.1/eps)  9. Last row and column zero\n    4. Random, cndnum = 0.1/eps       10. Middle row and column zero\n    5. Scaled near underflow          11. Scaled near underflow\n    6. Scaled near overflow           12. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U'*D*U - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SY"); return &y }()) {
		//
		//        SY: Symmetric indefinite full,
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		if blas.Lsame(c3, func() *byte {y := byte('Y'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite packed matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
				return &y
			}())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        7. Random, cndnum = sqrt(0.1/eps)\n    2. Random, cndnum = 2              8. Random, cndnum = 0.1/eps\n    3. First row and column zero       9. Scaled near underflow\n    4. Last row and column zero       10. Scaled near overflow\n    5. Middle row and column zero     11. Block diagonal matrix\n    6. Last n/2 rows and columns zero\n")
				return &y
			}())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 9)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SR"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SK"); return &y }())) {
		//
		//        SR: Symmetric indefinite full,
		//            with rook (bounded Bunch-Kaufman) pivoting algorithm
		//
		//        SK: Symmetric indefinite full,
		//            with rook (bounded Bunch-Kaufman) pivoting algorithm,
		//            ( new storage format for factors:
		//              L and diagonal of D is stored in a,
		//              subdiagonal of D is stored in E)
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
			return &y
		}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
				return &y
			}())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        7. Random, cndnum = sqrt(0.1/eps)\n    2. Random, cndnum = 2              8. Random, cndnum = 0.1/eps\n    3. First row and column zero       9. Scaled near underflow\n    4. Last row and column zero       10. Scaled near overflow\n    5. Middle row and column zero     11. Block diagonal matrix\n    6. Last n/2 rows and columns zero\n")
				return &y
			}())
		}
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: ABS( Largest element in L)\n             - ( 1 / ( 1 - alpha)) + thresh\n")
			return &y
		}(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("       where alpha = ( 1 + SQRt( 17)) / 8\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: Largest 2-Norm of 2-by-2 pivots\n             - ( ( 1 + alpha) / ( 1 - alpha)) + thresh\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("       where alpha = ( 1 + SQRt( 17)) / 8\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SP"); return &y }()) {
		//
		//        SP: Symmetric indefinite packed,
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		if blas.Lsame(c3, func() *byte {y := byte('Y'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite packed matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
				return &y
			}())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("    1. diagonal                        7. Random, cndnum = sqrt(0.1/eps)\n    2. Random, cndnum = 2              8. Random, cndnum = 0.1/eps\n    3. First row and column zero       9. Scaled near underflow\n    4. Last row and column zero       10. Scaled near overflow\n    5. Middle row and column zero     11. Block diagonal matrix\n    6. Last n/2 rows and columns zero\n")
				return &y
			}())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HA"); return &y }()) {
		//
		//        HA: Hermitian,
		//            with Assen Algorithm
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  %9s indefinite matrices, partial (Bunch-Kaufman) pivoting\n")
			return &y
		}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
			return &y
		}())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 9)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HE"); return &y }()) {
		//
		//        HE: Hermitian indefinite full,
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  %9s indefinite matrices, partial (Bunch-Kaufman) pivoting\n")
			return &y
		}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
			return &y
		}())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 9)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HR"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HR"); return &y }())) {
		//
		//        HR: Hermitian indefinite full,
		//            with rook (bounded Bunch-Kaufman) pivoting algorithm
		//
		//        HK: Hermitian indefinite full,
		//            with rook (bounded Bunch-Kaufman) pivoting algorithm,
		//            ( new storage format for factors:
		//              L and diagonal of D is stored in a,
		//              subdiagonal of D is stored in E)
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
			return &y
		}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
			return &y
		}())
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: ABS( Largest element in L)\n             - ( 1 / ( 1 - alpha)) + thresh\n")
			return &y
		}(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("       where alpha = ( 1 + SQRt( 17)) / 8\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: Largest 2-Norm of 2-by-2 pivots\n             - ( ( 1 + alpha) / ( 1 - alpha)) + thresh\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("       where alpha = ( 1 + SQRt( 17)) / 8\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HP"); return &y }()) {
		//
		//        HP: Hermitian indefinite packed,
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		if blas.Lsame(c3, func() *byte {y := byte('E'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s:  %9s indefinite packed matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		}
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Last n/2 rows and columns zero\n    2. Random, cndnum = 2              7. Random, cndnum = sqrt(0.1/eps)\n    3. First row and column zero       8. Random, cndnum = 0.1/eps\n    4. Last row and column zero        9. Scaled near underflow\n    5. Middle row and column zero     10. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U*D*U' - A) / ( N * norm(a) * eps), or\n       norm( L*D*l' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TR"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TP"); return &y }())) {
		//
		//        TR: Triangular full
		//        TP: Triangular packed
		//
		if blas.Lsame(c3, func() *byte {y := byte('R'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  Triangular matrices\n"); return &y }(), (*(path)))
			(*subnam) = append(append([]byte{}, (*(path))[0]), func() []byte {y := []byte("LATRS"); return y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  Triangular packed matrices\n"); return &y }(), (*(path)))
			(*subnam) = append(append([]byte{}, (*(path))[0]), func() []byte {y := []byte("LATPS"); return y }())
		}
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types for %3s routines:\n    1. diagonal                        6. Scaled near overflow\n    2. Random, cndnum = 2              7. Identity\n    3. Random, cndnum = sqrt(0.1/eps)  8. Unit triangular, cndnum = 2\n    4. Random, cndnum = 0.1/eps        9. Unit, cndnum = sqrt(0.1/eps)\n    5. Scaled near underflow          10. Unit, cndnum = 0.1/eps\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Special types for testing %s:\n   11. Matrix elements are O1, large right hand side\n   12. First diagonal causes overflow, offdiagonal column norms < 1\n   13. First diagonal causes overflow, offdiagonal column norms > 1\n   14. Growth factor underflows, solution does not overflow\n   15. Small diagonal causes gradual overflow\n   16. one zero diagonal element\n   17. Large offdiagonals cause overflow when adding a column\n   18. Unit triangular with large right hand side\n")
			return &y
		}(), (*subnam)[1:lenTrim(subnam)-1])
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( I - A*ainv) / ( N * norm(a) * norm(ainv) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Test ratio for %s:\n   %2d: norm( s*b - A*x)  / ( norm(a) * norm(x) * eps)\n")
			return &y
		}(), (*subnam)[1:lenTrim(subnam)-1], 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TB"); return &y }()) {
		//
		//        TB: Triangular band
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  Triangular band matrices\n"); return &y }(), (*(path)))
		(*subnam) = append(append([]byte{}, (*(path))[0]), func() []byte {y := []byte("LATBS"); return y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types for %3s routines:\n    1. Random, cndnum = 2              6. Identity\n    2. Random, cndnum = sqrt(0.1/eps)  7. Unit triangular, cndnum = 2\n    3. Random, cndnum = 0.1/eps        8. Unit, cndnum = sqrt(0.1/eps)\n    4. Scaled near underflow           9. Unit, cndnum = 0.1/eps\n    5. Scaled near overflow\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Special types for testing %s:\n   10. Matrix elements are O1, large right hand side\n   11. First diagonal causes overflow, offdiagonal column norms < 1\n   12. First diagonal causes overflow, offdiagonal column norms > 1\n   13. Growth factor underflows, solution does not overflow\n   14. Small diagonal causes gradual overflow\n   15. one zero diagonal element\n   16. Large offdiagonals cause overflow when adding a column\n   17. Unit triangular with large right hand side\n")
			return &y
		}(), (*subnam)[1:lenTrim(subnam)-1])
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps), refined\n")
			return &y
		}(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Test ratio for %s:\n   %2d: norm( s*b - A*x)  / ( norm(a) * norm(x) * eps)\n")
			return &y
		}(), (*subnam)[1:lenTrim(subnam)-1], 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QR"); return &y }()) {
		//
		//        QR decomposition of rectangular matrices
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %2s factorization of general matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("QR"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        5. Random, cndnum = sqrt(0.1/eps)\n    2. Upper triangular                6. Random, cndnum = 0.1/eps\n    3. Lower triangular                7. Scaled near underflow\n    4. Random, cndnum = 2              8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - Q' * A) / ( m * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( R - Q' * A) / ( m * norm(a) * eps)        [RFPG]\n")
			return &y
		}(), 8)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q)   / ( m * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C)  / ( %c * norm(c) * eps)\n"); return &y }(), 3, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q)  / ( %c * norm(c) * eps)\n"); return &y }(), 4, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C)/ ( %c * norm(c) * eps)\n"); return &y }(), 5, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q')/ ( %c * norm(c) * eps)\n"); return &y }(), 6, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: diagonal is not non-negative\n"); return &y }(), 9)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("LQ"); return &y }()) {
		//
		//        LQ decomposition of rectangular matrices
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %2s factorization of general matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("LQ"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        5. Random, cndnum = sqrt(0.1/eps)\n    2. Upper triangular                6. Random, cndnum = 0.1/eps\n    3. Lower triangular                7. Scaled near underflow\n    4. Random, cndnum = 2              8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L - A * Q') / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q*Q')   / ( N * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C)  / ( %c * norm(c) * eps)\n"); return &y }(), 3, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q)  / ( %c * norm(c) * eps)\n"); return &y }(), 4, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C)/ ( %c * norm(c) * eps)\n"); return &y }(), 5, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q')/ ( %c * norm(c) * eps)\n"); return &y }(), 6, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QL"); return &y }()) {
		//
		//        QL decomposition of rectangular matrices
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %2s factorization of general matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("QL"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        5. Random, cndnum = sqrt(0.1/eps)\n    2. Upper triangular                6. Random, cndnum = 0.1/eps\n    3. Lower triangular                7. Scaled near underflow\n    4. Random, cndnum = 2              8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L - Q' * A) / ( m * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q)   / ( m * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C)  / ( %c * norm(c) * eps)\n"); return &y }(), 3, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q)  / ( %c * norm(c) * eps)\n"); return &y }(), 4, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C)/ ( %c * norm(c) * eps)\n"); return &y }(), 5, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q')/ ( %c * norm(c) * eps)\n"); return &y }(), 6, 'M')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("RQ"); return &y }()) {
		//
		//        RQ decomposition of rectangular matrices
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  %2s factorization of general matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("RQ"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        5. Random, cndnum = sqrt(0.1/eps)\n    2. Upper triangular                6. Random, cndnum = 0.1/eps\n    3. Lower triangular                7. Scaled near underflow\n    4. Random, cndnum = 2              8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - A * Q') / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q*Q')   / ( N * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C)  / ( %c * norm(c) * eps)\n"); return &y }(), 3, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q)  / ( %c * norm(c) * eps)\n"); return &y }(), 4, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C)/ ( %c * norm(c) * eps)\n"); return &y }(), 5, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q')/ ( %c * norm(c) * eps)\n"); return &y }(), 6, 'N')
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QP"); return &y }()) {
		//
		//        QR decomposition with column pivoting
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  QR factorization with column pivoting\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (2-6 have condition 1/eps):\n    1. zero matrix                     4. First n/2 columns fixed\n    2. one small eigenvalue            5. Last n/2 columns fixed\n    3. Geometric distribution          6. Every second column fixed\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm(svd(a) - svd(r)) / ( m * norm(svd(r)) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( A*P - Q*R)     / ( m * norm(a) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q)      / ( m * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TZ"); return &y }()) {
		//
		//        TZ:  Trapezoidal
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  RQ factorization of trapezoidal matrix\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (2-3 have condition 1/eps):\n    1. zero matrix\n    2. one small eigenvalue\n    3. Geometric distribution\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios (1-3: %cTZRZF):\n"); return &y }(), (*c1))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm(svd(a) - svd(r)) / ( m * norm(svd(r)) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( A - R*Q)       / ( m * norm(a) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q)      / ( m * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("LS"); return &y }()) {
		//
		//        LS:  Least Squares driver routines for
		//             LS, LSD, LSS, LSX and LSY.
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  Least squares driver routines\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (1-3: full rank, 4-6: rank deficient):\n    1 and 4. Normal scaling\n    2 and 5. Scaled near overflow\n    3 and 6. Scaled near underflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Test ratios:\n    (1-2: %cGELS, 3-6: %cGELSY, 7-10: %cGELSS, 11-14: %cGELSD, 15-16: %cGEtsLS)\n")
			return &y
		}(), (*c1), (*c1), (*c1), (*c1))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( B - A * X)   / ( max(m,N) * norm(a) * norm(x) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( (A*X-B)' *a) / ( max(m,N,nrhs) * norm(a) * norm(b) * eps)\n       if trans='N' and M.GE.N or trans='T' and M.LT.N, otherwise\n       check if X is in the row space of A or A' (overdetermined case)\n")
			return &y
		}(), 2)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm(svd(a)-svd(r)) / ( min(m,N) * norm(svd(r)) * eps)\n")
			return &y
		}(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( B - A * X)   / ( max(m,N) * norm(a) * norm(x) * eps)\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( (A*X-B)' *a) / ( max(m,N,nrhs) * norm(a) * norm(b) * eps)\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: Check if X is in the row space of A or A'\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("    7-10: same as 3-6    11-14: same as 3-6\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("LU"); return &y }()) {
		//
		//        LU factorization variants
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  LU factorization variants\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        7. Last n/2 columns zero\n    2. Upper triangular                8. Random, cndnum = sqrt(0.1/eps)\n    3. Lower triangular                9. Random, cndnum = 0.1/eps\n    4. Random, cndnum = 2             10. Scaled near underflow\n    5. First column zero              11. Scaled near overflow\n    6. Last column zero\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratio:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("CH"); return &y }()) {
		//
		//        Cholesky factorization variants
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  Cholesky factorization variants\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        6. Random, cndnum = sqrt(0.1/eps)\n    2. Random, cndnum = 2              7. Random, cndnum = 0.1/eps\n   *3. First row and column zero       8. Scaled near underflow\n   *4. Last row and column zero        9. Scaled near overflow\n   *5. Middle row and column zero\n   (* - tests error exits, no test ratios are computed)\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratio:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( U' * U - A) / ( N * norm(a) * eps), or\n       norm( L * L' - A) / ( N * norm(a) * eps)\n")
			return &y
		}(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QS"); return &y }()) {
		//
		//        QR factorization variants
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  QR factorization variants\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        5. Random, cndnum = sqrt(0.1/eps)\n    2. Upper triangular                6. Random, cndnum = 0.1/eps\n    3. Lower triangular                7. Scaled near underflow\n    4. Random, cndnum = 2              8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QT"); return &y }()) {
		//
		//        QRT (general matrices)
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  QRT factorization for general matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - Q'*a) / ( m * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q) / ( m * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( m * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( m * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( m * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( m * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("QX"); return &y }()) {
		//
		//        QRT (triangular-pentagonal)
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  QRT factorization for triangular-pentagonal matrices\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - Q'*a) / ( (M+N) * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q) / ( (M+N) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( (M+N) * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("TQ"); return &y }()) {
		//
		//        QRT (triangular-pentagonal)
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  LQT factorization for general matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L - A*Q') / ( (M+N) * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q*Q') / ( (M+N) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( (M+N) * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("XQ"); return &y }()) {
		//
		//        QRT (triangular-pentagonal)
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  LQT factorization for triangular-pentagonal matrices\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L - A*Q') / ( (M+N) * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q*Q') / ( (M+N) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( (M+N) * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("ts"); return &y }()) {
		//
		//        ts:  QR routines for tall-skinny and short-wide matrices
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  ts factorization for tall-skinny or short-wide matrices\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - Q'*a) / ( (M+N) * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q) / ( (M+N) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( (M+N) * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( (M+N) * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HH"); return &y }()) {
		//
		//        HH:  Householder re_construction for tall-skinny matrices
		//
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s:  Householder recostruction from tsQR factorization output \n for tall-skinny matrices.\n")
			return &y
		}(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( R - Q'*a) / ( m * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( I - Q'*Q) / ( m * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q*C - Q*C) / ( m * norm(c) * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( Q'*C - Q'*C) / ( m * norm(c) * eps)\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q - C*Q) / ( m * norm(c) * eps)\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( C*Q' - C*Q') / ( m * norm(c) * eps)\n"); return &y }(), 6)
		//
	} else {
		//
		//        Print error message if no header is available.
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s:  No header available\n"); return &y }(), (*(path)))
	}
	//
	//     First line of header
	//
	//Label9891:

	// Unused by f4go :  9891 FORMAT ( / 1 X, A3, ":  ", A9, " indefinite packed matrices", ", 'rook' (bounded Bunch-Kaufman) pivoting")
	//
	//     GE matrix types
	//
	//
	//     GB matrix types
	//
	//
	//     GT matrix types
	//
	//
	//     PT matrix types
	//
	//
	//     PO, PP matrix types
	//
	//
	//     CH matrix types
	//
	//
	//     PS matrix types
	//
	//
	//     PB matrix types
	//
	//
	//     SSY, SSR, SSP, CHE, CHR, CHP matrix types
	//
	//
	//     CSY, CSR, CSP matrix types
	//
	//
	//     QR matrix types
	//
	//
	//     QP matrix types
	//
	//
	//     TZ matrix types
	//
	//
	//     LS matrix types
	//
	//
	//     TR, TP matrix types
	//
	//
	//     TB matrix types
	//
	//
	//     Test ratios
	//
	//
	//
	return
	//
	//     End of Alahd
	//
}
