package goblas

import (
	"github.com/whipstein/golapack/blas"
)

// Aladhd prints header information for the driver routines test paths.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Aladhd( iounit, path)
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
// Aladhd prints header information for the driver routines test paths.
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
//                     Assen Algorithm
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
// \date December 2016
//
// \ingroup aux_lin
//
//  =====================================================================
func Aladhd(iounit *int, path *[]byte) {
	corz := new(bool)
	sord := new(bool)
	c1 := new(byte)
	c3 := new(byte)
	p2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	sym := func() *[]byte {
		arr := make([]byte, 9)
		return &arr
	}()
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  General dense matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. diagonal                        7. Last n/2 columns zero\n    2. Upper triangular                8. Random, cndnum = sqrt(0.1/eps)\n    3. Lower triangular                9. Random, cndnum = 0.1/eps\n    4. Random, cndnum = 2             10. Scaled near underflow\n    5. First column zero              11. Scaled near overflow\n    6. Last column zero\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: ABS( work1 - RPVGRW) / ( max( work1, RPVGRW) * eps)\n")
			return &y
		}(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GB"); return &y }()) {
		//
		//        GB: General band
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  General band matrices\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Matrix types:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("    1. Random, cndnum = 2              5. Random, cndnum = sqrt(0.1/eps)\n    2. First column zero               6. Random, cndnum = 0.1/eps\n    3. Last column zero                7. Scaled near underflow\n    4. Last n/2 columns zero           8. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: ABS( work1 - RPVGRW) / ( max( work1, RPVGRW) * eps)\n")
			return &y
		}(), 7)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("GT"); return &y }()) {
		//
		//        GT: General tridiagonal
		//
		WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  General tridiagonal\n"); return &y }(), (*(path)))
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte(" Matrix types (1-6 have specified condition numbers):\n    1. diagonal                        7. Random, unspecified cndnum\n    2. Random, cndnum = 2              8. First column zero\n    3. Random, cndnum = sqrt(0.1/eps)  9. Last column zero\n    4. Random, cndnum = 0.1/eps       10. Last n/2 columns zero\n    5. Scaled near underflow          11. Scaled near underflow\n    6. Scaled near overflow           12. Scaled near overflow\n")
			return &y
		}())
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Test ratios:\n"); return &y }())
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( L * U - A)  / ( N * norm(a) * eps)\n"); return &y }(), 1)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PO"); return &y }()) || Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PP"); return &y }()) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PS"); return &y }())) {
		//
		//        PO: Positive definite full
		//        PS: Positive definite full
		//        PP: Positive definite packed
		//
		if *sord {
			(*sym) = *func() *[]byte {y := []byte("Symmetric"); return &y }()
		} else {
			(*sym) = *func() *[]byte {y := []byte("Hermitian"); return &y }()
		}
		if blas.Lsame(c3, func() *byte {y := byte('O'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite matrices\n"); return &y }(), (*(path)), (*sym))
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite packed matrices\n"); return &y }(), (*(path)), (*sym))
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PB"); return &y }()) {
		//
		//        PB: Positive definite band
		//
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite band matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite band matrices\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
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
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("PT"); return &y }()) {
		//
		//        PT: Positive definite tridiagonal
		//
		if *sord {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite tridiagonal\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {y := []byte("\n %3s drivers:  %9s positive definite tridiagonal\n"); return &y }(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
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
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 4)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SY"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("SP"); return &y }())) {
		//
		//        SY: Symmetric indefinite full
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//        SP: Symmetric indefinite packed
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		if blas.Lsame(c3, func() *byte {y := byte('Y'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s drivers:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Symmetric"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s drivers:  %9s indefinite packed matrices, partial (Bunch-Kaufman) pivoting\n")
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
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
			y := []byte("\n %3s drivers:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HA"); return &y }()) {
		//
		//        HA: Hermitian
		//            Aasen algorithm
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("\n %3s drivers:  %9s indefinite matrices, 'Aasen' Algorithm\n")
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HE"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HP"); return &y }())) {
		//
		//        HE: Hermitian indefinite full
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//        HP: Hermitian indefinite packed
		//            with partial (Bunch-Kaufman) pivoting algorithm
		//
		if blas.Lsame(c3, func() *byte {y := byte('E'); return &y }()) {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s drivers:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		} else {
			WRITE((*(iounit)), *func() *[]byte {
				y := []byte("\n %3s drivers:  %9s indefinite packed matrices, partial (Bunch-Kaufman) pivoting\n")
				return &y
			}(), (*(path)), *func() *[]byte {y := []byte("Hermitian"); return &y }())
		}
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: (backward error)   / eps\n"); return &y }(), 4)
		WRITE((*(iounit)), *func() *[]byte {
			y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * (error bound))\n")
			return &y
		}(), 5)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: rcond * cndnum - 1.0\n"); return &y }(), 6)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
		//
	} else if (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HR"); return &y }())) || (Lsamen(func() *int {y := 2; return &y }(), p2, func() *[]byte {y := []byte("HK"); return &y }())) {
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
			y := []byte("\n %3s drivers:  %9s indefinite matrices, 'rook' (bounded Bunch-Kaufman) pivoting\n")
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
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( B - A * X)  / ( norm(a) * norm(x) * eps)\n"); return &y }(), 2)
		WRITE((*(iounit)), *func() *[]byte {y := []byte("   %2d: norm( X - xact)   / ( norm(xact) * cndnum * eps)\n"); return &y }(), 3)
		WRITE((*(iounit)), *func() *[]byte {y := []byte(" Messages:\n"); return &y }())
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

	// Unused by f4go :  9891 FORMAT ( / 1 X, A3, " drivers:  ", A9, " indefinite packed matrices", ", 'rook' (bounded Bunch-Kaufman) pivoting")
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
	//     PB matrix types
	//
	//
	//     SSY, SSP, CHE, CHP matrix types
	//
	//
	//     CSY, CSP matrix types
	//
	//
	//     Test ratios
	//
	//
	return
	//
	//     End of Aladhd
	//
}
