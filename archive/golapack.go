// Copyright (c) 2020 Matt Braunstein

package archive

/*
#cgo cfLAGS: -g -O2
#cgo ldfLAGS: -L/opt/OpenBLAS/lib -lopenblas
#include </opt/OpenBLAS/include/lapacke.h>
*/
import "C"

// Caxpy _constant times a vector plus a vector.
// results stored in cy
// ca * cx + cy
func Caxpy(n int, ca complex64, cx[]complex64, incX int, cy[]complex64, incY int) {
	var _cx *complex64
	if len(cx) > 0 {
		_cx = &cx[0]
	}
	var _cy *complex64
	if len(cy) > 0 {
		_cy = &cy[0]
	}

	C.lapACKE_caxpy(C.lapACK_ROW_MAJOR, (C.lapack_int)(n), (C.complex)(ca), (*C.complex)(_cx), (C.lapack_int)(incX), (*C.complex)(_cy), (C.lapack_int)(incY))
}

// Caxpywork _constant times a vector plus a vector.
// results stored in cy
// ca * cx + cy
func Caxpywork(n int, ca complex64, cx[]complex64, incX int, cy[]complex64, incY int) {
	var _cx *complex64
	if len(cx) > 0 {
		_cx = &cx[0]
	}
	var _cy *complex64
	if len(cy) > 0 {
		_cy = &cy[0]
	}

	C.lapACKE_caxpy_work(C.lapACK_ROW_MAJOR, (C.lapack_int)(n), (C.complex)(ca), (*C.complex)(_cx), (C.lapack_int)(incX), (*C.complex)(_cy), (C.lapack_int)(incY))
}
