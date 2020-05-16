package goblas

import (
	"fmt"
	"math/cmplx"
	"testing"
)

var epsf32 = float32(1.1920929e-07)
var epsf64 = 1.1920929e-07
var roguef32 = float32(-1e10)
var roguef64 = -1e10
var roguec64 = complex64(-1e10 + 1e10i)
var roguec128 = -1e10 + 1e10i
var sfac = float32(9.765625e-4)
var dfac = 9.765625e-4
var nmax = 65
var incmax = 2
var zerof32 = func() *float32 { y := float32(0.0); return &y }()
var onef32 = func() *float32 { y := float32(1.0); return &y }()
var zerof64 = func() *float64 { y := 0.0; return &y }()
var onef64 = func() *float64 { y := 1.0; return &y }()
var zeroc64 = func() *complex64 { y := complex64(0.0); return &y }()
var onec64 = func() *complex64 { y := complex64(1.0); return &y }()
var zeroc128 = func() *complex128 { y := complex128(0.0); return &y }()
var onec128 = func() *complex128 { y := complex128(1.0); return &y }()
var zeroi = func() *int { y := 0; return &y }()
var onei = func() *int { y := 1; return &y }()
var twoi = func() *int { y := 2; return &y }()
var negonei = func() *int { y := -1; return &y }()
var _C = func() *byte { y := byte('C'); return &y }()
var _L = func() *byte { y := byte('L'); return &y }()
var _N = func() *byte { y := byte('N'); return &y }()
var _R = func() *byte { y := byte('R'); return &y }()
var _T = func() *byte { y := byte('T'); return &y }()
var _U = func() *byte { y := byte('U'); return &y }()
var _X = func() *byte { y := byte(' '); return &y }()

// BLAS Level 2 Test Conditions
var ichdL2 = []byte{'U', 'N'}
var ichtL2 = []byte{'N', 'T', 'C'}
var ichuL2 = []byte{'U', 'L'}
var idimL2 = []int{0, 1, 2, 3, 5, 9}
var kbL2 = []int{0, 1, 2, 4}
var incL2 = []int{1, 2, -1, -2}
var alff32L2 = []float32{0.0, 1.0, 0.7}
var alff64L2 = []float64{0.0, 1.0, 0.7}
var alfc64L2 = []complex64{(0.0 + 0.0i), (1.0 + 0.0i), (0.7 - 0.9i)}
var alfc128L2 = []complex128{(0.0 + 0.0i), (1.0 + 0.0i), (0.7 - 0.9i)}
var betf32L2 = []float32{0.0, 1.0, 0.9}
var betf64L2 = []float64{0.0, 1.0, 0.9}
var betc64L2 = []complex64{(0.0 + 0.0i), (1.0 + 0.0i), (1.3 - 1.1i)}
var betc128L2 = []complex128{(0.0 + 0.0i), (1.0 + 0.0i), (1.3 - 1.1i)}

// BLAS Level 3 Test Conditions
var ichdL3 = []byte{'U', 'N'}
var ichsL3 = []byte{'L', 'R'}
var ichtL3 = []byte{'N', 'T', 'C'}
var ichuL3 = []byte{'U', 'L'}
var idimL3 = []int{0, 1, 2, 3, 5, 9}
var alff32L3 = []float32{0.0, 1.0, 0.7}
var alff64L3 = []float64{0.0, 1.0, 0.7}
var alfc64L3 = []complex64{(0.0 + 0.0i), (1.0 + 0.0i), (0.7 - 0.9i)}
var alfc128L3 = []complex128{(0.0 + 0.0i), (1.0 + 0.0i), (0.7 - 0.9i)}
var betf32L3 = []float32{0.0, 1.0, 1.3}
var betf64L3 = []float64{0.0, 1.0, 1.3}
var betc64L3 = []complex64{(0.0 + 0.0i), (1.0 + 0.0i), (1.3 - 1.1i)}
var betc128L3 = []complex128{(0.0 + 0.0i), (1.0 + 0.0i), (1.3 - 1.1i)}

func printTestsRun(name string, num int) {
	fmt.Printf("\t%10s: %9d tests run\n", name, num)
}

func cReturn(a *[][]complex64, aa *[]complex64) {
	ncols := len(*a)
	for y := range *a {
		for x := range (*a)[0] {
			(*aa)[y+x*ncols] = (*a)[y][x]
		}
	}
}

func cExpand(a *[]complex64, nrows, ncols *int) *[][]complex64 {
	aa := make([][]complex64, *ncols)
	for y := range aa {
		aa[y] = make([]complex64, *nrows)
		for x := range aa[y] {
			aa[y][x] = (*a)[y+x*(*ncols)]
		}
	}

	return &aa
}

func dReturn(a *[][]float64, aa *[]float64) {
	ncols := len(*a)
	for y := range *a {
		for x := range (*a)[0] {
			(*aa)[y+x*ncols] = (*a)[y][x]
		}
	}
}

func dExpand(a *[]float64, nrows, ncols *int) *[][]float64 {
	aa := make([][]float64, *ncols)
	for y := range aa {
		aa[y] = make([]float64, *nrows)
		for x := range aa[y] {
			aa[y][x] = (*a)[y+x*(*ncols)]
		}
	}

	return &aa
}

func sReturn(a *[][]float32, aa *[]float32) {
	ncols := len(*a)
	for y := range *a {
		for x := range (*a)[0] {
			(*aa)[y+x*ncols] = (*a)[y][x]
		}
	}
}

func sExpand(a *[]float32, nrows, ncols *int) *[][]float32 {
	aa := make([][]float32, *ncols)
	for y := range aa {
		aa[y] = make([]float32, *nrows)
		for x := range aa[y] {
			aa[y][x] = (*a)[y+x*(*ncols)]
		}
	}

	return &aa
}

func zReturn(a *[][]complex128, aa *[]complex128) {
	ncols := len(*a)
	for y := range *a {
		for x := range (*a)[0] {
			(*aa)[y+x*ncols] = (*a)[y][x]
		}
	}
}

func zExpand(a *[]complex128, nrows, ncols *int) *[][]complex128 {
	aa := make([][]complex128, *ncols)
	for y := range aa {
		aa[y] = make([]complex128, *nrows)
		for x := range aa[y] {
			aa[y][x] = (*a)[y+x*(*ncols)]
		}
	}

	return &aa
}

func cTest(t *testing.T, nc, len int, ccomp, ctrue, csize *[]complex64, sfac *float32, name string) {
	scomp := make([]float32, 20)
	ssize := make([]float32, 20)
	strue := make([]float32, 20)

	for i := 1; i <= len; i++ {
		scomp[2*i-0] = real((*ccomp)[i-1])
		scomp[2*i-1] = imag((*ccomp)[i-1])
		strue[2*i-0] = real((*ctrue)[i-1])
		strue[2*i-1] = imag((*ctrue)[i-1])
		ssize[2*i-0] = real((*csize)[i-1])
		ssize[2*i-1] = imag((*csize)[i-1])
	}

	sTest(t, nc, 2*len, &scomp, &strue, &ssize, sfac, name)
}

func zTest(t *testing.T, nc, len int, ccomp, ctrue, csize *[]complex128, sfac *float64, name string) {
	scomp := make([]float64, 20)
	ssize := make([]float64, 20)
	strue := make([]float64, 20)

	for i := 1; i <= len; i++ {
		scomp[2*i-0] = real((*ccomp)[i-1])
		scomp[2*i-1] = imag((*ccomp)[i-1])
		strue[2*i-0] = real((*ctrue)[i-1])
		strue[2*i-1] = imag((*ctrue)[i-1])
		ssize[2*i-0] = real((*csize)[i-1])
		ssize[2*i-1] = imag((*csize)[i-1])
	}

	dTest(t, nc, 2*len, &scomp, &strue, &ssize, sfac, name)
}

func sTest(t *testing.T, nc, len int, scomp, strue, ssize *[]float32, sfac *float32, name string) {
	for i := 0; i < len; i++ {
		if absf32((*sfac)*((*scomp)[i]-(*strue)[i])) >= absf32((*ssize)[i])*epsilonf32()+1e-9 {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i%2, nc, absf32((*sfac)*((*scomp)[i]-(*strue)[i])), absf32((*ssize)[i])*epsilonf32()+1e-9)
		}
	}
}

func dTest(t *testing.T, nc, len int, scomp, strue, ssize *[]float64, sfac *float64, name string) {
	for i := 0; i < len; i++ {
		if absf64((*sfac)*((*scomp)[i]-(*strue)[i])) >= absf64((*ssize)[i])*epsilonf64()+1e-9 {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i%2, nc, absf64((*sfac)*((*scomp)[i]-(*strue)[i])), absf64((*ssize)[i])*epsilonf64()+1e-9)
		}
	}
}

func checkByte(t *testing.T, iter int, val, valtrue *byte, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {'%c'} output, {'%c'} expected", name, iter, *val, *valtrue)
	}
}

func checkBool(t *testing.T, iter int, val, valtrue *bool, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%t} output, {%t} expected", name, iter, *val, *valtrue)
	}
}

func checkComplex64(t *testing.T, iter int, val, valtrue *complex64, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%10.8f} output, {%10.8f} expected", name, iter, *val, *valtrue)
	}
}

func checkComplex128(t *testing.T, iter int, val, valtrue *complex128, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%10.8f} output, {%10.8f} expected", name, iter, *val, *valtrue)
	}
}

func checkFloat32(t *testing.T, iter int, val, valtrue *float32, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%10.8f} output, {%10.8f} expected", name, iter, *val, *valtrue)
	}
}

func checkFloat64(t *testing.T, iter int, val, valtrue *float64, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%10.8f} output, {%10.8f} expected", name, iter, *val, *valtrue)
	}
}

func checkInt(t *testing.T, iter int, val, valtrue *int, name string) {
	if *val != *valtrue {
		t.Errorf("Test Failed: %s: iteration %d: {%d} output, {%d} expected", name, iter, *val, *valtrue)
	}
}

func checkIntArray1D(t *testing.T, iter int, rows *int, val *[]int, valtrue *[]int, name string) {
	for i := 0; i < *rows; i++ {
		if (*val)[i] != (*valtrue)[i] {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%d} output, {%d} expected", name, i, iter, (*val)[i], (*valtrue)[i])
		}
	}
}

func checkIntArray2D(t *testing.T, iter int, rows, cols *int, val, valtrue *[][]int, name string) {
	for j := 0; j < *cols; j++ {
		for i := 0; i < *rows; i++ {
			if (*val)[i][j] != (*valtrue)[i][j] {
				t.Errorf("Test Failed: %s[%d][%d]: iteration %d: {%d} output, {%d} expected", name, i, j, iter, (*val)[i][j], (*valtrue)[i][j])
			}
		}
	}
}

func checkFloat32Array1D(t *testing.T, iter int, rows *int, val, valtrue *[]float32, name string) {
	for i := 0; i < *rows; i++ {
		if (*val)[i] != (*valtrue)[i] {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, iter, (*val)[i], (*valtrue)[i])
		}
	}
}

func checkComplex64Array1D(t *testing.T, iter int, rows *int, val, valtrue *[]complex64, name string) {
	for i := 0; i < *rows; i++ {
		if (*val)[i] != (*valtrue)[i] {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, iter, (*val)[i], (*valtrue)[i])
		}
	}
}

func checkComplex128Array1D(t *testing.T, iter int, rows *int, val, valtrue *[]complex128, name string) {
	for i := 0; i < *rows; i++ {
		if (*val)[i] != (*valtrue)[i] {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, iter, (*val)[i], (*valtrue)[i])
		}
	}
}

func checkFloat64Array1D(t *testing.T, iter int, rows *int, val, valtrue *[]float64, name string) {
	for i := 0; i < *rows; i++ {
		if (*val)[i] != (*valtrue)[i] {
			t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, iter, (*val)[i], (*valtrue)[i])
		}
	}
}

func checkFloat32Array2D(t *testing.T, iter int, rows, cols *int, val, valtrue *[][]float32, name string) {
	for j := 0; j < *cols; j++ {
		for i := 0; i < *rows; i++ {
			if (*val)[i][j] != (*valtrue)[i][j] {
				t.Errorf("Test Failed: %s[%d][%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, j, iter, (*val)[i][j], (*valtrue)[i][j])
			}
		}
	}
}

func checkFloat64Array2D(t *testing.T, iter int, rows, cols *int, val, valtrue *[][]float64, name string) {
	for j := 0; j < *cols; j++ {
		for i := 0; i < *rows; i++ {
			if (*val)[i][j] != (*valtrue)[i][j] {
				t.Errorf("Test Failed: %s[%d][%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, j, iter, (*val)[i][j], (*valtrue)[i][j])
			}
		}
	}
}

func checkComplex64Array2D(t *testing.T, iter int, rows, cols *int, val, valtrue *[][]complex64, name string) {
	for j := 0; j < *cols; j++ {
		for i := 0; i < *rows; i++ {
			if (*val)[i][j] != (*valtrue)[i][j] {
				t.Errorf("Test Failed: %s[%d][%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, j, iter, (*val)[i][j], (*valtrue)[i][j])
			}
		}
	}
}

func checkComplex128Array2D(t *testing.T, iter int, rows, cols *int, val, valtrue *[][]complex128, name string) {
	for j := 0; j < *cols; j++ {
		for i := 0; i < *rows; i++ {
			if (*val)[i][j] != (*valtrue)[i][j] {
				t.Errorf("Test Failed: %s[%d][%d]: iteration %d: {%10.8f} output, {%10.8f} expected", name, i, j, iter, (*val)[i][j], (*valtrue)[i][j])
			}
		}
	}
}

func sbegTest(reset *bool, mi, ic, i *int) float32 {
	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++

	for {
		*i *= *mi
		*i -= 1000 * ((*i) / 1000)
		if *ic >= 5 {
			*ic = 0
		} else {
			break
		}
	}

	return float32((*i)-500) / 1001.0
}

func smakeGBL2(m, n *int, a *[][]float32, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, i1, i2, i3, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i1 = 1; i1 <= (*ku)+1-j; i1++ {
			(*aa)[i1+(j-1)*_lda-1] = roguef32
		}
		for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
			(*aa)[i2+(j-1)*_lda-1] = (*a)[i2+j-(*ku)-2][j-1]
		}
		for i3 = i2; i3 <= *lda; i3++ {
			(*aa)[i3+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeGEL2(m, n *int, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)
	a := func() *[][]float32 {
		arr := make([][]float32, *m)
		for u := range arr {
			arr[u] = make([]float32, *n)
		}
		return &arr
	}()

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeGE2L2(m, n *int, a, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i+j-2] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i+j-2] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-2]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeGEL3(m, n *int, a *[][]float32, aa *[]float32, lda *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			if i != j {
				if *n > 3 && j == (*n)/2 {
					(*a)[i-1][j-1] = 0.0
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef32
		}
	}
}

func smakeSBL2(uplo *byte, m, n *int, a *[][]float32, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			iend = (*kl) + 1
		} else {
			kk = 1
			ibeg = 1
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeSPL2(uplo *byte, m, n *int, a *[][]float32, aa *[]float32, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, ioff, j int

	//
	//     Generate data in array a.
	//
	ioff = 0
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}

		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
		}
	}
}

func smakeSYL2(uplo *byte, n *int, a *[][]float32, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, _lda int

	_lda = absint(*lda)
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *n; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}

		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeSYL3(uplo *byte, m, n *int, a *[][]float32, aa *[]float32, lda *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j int
	var upper bool

	upper = *uplo == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= *lda; i++ {
			if i >= ibeg && i <= iend {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			} else {
				(*aa)[i+(j-1)*(*lda)-1] = roguef32
			}
		}
	}
}

func smakeTBL2(uplo, diag *byte, m, n *int, a *[][]float32, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			if *diag == 'U' {
				iend = *kl
			} else {
				iend = (*kl) + 1
			}
		} else {
			kk = 1
			if *diag == 'U' {
				ibeg = 2
			} else {
				ibeg = 1
			}
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeTPL2(uplo, diag *byte, m, n *int, a *[][]float32, aa *[]float32, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, ioff, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j && *diag == 'U' {
				(*aa)[ioff-1] = roguef32
			}
		}
	}
}

func smakeTRL2(uplo, diag *byte, m, n *int, a *[][]float32, aa *[]float32, lda, kl, ku *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef32
		}
	}
}

func smakeTRL3(uplo, diag *byte, m, n *int, a *[][]float32, aa *[]float32, lda *int, reset *bool, transl *float32, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = sbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] += 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef32
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef32
		}
	}
}

func smmchTest(transa, transb *byte, m, n, kk *int, alpha *float32, a, b *[][]float32, beta *float32, c *[][]float32, ct, g, cc *[]float32, ldc *int, iter int, name string, t *testing.T) {
	var i, j, k int
	var err, erri float32

	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = 0.0
			(*g)[i-1] = 0.0
		}
		if (*transa != 'T' && *transa != 'C') && (*transb != 'T' && *transb != 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf32((*a)[i-1][k-1] * (*b)[k-1][j-1])
				}
			}
		} else if (*transa == 'T' || *transa == 'C') && (*transb != 'T' && *transb != 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf32((*a)[k-1][i-1] * (*b)[k-1][j-1])
				}
			}
		} else if (*transa != 'T' && *transa != 'C') && (*transb == 'T' || *transb == 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf32((*a)[i-1][k-1] * (*b)[j-1][k-1])
				}
			}
		} else if (*transa == 'T' || *transa == 'C') && (*transb == 'T' || *transb == 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf32((*a)[k-1][i-1] * (*b)[j-1][k-1])
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = absf32(*alpha)*(*g)[i-1] + absf32((*beta)*(*c)[i-1][j-1])
		}

		err = 0
		for i = 1; i <= *m; i++ {
			erri = absf32((*ct)[i-1]-(*cc)[i+(j-1)*(*ldc)-1]) / epsf32
			if (*g)[i-1] != 0.0 {
				erri /= (*g)[i-1]
			}
			err = maxf32(err, erri)
			if err*sqrtf32(epsf32) >= 1.0 {
				t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected, {%10.8e} error", name, i+(j-1)*(*ldc)-1, iter, (*cc)[i+(j-1)*(*ldc)-1], (*ct)[i-1], err*sqrtf32(epsf32))
			}
		}
	}
}

func smvchTest(trans *byte, m, n *int, alpha *float32, a *[][]float32, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int, yt, g *[]float32) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int
	_x := make([]float32, len(*x)/absint(*incx))

	for i = range _x {
		_x[i] = (*x)[i*absint(*incx)]
	}
	if *trans == 'T' || *trans == 'C' {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if *trans == 'T' || *trans == 'C' {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * _x[jx-1]
				(*g)[iy-1] += absf32((*a)[j-1][i-1] * _x[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * _x[jx-1]
				(*g)[iy-1] += absf32((*a)[i-1][j-1] * _x[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[(iy-1)*absint(*incy)]
		(*g)[iy-1] = absf32(*alpha)*(*g)[iy-1] + absf32((*beta)*(*y)[(iy-1)*absint(*incy)])
		iy += incyl
	}
}

func smvch2Test(trans *byte, m, n *int, alpha *float32, a *[][]float32, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int, yt, g *[]float32) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int

	if *trans == 'T' || *trans == 'C' {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if *trans == 'T' || *trans == 'C' {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = absf32(*alpha)*(*g)[iy-1] + absf32((*beta)*(*y)[iy-1])
		iy += incyl
	}
}

func dbegTest(reset *bool, mi, ic, i *int) float64 {
	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++

	for {
		*i *= *mi
		*i -= 1000 * ((*i) / 1000)
		if *ic >= 5 {
			*ic = 0
		} else {
			break
		}
	}

	return float64((*i)-500) / 1001.0
}

func dmakeGBL2(m, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, i1, i2, i3, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i1 = 1; i1 <= (*ku)+1-j; i1++ {
			(*aa)[i1+(j-1)*_lda-1] = roguef64
		}
		for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
			(*aa)[i2+(j-1)*_lda-1] = (*a)[i2+j-(*ku)-2][j-1]
		}
		for i3 = i2; i3 <= *lda; i3++ {
			(*aa)[i3+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeGEL2(m, n *int, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)
	a := func() *[][]float64 {
		arr := make([][]float64, *m)
		for u := range arr {
			arr[u] = make([]float64, *n)
		}
		return &arr
	}()

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeGE2L2(m, n *int, a, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i+j-2] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i+j-2] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-2]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmake2GEL2(m, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeGEL3(m, n *int, a *[][]float64, aa *[]float64, lda *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
			if i != j {
				if *n > 3 && j == (*n)/2 {
					(*a)[i-1][j-1] = 0.0
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef64
		}
	}
}

func dmakeSBL2(uplo *byte, m, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			iend = (*kl) + 1
		} else {
			kk = 1
			ibeg = 1
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeSPL2(uplo *byte, m, n *int, a *[][]float64, aa *[]float64, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, ioff, j int

	//
	//     Generate data in array a.
	//
	ioff = 0
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}

		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
		}
	}
}

func dmakeSYL2(uplo *byte, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, _lda int

	_lda = absint(*lda)
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *n; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}

		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeSYL3(uplo *byte, m, n *int, a *[][]float64, aa *[]float64, lda *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j int
	var upper bool

	upper = *uplo == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= *lda; i++ {
			if i >= ibeg && i <= iend {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			} else {
				(*aa)[i+(j-1)*(*lda)-1] = roguef64
			}
		}
	}
}

func dmakeTBL2(uplo, diag *byte, m, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			if *diag == 'U' {
				iend = *kl
			} else {
				iend = (*kl) + 1
			}
		} else {
			kk = 1
			if *diag == 'U' {
				ibeg = 2
			} else {
				ibeg = 1
			}
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeTPL2(uplo, diag *byte, m, n *int, a *[][]float64, aa *[]float64, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, ioff, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j && *diag == 'U' {
				(*aa)[ioff-1] = roguef64
			}
		}
	}
}

func dmakeTRL2(uplo, diag *byte, m, n *int, a *[][]float64, aa *[]float64, lda, kl, ku *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguef64
		}
	}
}

func dmakeTRL3(uplo, diag *byte, m, n *int, a *[][]float64, aa *[]float64, lda *int, reset *bool, transl *float64, miSbeg, icSbeg, iSbeg *int) {
	var i, ibeg, iend, j int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = dbegTest(reset, miSbeg, icSbeg, iSbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] += 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguef64
		}
	}
}

func dmmchTest(transa, transb *byte, m, n, kk *int, alpha *float64, a, b *[][]float64, beta *float64, c *[][]float64, ct, g, cc *[]float64, ldc *int, iter int, name string, t *testing.T) {
	var i, j, k int
	var err, erri float64

	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = 0.0
			(*g)[i-1] = 0.0
		}
		if (*transa != 'T' && *transa != 'C') && (*transb != 'T' && *transb != 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf64((*a)[i-1][k-1] * (*b)[k-1][j-1])
				}
			}
		} else if (*transa == 'T' || *transa == 'C') && (*transb != 'T' && *transb != 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf64((*a)[k-1][i-1] * (*b)[k-1][j-1])
				}
			}
		} else if (*transa != 'T' && *transa != 'C') && (*transb == 'T' || *transb == 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf64((*a)[i-1][k-1] * (*b)[j-1][k-1])
				}
			}
		} else if (*transa == 'T' || *transa == 'C') && (*transb == 'T' || *transb == 'C') {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf64((*a)[k-1][i-1] * (*b)[j-1][k-1])
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = absf64(*alpha)*(*g)[i-1] + absf64((*beta)*(*c)[i-1][j-1])
		}

		err = 0
		for i = 1; i <= *m; i++ {
			erri = absf64((*ct)[i-1]-(*cc)[i+(j-1)*(*ldc)-1]) / epsf64
			if (*g)[i-1] != 0.0 {
				erri /= (*g)[i-1]
			}
			err = maxf64(err, erri)
			if err*sqrtf64(epsf64) >= 1.0 {
				t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected, {%10.8e} error", name, i+(j-1)*(*ldc)-1, iter, (*cc)[i+(j-1)*(*ldc)-1], (*ct)[i-1], err*sqrtf64(epsf64))
			}
		}
	}
}

func dmvchTest(trans *byte, m, n *int, alpha *float64, a *[][]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int, yt, g *[]float64) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int
	_x := make([]float64, len(*x)/absint(*incx))

	for i = range _x {
		_x[i] = (*x)[i*absint(*incx)]
	}
	if *trans == 'T' || *trans == 'C' {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if *trans == 'T' || *trans == 'C' {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * _x[jx-1]
				(*g)[iy-1] += absf64((*a)[j-1][i-1] * _x[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * _x[jx-1]
				(*g)[iy-1] += absf64((*a)[i-1][j-1] * _x[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[(iy-1)*absint(*incy)]
		(*g)[iy-1] = absf64(*alpha)*(*g)[iy-1] + absf64((*beta)*(*y)[(iy-1)*absint(*incy)])
		iy += incyl
	}
}

func dmvch2Test(trans *byte, m, n *int, alpha *float64, a *[][]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int, yt, g *[]float64) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int

	if *trans == 'T' || *trans == 'C' {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if *trans == 'T' || *trans == 'C' {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += absf64((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += absf64((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[(iy-1)*absint(*incy)]
		(*g)[iy-1] = absf64(*alpha)*(*g)[iy-1] + absf64((*beta)*(*y)[(iy-1)*absint(*incy)])
		iy += incyl
	}
}

func cbegTest(reset *bool, mi, mj, ic, i, j *int) complex64 {
	if *reset {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*ic = 0
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++

	for {
		*i *= *mi
		*j *= *mj
		*i -= 1000 * ((*i) / 1000)
		*j -= 1000 * ((*j) / 1000)
		if *ic >= 5 {
			*ic = 0
		} else {
			break
		}
	}

	return complex(float32((*i)-500)/1001.0, float32((*j)-500)/1001.0)
}

func cmakeGBL2(m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, i1, i2, i3, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
	}

	for j = 1; j <= *n; j++ {
		for i1 = 1; i1 <= (*ku)+1-j; i1++ {
			(*aa)[i1+(j-1)*_lda-1] = roguec64
		}
		for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
			(*aa)[i2+(j-1)*_lda-1] = (*a)[i2+j-(*ku)-2][j-1]
		}
		for i3 = i2; i3 <= *lda; i3++ {
			(*aa)[i3+(j-1)*_lda-1] = roguec64
		}
	}
}

func cmakeGEL2(m, n *int, a, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i+(j-1)*(*m)-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i+(j-1)*(*m)-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+(j-1)*(*m)-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
	}
}

func cmakeGE2L2(m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
	}
}

func cmakeGEL3(m, n *int, a *[][]complex64, aa *[]complex64, lda *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			if i != j {
				if *n > 3 && j == (*n)/2 {
					(*a)[i-1][j-1] = 0.0
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
	}
}

func cmakeHBL2(uplo *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj, kk, _lda int
	var upper bool

	upper = *uplo == 'U'
	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = complex64(cmplx.Conj(complex128((*a)[i-1][j-1])))
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)

		if upper {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			iend = (*kl) + 1
		} else {
			kk = 1
			ibeg = 1
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		jj = kk + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef32)
	}
}

func cmakeHPL2(uplo *byte, m, n *int, a *[][]complex64, aa *[]complex64, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, ioff, j int
	var upper bool

	upper = *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = complex64(cmplx.Conj(complex128((*a)[i-1][j-1])))
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j {
				(*aa)[ioff-1] = complex(real((*aa)[ioff-1]), roguef32)
			}
		}
	}
}

func cmakeHEL2(uplo, diag *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = complex64(cmplx.Conj(complex128((*a)[i-1][j-1])))
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		jj = j + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef32)
	}
}

func cmakeHEL3(uplo *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj int
	upper := *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = complex64(cmplx.Conj(complex128((*a)[i-1][j-1])))
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
		jj = j + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef32)
	}
}

func cmakeSYL3(uplo *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j int
	upper := *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
	}
}

func cmakeTBL2(uplo, diag *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			if *diag == 'U' {
				iend = *kl
			} else {
				iend = (*kl) + 1
			}
		} else {
			kk = 1
			if *diag == 'U' {
				ibeg = 2
			} else {
				ibeg = 1
			}
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
	}
}

func cmakeTPL2(uplo, diag *byte, m, n *int, a *[][]complex64, aa *[]complex64, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, ioff, j int
	var upper bool

	upper = *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j && *diag == 'U' {
				(*aa)[ioff-1] = roguec64
			}
		}
	}
}

func cmakeTRL2(uplo, diag *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda, kl, ku *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, _lda int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'
	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if upper {
			ibeg = 1
			if unit {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if unit {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec64
		}
	}
}

func cmakeTRL3(uplo, diag *byte, m, n *int, a *[][]complex64, aa *[]complex64, lda *int, reset *bool, transl *complex64, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = cbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] += 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec64
		}
	}
}

func cmmchTest(transa, transb *byte, m, n, kk *int, alpha *complex64, a, b *[][]complex64, beta *complex64, c *[][]complex64, ct *[]complex64, g *[]float32, cc *[]complex64, ldc *int, iter int, name string, t *testing.T) {
	var i, j, k int
	var err, erri float32

	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = 0.0
			(*g)[i-1] = 0.0
		}
		switch {
		case (*transa != 'T' && *transa != 'C') && (*transb != 'T' && *transb != 'C'):
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abssumf32((*a)[i-1][k-1]) * abssumf32((*b)[k-1][j-1])
				}
			}
		case (*transa == 'T' || *transa == 'C') && (*transb != 'T' && *transb != 'C'):
			if *transa == 'C' {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += complex64(cmplx.Conj(complex128((*a)[k-1][i-1]))) * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[k-1][j-1])
					}
				}
			}
		case (*transa != 'T' && *transa != 'C') && (*transb == 'T' || *transb == 'C'):
			if *transb == 'C' {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * complex64(cmplx.Conj(complex128((*b)[j-1][k-1])))
						(*g)[i-1] += abssumf32((*a)[i-1][k-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf32((*a)[i-1][k-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			}
		case (*transa == 'T' || *transa == 'C') && (*transb == 'T' || *transb == 'C'):
			switch {
			case *transa == 'C' && *transb == 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += complex64(cmplx.Conj(complex128((*a)[k-1][i-1]))) * complex64(cmplx.Conj(complex128((*b)[j-1][k-1])))
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			case *transa == 'C' && *transb != 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += complex64(cmplx.Conj(complex128((*a)[k-1][i-1]))) * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			case *transa != 'C' && *transb == 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * complex64(cmplx.Conj(complex128((*b)[j-1][k-1])))
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			case *transa != 'C' && *transb != 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf32((*a)[k-1][i-1]) * abssumf32((*b)[j-1][k-1])
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumf32(*alpha)*(*g)[i-1] + abssumf32(*beta)*abssumf32((*c)[i-1][j-1])
		}

		err = 0
		for i = 1; i <= *m; i++ {
			erri = abssumf32((*ct)[i-1]-(*cc)[i+(j-1)*(*ldc)-1]) / epsf32
			if (*g)[i-1] != 0.0 {
				erri /= (*g)[i-1]
			}
			err = maxf32(err, erri)
			if err*sqrtf32(epsf32) >= 1.0 {
				t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected, {%10.8e} error", name, i+(j-1)*(*ldc)-1, iter, (*cc)[i+(j-1)*(*ldc)-1], (*ct)[i-1], err*sqrtf32(epsf32))
			}
		}
	}
}

func cmvchTest(trans *byte, m, n *int, alpha *complex64, a *[][]complex64, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int, yt *[]complex64, g *[]float32) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int
	var tran, ctran bool

	tran = *trans == 'T'
	ctran = *trans == 'C'

	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	// _x := make([]complex64, len(*x)/absint(*incx))
	// for i = range _x {
	// 	_x[i] = (*x)[i*absint(*incx)]
	// }
	// _y := make([]complex64, len(*y)/absint(*incy))
	// for i = range _y {
	// 	_y[i] = (*y)[i*absint(*incy)]
	// }
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[j-1][i-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += complex64(cmplx.Conj(complex128((*a)[j-1][i-1]))) * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[j-1][i-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[i-1][j-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumf32(*alpha)*(*g)[iy-1] + abssumf32(*beta)*abssumf32((*y)[iy-1])
		iy += incyl
	}
}

func zbegTest(reset *bool, mi, mj, ic, i, j *int) complex128 {
	if *reset {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*ic = 0
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++

	for {
		*i *= *mi
		*j *= *mj
		*i -= 1000 * ((*i) / 1000)
		*j -= 1000 * ((*j) / 1000)
		if *ic >= 5 {
			*ic = 0
		} else {
			break
		}
	}

	return complex(float64((*i)-500)/1001.0, float64((*j)-500)/1001.0)
}

func zmakeGBL2(m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, i1, i2, i3, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
	}

	for j = 1; j <= *n; j++ {
		for i1 = 1; i1 <= (*ku)+1-j; i1++ {
			(*aa)[i1+(j-1)*_lda-1] = roguec128
		}
		for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
			(*aa)[i2+(j-1)*_lda-1] = (*a)[i2+j-(*ku)-2][j-1]
		}
		for i3 = i2; i3 <= *lda; i3++ {
			(*aa)[i3+(j-1)*_lda-1] = roguec128
		}
	}
}

func zmakeGEL2(m, n *int, a, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i+(j-1)*(*m)-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i+(j-1)*(*m)-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+(j-1)*(*m)-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
	}
}

func zmakeGE2L2(m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
				(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			} else {
				(*a)[i-1][j-1] = 0.0
			}
		}
		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
	}
}

func zmakeGEL3(m, n *int, a *[][]complex128, aa *[]complex128, lda *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, j int

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
			if i != j {
				if *n > 3 && j == (*n)/2 {
					(*a)[i-1][j-1] = 0.0
				}
			}
		}

		for i = 1; i <= *m; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = (*m) + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
	}
}

func zmakeHBL2(uplo *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj, kk, _lda int
	var upper bool

	upper = *uplo == 'U'
	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = cmplx.Conj((*a)[i-1][j-1])
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)

		if upper {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			iend = (*kl) + 1
		} else {
			kk = 1
			ibeg = 1
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		jj = kk + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef64)
	}
}

func zmakeHPL2(uplo *byte, m, n *int, a *[][]complex128, aa *[]complex128, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, ioff, j int
	var upper bool

	upper = *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = cmplx.Conj((*a)[i-1][j-1])
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j {
				(*aa)[ioff-1] = complex(real((*aa)[ioff-1]), roguef64)
			}
		}
	}
}

func zmakeHEL2(uplo, diag *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = cmplx.Conj((*a)[i-1][j-1])
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		jj = j + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef64)
	}
}

func zmakeHEL3(uplo *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, jj int
	upper := *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = cmplx.Conj((*a)[i-1][j-1])
				}
			}
		}
		(*a)[j-1][j-1] = complex(real((*a)[j-1][j-1]), 0.0)
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
		jj = j + (j-1)*(*lda)
		(*aa)[jj-1] = complex(real((*aa)[jj-1]), roguef64)
	}
}

func zmakeSYL3(uplo *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j int
	upper := *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = (*a)[i-1][j-1]
				}
			}
		}
	}

	for j = 1; j <= *n; j++ {
		if upper {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
	}
}

func zmakeTBL2(uplo, diag *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, kk, _lda int

	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (*uplo == 'U' && i <= j) || (*uplo == 'L' && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			kk = (*kl) + 1
			ibeg = max(1, (*kl)+2-j)
			if *diag == 'U' {
				iend = *kl
			} else {
				iend = (*kl) + 1
			}
		} else {
			kk = 1
			if *diag == 'U' {
				ibeg = 2
			} else {
				ibeg = 1
			}
			iend = min((*kl)+1, 1+(*m)-j)
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i+j-kk-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
	}
}

func zmakeTPL2(uplo, diag *byte, m, n *int, a *[][]complex128, aa *[]complex128, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, ioff, j int
	var upper bool

	upper = *uplo == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if *diag == 'U' {
			(*a)[j-1][j-1] = 1.0
		}
	}

	ioff = 0
	for j = 1; j <= *n; j++ {
		if *uplo == 'U' {
			ibeg = 1
			iend = j
		} else {
			ibeg = j
			iend = *n
		}
		for i = ibeg; i <= iend; i++ {
			ioff++
			(*aa)[ioff-1] = (*a)[i-1][j-1]
			if i == j && *diag == 'U' {
				(*aa)[ioff-1] = roguec128
			}
		}
	}
}

func zmakeTRL2(uplo, diag *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda, kl, ku *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j, _lda int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'
	_lda = absint(*lda)

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				} else {
					(*a)[i-1][j-1] = 0.0
				}
				if i != j {
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] = (*a)[j-1][j-1] + 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if upper {
			ibeg = 1
			if unit {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if unit {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i <= ibeg-1; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*_lda-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= _lda; i++ {
			(*aa)[i+(j-1)*_lda-1] = roguec128
		}
	}
}

func zmakeTRL3(uplo, diag *byte, m, n *int, a *[][]complex128, aa *[]complex128, lda *int, reset *bool, transl *complex128, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg *int) {
	var i, ibeg, iend, j int
	var upper, unit bool

	upper = *uplo == 'U'
	unit = *diag == 'U'

	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if (upper && i <= j) || (!upper && i >= j) {
				(*a)[i-1][j-1] = zbegTest(reset, miCbeg, mjCbeg, icSbeg, iCbeg, jCbeg) + (*transl)
				if i != j {
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = 0.0
					}
					(*a)[j-1][i-1] = 0.0
				}
			}
		}
		(*a)[j-1][j-1] += 1.0
		if unit {
			(*a)[j-1][j-1] = 1.0
		}

		if *uplo == 'U' {
			ibeg = 1
			if *diag == 'U' {
				iend = j - 1
			} else {
				iend = j
			}
		} else {
			if *diag == 'U' {
				ibeg = j + 1
			} else {
				ibeg = j
			}
			iend = *n
		}
		for i = 1; i < ibeg; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
		for i = ibeg; i <= iend; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
		}
		for i = iend + 1; i <= *lda; i++ {
			(*aa)[i+(j-1)*(*lda)-1] = roguec128
		}
	}
}

func zmmchTest(transa, transb *byte, m, n, kk *int, alpha *complex128, a, b *[][]complex128, beta *complex128, c *[][]complex128, ct *[]complex128, g *[]float64, cc *[]complex128, ldc *int, iter int, name string, t *testing.T) {
	var i, j, k int
	var err, erri float64

	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = 0.0
			(*g)[i-1] = 0.0
		}
		switch {
		case (*transa != 'T' && *transa != 'C') && (*transb != 'T' && *transb != 'C'):
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[k-1][j-1])
				}
			}
		case (*transa == 'T' || *transa == 'C') && (*transb != 'T' && *transb != 'C'):
			if *transa == 'C' {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += cmplx.Conj((*a)[k-1][i-1]) * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[k-1][j-1])
					}
				}
			}
		case (*transa != 'T' && *transa != 'C') && (*transb == 'T' || *transb == 'C'):
			if *transb == 'C' {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * cmplx.Conj((*b)[j-1][k-1])
						(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			}
		case (*transa == 'T' || *transa == 'C') && (*transb == 'T' || *transb == 'C'):
			switch {
			case *transa == 'C' && *transb == 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += cmplx.Conj((*a)[k-1][i-1]) * cmplx.Conj((*b)[j-1][k-1])
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			case *transa == 'C' && *transb != 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += cmplx.Conj((*a)[k-1][i-1]) * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			case *transa != 'C' && *transb == 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * cmplx.Conj((*b)[j-1][k-1])
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			case *transa != 'C' && *transb != 'C':
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumf64(*alpha)*(*g)[i-1] + abssumf64(*beta)*abssumf64((*c)[i-1][j-1])
		}

		err = 0
		for i = 1; i <= *m; i++ {
			erri = abssumf64((*ct)[i-1]-(*cc)[i+(j-1)*(*ldc)-1]) / epsf64
			if (*g)[i-1] != 0.0 {
				erri /= (*g)[i-1]
			}
			err = maxf64(err, erri)
			if err*sqrtf64(epsf64) >= 1.0 {
				t.Errorf("Test Failed: %s[%d]: iteration %d: {%10.8f} output, {%10.8f} expected, {%10.8e} error", name, i+(j-1)*(*ldc)-1, iter, (*cc)[i+(j-1)*(*ldc)-1], (*ct)[i-1], err*sqrtf64(epsf64))
			}
		}
	}
}

func zmvchTest(trans *byte, m, n *int, alpha *complex128, a *[][]complex128, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int, yt *[]complex128, g *[]float64) {
	var i, j, ml, nl, kx, ky, jx, iy, incxl, incyl int
	var tran, ctran bool

	tran = *trans == 'T'
	ctran = *trans == 'C'

	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0.0
		(*g)[iy-1] = 0.0
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[j-1][i-1]) * abssumf64((*x)[jx-1])
				jx += incxl
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += cmplx.Conj((*a)[j-1][i-1]) * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[j-1][i-1]) * abssumf64((*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[i-1][j-1]) * abssumf64((*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumf64(*alpha)*(*g)[iy-1] + abssumf64(*beta)*abssumf64((*y)[iy-1])
		iy += incyl
	}
}
