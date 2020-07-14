package golapack

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestAbsf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valf32 := 0.25 * float32(i*mult)
		if got, want := absf32(valf32), float32(math.Abs(float64(valf32))); got != want {
			t.Errorf("absf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestAbsint(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valint := i * mult
		if got, want := absint(valint), int(math.Abs(float64(valint))); got != want {
			t.Errorf("absint: values do not match: expected %d got %d", want, got)
		}
	}
}

func TestAbssum(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valc64 := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		valc128 := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := abssum(valc64).(float32), float32(math.Abs(float64(real(valc64)))+math.Abs(float64(imag(valc64)))); got != want {
			t.Errorf("abssum: complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := abssum(&valc64).(float32), float32(math.Abs(float64(real(valc64)))+math.Abs(float64(imag(valc64)))); got != want {
			t.Errorf("abssum: *complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := abssum(valc128).(float64), math.Abs(real(valc128))+math.Abs(imag(valc128)); got != want {
			t.Errorf("abssum: complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := abssum(&valc128).(float64), math.Abs(real(valc128))+math.Abs(imag(valc128)); got != want {
			t.Errorf("abssum: *complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestConj(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valc64 := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		valc128 := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := conj(valc64).(complex64), complex64(cmplx.Conj(complex128(valc64))); got != want {
			t.Errorf("conj: complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := conj(&valc64).(complex64), complex64(cmplx.Conj(complex128(valc64))); got != want {
			t.Errorf("conj: *complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := conj(valc128).(complex128), cmplx.Conj(valc128); got != want {
			t.Errorf("conj: complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := conj(&valc128).(complex128), cmplx.Conj(valc128); got != want {
			t.Errorf("conj: *complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestComplx(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valf32 := 0.25 * float32(i*mult)
		valf64 := 0.25 * float64(i*mult)
		if got, want := complx(valf32, valf32).(complex64), complex64(complex(float64(valf32), float64(valf32))); got != want {
			t.Errorf("complx: float32 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := complx(&valf32, &valf32).(complex64), complex64(complex(float64(valf32), float64(valf32))); got != want {
			t.Errorf("complx: *float32 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := complx(valf64, valf64).(complex128), complex(valf64, valf64); got != want {
			t.Errorf("complx: float64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := complx(&valf64, &valf64).(complex128), complex(valf64, valf64); got != want {
			t.Errorf("complx: *float64 values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestCosf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valf32 := 0.25 * float32(i*mult)
		if got, want := cosf32(valf32), float32(math.Cos(float64(valf32))); got != want {
			t.Errorf("cosf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestExpf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		valf32 := 0.25 * float32(i*mult)
		if got, want := expf32(valf32), float32(math.Exp(float64(valf32))); got != want {
			t.Errorf("expf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestIm(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valc64 := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		valc128 := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := im(valc64).(float32), imag(valc64); got != want {
			t.Errorf("im: complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := im(&valc64).(float32), imag(valc64); got != want {
			t.Errorf("im: *complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := im(valc128).(float64), imag(valc128); got != want {
			t.Errorf("im: complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := im(&valc128).(float64), imag(valc128); got != want {
			t.Errorf("im: *complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestLogf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		valf32 := 0.25 * float32(i*mult)
		if got, want := logf32(valf32), float32(math.Log(float64(valf32))); got != want {
			t.Errorf("logf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestLog10f32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		valf32 := 0.25 * float32(i*mult)
		if got, want := log10f32(valf32), float32(math.Log10(float64(valf32))); got != want {
			t.Errorf("log10f32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestMaxf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			var wantf32 float32
			multi := 1
			if i%2 == 0 {
				multi = -1
			}
			valf32i := 0.25 * float32(i*multi)
			multj := 1
			if j%2 == 0 {
				multj = -1
			}
			valf32j := 0.25 * float32(j*multj)
			if valf32i < valf32j {
				wantf32 = valf32j
			} else {
				wantf32 = valf32i
			}
			if got := maxf32(valf32i, valf32j); got != wantf32 {
				t.Errorf("maxf32: values do not match: expected %6.4f got %6.4f", wantf32, got)
			}
		}
	}
}

func TestMaxint(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			var wantint int
			multi := 1
			if i%2 == 0 {
				multi = -1
			}
			valinti := i * multi
			multj := 1
			if j%2 == 0 {
				multj = -1
			}
			valintj := j * multj
			if valinti < valintj {
				wantint = valintj
			} else {
				wantint = valinti
			}
			if got := maxint(valinti, valintj); got != wantint {
				t.Errorf("maxint: values do not match: expected %d got %d", wantint, got)
			}
		}
	}
}

func TestModint(t *testing.T) {
	for i := 1; i < 100; i++ {
		for j := 1; j < 100; j++ {
			if got, want := modint(i, j), i%j; got != want {
				t.Errorf("modint: values do not match: expected %d got %d", want, got)
			}
		}
	}
}

func TestNintf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valf32 := 0.25 * float32(i*mult)
		if got, want := nintf32(valf32), int(math.Round(float64(valf32))); got != want {
			t.Errorf("nintf32: values do not match: expected %d got %d", want, got)
		}
	}
}

func TestPowf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			valf32i := 0.25 * float32(i)
			valf32j := 0.25 * float32(j)
			if got, want := powf32(valf32i, valf32j), float32(math.Pow(float64(valf32i), float64(valf32j))); got != want {
				t.Errorf("powf32: values do not match: expected %6.4f got %6.4f", want, got)
			}
		}
	}
}

func TestPowint(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			if got, want := powint(i, j), int(math.Pow(float64(i), float64(j))); got != want {
				t.Errorf("powint: values do not match: expected %d got %d", want, got)
			}
		}
	}
}

func TestRe(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valc64 := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		valc128 := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := re(valc64).(float32), real(valc64); got != want {
			t.Errorf("re: complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := re(&valc64).(float32), real(valc64); got != want {
			t.Errorf("re: *complex64 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := re(valc128).(float64), real(valc128); got != want {
			t.Errorf("re: complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
		if got, want := re(&valc128).(float64), real(valc128); got != want {
			t.Errorf("re: *complex128 values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestSignf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			multi := 1
			if i%2 == 0 {
				multi = -1
			}
			valf32i := 0.25 * float32(i*multi)
			multj := 1
			if j%2 == 0 {
				multj = -1
			}
			valf32j := 0.25 * float32(i*multj)
			if got, want := signf32(valf32i, valf32j), float32(math.Copysign(float64(valf32i), float64(valf32j))); got != want {
				t.Errorf("signf32: values do not match: expected %6.4f got %6.4f", want, got)
			}
		}
	}
}

func TestSignint(t *testing.T) {
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			multi := 1
			if i%2 == 0 {
				multi = -1
			}
			valinti := i * multi
			multj := 1
			if j%2 == 0 {
				multj = -1
			}
			valintj := j * multj
			if got, want := signint(valinti, valintj), int(math.Copysign(float64(valinti), float64(valintj))); got != want {
				t.Errorf("signint: values do not match: expected %d got %d", want, got)
			}
		}
	}
}

func TestSinf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		valf32 := 0.25 * float32(i*mult)
		if got, want := sinf32(valf32), float32(math.Sin(float64(valf32))); got != want {
			t.Errorf("sinf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestSqrtf32(t *testing.T) {
	var valf32 float32
	for i := 0; i < 100; i++ {
		valf32 = 0.25 * float32(i)
		if got, want := sqrtf32(valf32), float32(math.Sqrt(float64(valf32))); got != want {
			t.Errorf("sqrt: float32 values do not match: expected %12.10e got %12.10e", want, got)
		}
	}
}
