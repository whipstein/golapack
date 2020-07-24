package golapack

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestAbsc64(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		if got, want := absc64(val), float32(cmplx.Abs(complex128(val))); got != want {
			t.Errorf("absc64: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestAbsc128(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := absc128(val), cmplx.Abs(val); got != want {
			t.Errorf("absc128: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestAbsf32(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := 0.25 * float32(i*mult)
		if got, want := absf32(val), float32(math.Abs(float64(val))); got != want {
			t.Errorf("absf32: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestAbsf64(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := 0.25 * float64(i*mult)
		if got, want := absf64(val), math.Abs(val); got != want {
			t.Errorf("absf64: values do not match: expected %6.4f got %6.4f", want, got)
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

func TestAbssumc64(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		if got, want := abssumc64(val), float32(math.Abs(float64(real(val)))+math.Abs(float64(imag(val)))); got != want {
			t.Errorf("abssumc64: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestAbssumc128(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := abssumc128(val), math.Abs(real(val))+math.Abs(imag(val)); got != want {
			t.Errorf("abssumc128: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestConjc64(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float32(i*mult), 0.25*float32(i*mult))
		if got, want := conjc64(val), complex64(cmplx.Conj(complex128(val))); got != want {
			t.Errorf("conjc64: values do not match: expected %6.4f got %6.4f", want, got)
		}
	}
}

func TestConjc128(t *testing.T) {
	for i := 0; i < 100; i++ {
		mult := 1
		if i%2 == 0 {
			mult = -1
		}
		val := complex(-0.25*float64(i*mult), 0.25*float64(i*mult))
		if got, want := conjc128(val), cmplx.Conj(val); got != want {
			t.Errorf("conjc128: values do not match: expected %6.4f got %6.4f", want, got)
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
