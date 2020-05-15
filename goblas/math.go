package goblas

import "C"
import (
	"math"
	"math/cmplx"
)

func absc64(a complex64) float32 {
	return float32(cmplx.Abs(complex128(a)))
}

func absc128(a complex128) float64 {
	return cmplx.Abs(a)
}

func absccdouble(a C.complexdouble) C.double {
	return C.double(cmplx.Abs(complex128(a)))
}

func absccfloat(a C.complexfloat) C.float {
	return C.float(cmplx.Abs(complex128(a)))
}

func abscdouble(a C.double) C.double {
	return C.double(math.Abs(float64(a)))
}

func abscfloat(a C.float) C.float {
	return C.float(math.Abs(float64(a)))
}

func abscint(a C.int) C.int {
	if a < 0 {
		return -a
	}
	return a
}

func absf32(a float32) float32 {
	return float32(math.Abs(float64(a)))
}

func absf64(a float64) float64 {
	return math.Abs(a)
}

func absint(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func abssumcdouble(a C.complexdouble) C.double {
	return (abscdouble(realcdouble(a)) + abscdouble(imagcdouble(a)))
}

func abssumcfloat(a C.complexfloat) C.float {
	return (abscfloat(realcfloat(a)) + abscfloat(imagcfloat(a)))
}

func abssumf32(a complex64) float32 {
	return float32(absf32(realc64(a)) + absf32(imagc64(a)))
}

func abssumf64(a complex128) float64 {
	return (absf64(real(a)) + absf64(imag(a)))
}

func conjgcdouble(c C.complexdouble) C.complexdouble {
	return C.complexdouble(cmplx.Conj(complex128(c)))
}

func conjgcfloat(c C.complexfloat) C.complexfloat {
	return C.complexfloat(cmplx.Conj(complex128(c)))
}

func conjgc64(c complex64) complex64 {
	return complex64(cmplx.Conj(complex128(c)))
}

func conjgc128(c complex128) complex128 {
	return cmplx.Conj(c)
}

func cabs(a complex128) float64 {
	return cmplx.Abs(a)
}

func cabscdouble(a C.complexdouble) C.double {
	return C.double(cmplx.Abs(complex128(a)))
}

func cabscfloat(a C.complexfloat) C.float {
	return C.float(cmplx.Abs(complex128(a)))
}

func cabsf32(a complex64) float32 {
	return float32(cmplx.Abs(complex128(a)))
}

func cmplxcdouble(a C.double, b C.double) C.complexdouble {
	return C.complexdouble(complex(float64(a), float64(b)))
}

func cmplxcfloat(a C.float, b C.float) C.complexfloat {
	return C.complexfloat(complex(float64(a), float64(b)))
}

func cmplxc64(a float32, b float32) complex64 {
	return complex64(complex(float64(a), float64(b)))
}

func cmplxc128(a float64, b float64) complex128 {
	return complex(a, b)
}

func epsiloncdouble() C.double {
	return C.double(epsilonf64())
}

func epsiloncfloat() C.float {
	return C.float(epsilonf32())
}

func epsilonf32() float32 {
	return 1.1920929e-07
}

func epsilonf64() float64 {
	return 0.22204460E-15
}

func imagcdouble(a C.complexdouble) C.double {
	return C.double(imag(complex128(a)))
}

func imagcfloat(a C.complexfloat) C.float {
	return C.float(imag(complex128(a)))
}

func imagc64(a complex64) float32 {
	return float32(imag(complex128(a)))
}

func imagc128(a complex128) float64 {
	return imag(a)
}

func int2cdouble(a C.int) C.complexdouble {
	return C.complexdouble(complex(float64(int(a)), 0))
}

func int2cfloat(a C.int) C.complexfloat {
	return C.complexfloat(complex(float64(int(a)), 0))
}

func int2c64(a int) complex64 {
	return complex64(complex(float64(a), 0))
}

func int2c128(a int) complex128 {
	return complex(float64(a), 0)
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func maxcdouble(a, b C.double) C.double {
	if a > b {
		return a
	}
	return b
}

func maxcfloat(a, b C.float) C.float {
	if a > b {
		return a
	}
	return b
}

func maxcint(a, b C.int) C.int {
	if a > b {
		return a
	}
	return b
}

func maxf32(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}

func maxf64(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func mincdouble(a, b C.double) C.double {
	if a < b {
		return a
	}
	return b
}

func mincfloat(a, b C.float) C.float {
	if a < b {
		return a
	}
	return b
}

func minf32(a, b float32) float32 {
	if a < b {
		return a
	}
	return b
}

func minf64(a, b float64) float64 {
	if a < b {
		return a
	}
	return b
}

func mincint(a, b C.int) C.int {
	if a < b {
		return a
	}
	return b
}

func powcdouble(a, b C.double) C.double {
	return C.double(math.Pow(float64(a), float64(b)))
}

func powcfloat(a, b C.float) C.float {
	return C.float(math.Pow(float64(a), float64(b)))
}

func powf32(a, b float32) float32 {
	return float32(math.Pow(float64(a), float64(b)))
}

func powf64(a, b float64) float64 {
	return math.Pow(a, b)
}

func realcdouble(a C.complexdouble) C.double {
	return C.double(real(complex128(a)))
}

func realcfloat(a C.complexfloat) C.float {
	return C.float(real(complex128(a)))
}

func realc64(a complex64) float32 {
	return float32(real(complex128(a)))
}

func realc128(a complex128) float64 {
	return real(a)
}

func signf32(a float32, b float32) float32 {
	if (a < 0.0 && b < 0.0) || (a > 0.0 && b > 0.0) || (a == 0.0 && b == 0.0) {
		return a
	}
	return -a
}

func signf64(a float64, b float64) float64 {
	if (a < 0.0 && b < 0.0) || (a > 0.0 && b > 0.0) || (a == 0.0 && b == 0.0) {
		return a
	}
	return -a
}

func signint(a int, b int) int {
	if (a < 0 && b < 0) || (a > 0 && b > 0) || (a == 0 && b == 0) {
		return a
	}
	return -a
}

func sqrtcdouble(a C.double) C.double {
	return C.double(math.Sqrt(float64(a)))
}

func sqrtcfloat(a C.float) C.float {
	return C.float(math.Sqrt(float64(a)))
}

func sqrtf32(a float32) float32 {
	return float32(math.Sqrt(float64(a)))
}

func sqrtf64(a float64) float64 {
	return math.Sqrt(a)
}

func sqrtint(a int) int {
	return int(math.Sqrt(float64(a)))
}
