package golapack

import "C"
import (
	"log"
	"math"
	"math/cmplx"
	"reflect"
)

const (
	uvnan    = 0xFF800001
	uvinf    = 0x7F800000
	uvneginf = 0xFF800000
	uvone    = 0x3F800000
	mask     = 0xFF
	shift    = 32 - 8 - 1
	bias     = 127
	signMask = 1 << 31
	fracMask = 1<<shift - 1
)

func absc64(a complex64) float32 {
	return float32(cmplx.Abs(complex128(a)))
}

func absc128(a complex128) float64 {
	return cmplx.Abs(a)
}

func absf32(a float32) float32 {
	if a < 0 {
		return -a
	}
	return a
	// return float32(math.Abs(float64(a)))
}

func absf64(a float64) float64 {
	if a < 0 {
		return -a
	}
	return a
	// return float32(math.Abs(float64(a)))
}

func absint(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func abssum(a interface{}) interface{} {
	aval := reflect.ValueOf(a)

	switch a.(type) {
	case *complex64:
		x := reflect.Indirect(aval)
		return float32(math.Abs(real(x.Complex())) + math.Abs(imag(x.Complex())))
	case complex64:
		return float32(math.Abs(real(aval.Complex())) + math.Abs(imag(aval.Complex())))
	case *complex128:
		x := reflect.Indirect(aval)
		return math.Abs(real(x.Complex())) + math.Abs(imag(x.Complex()))
	case complex128:
		return math.Abs(real(aval.Complex())) + math.Abs(imag(aval.Complex()))
	case *C.complexfloat:
		x := reflect.Indirect(aval)
		return C.float(math.Abs(real(x.Complex())) + math.Abs(imag(x.Complex())))
	case C.complexfloat:
		return C.float(math.Abs(real(aval.Complex())) + math.Abs(imag(aval.Complex())))
	case *C.complexdouble:
		x := reflect.Indirect(aval)
		return C.double(math.Abs(real(x.Complex())) + math.Abs(imag(x.Complex())))
	case C.complexdouble:
		return C.double(math.Abs(real(aval.Complex())) + math.Abs(imag(aval.Complex())))
	default:
		log.Panic("abssum: unrecognized parameter type: ", reflect.TypeOf(a))
	}
	return nil
}

func abssumc64(a complex64) float32 {
	return float32(math.Abs(float64(real(a))) + math.Abs(float64(imag(a))))
}

func abssumc128(a complex128) float64 {
	return math.Abs(real(a)) + math.Abs(imag(a))
}

func ceilingf32(a float32) int {
	return int(math.Ceil(float64(a)))
}

func ceilingf64(a float64) int {
	return int(math.Ceil(a))
}

func conjc64(a complex64) complex64 {
	return complex64(cmplx.Conj(complex128(a)))
}

func conjc128(a complex128) complex128 {
	return cmplx.Conj(a)
}

func cosf32(a float32) float32 {
	return float32(math.Cos(float64(a)))
}

func cosf64(a float64) float64 {
	return math.Cos(a)
}

func epsilonf32() float32 {
	return powf32(2, -23)
	// return float32(1.1920929e-07)
}

func epsilonf64() float64 {
	return powf64(2, -52)
}

func expf32(a float32) float32 {
	return float32(math.Exp(float64(a)))
}

func expf64(a float64) float64 {
	return math.Exp(a)
}

func logf32(a float32) float32 {
	return float32(math.Log(float64(a)))
}

func logf64(a float64) float64 {
	return math.Log(a)
}

func log10f32(a float32) float32 {
	return float32(math.Log10(float64(a)))
}

func log10f64(a float64) float64 {
	return math.Log10(a)
}

func maxf32(a, b float32) float32 {
	if a >= b {
		return a
	}
	return b
}

func maxf64(a, b float64) float64 {
	if a >= b {
		return a
	}
	return b
}

func maxint(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

func minf32(a, b float32) float32 {
	if a <= b {
		return a
	}
	return b
}

func minf64(a, b float64) float64 {
	if a <= b {
		return a
	}
	return b
}

func minint(a, b int) int {
	if a <= b {
		return a
	}
	return b
}

func modint(a, b int) int {
	return a % b
}

func nintf32(a float32) int {
	return int(math.Round(float64(a)))
}

func nintf64(a float64) int {
	return int(math.Round(a))
}

func powf32(a, b float32) float32 {
	return float32(math.Pow(float64(a), float64(b)))
}

func powf64(a, b float64) float64 {
	return math.Pow(a, b)
}

func powint(a, b int) int {
	return int(math.Pow(float64(a), float64(b)))
}

func signf32(a, b float32) float32 {
	return float32(math.Copysign(float64(a), float64(b)))
}

func signf64(a, b float64) float64 {
	return math.Copysign(a, b)
}

func signint(a, b int) int {
	return int(math.Copysign(float64(a), float64(b)))
}

func sinf32(a float32) float32 {
	return float32(math.Sin(float64(a)))
}

func sinf64(a float64) float64 {
	return math.Sin(a)
}

func sqrtf32(a float32) float32 {
	return float32(math.Sqrt(float64(a)))
}

func sqrtf64(a float64) float64 {
	return math.Sqrt(a)
}
