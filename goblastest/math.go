package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func castToFloat64(w interface{}) float64 {
	switch v := w.(type) {
	case float64:
		return v
	case *float64:
		return *v
	case float32:
		return float64(v)
	case *float32:
		return float64(*v)
	case int:
		return float64(v)
	case *int:
		return float64(*v)
		// 	case complex128:
		// 		r := real(v)
		// 		i := imag(v)
		// 		return float64(math.Mod(r, 2) - math.Mod(i, 2))
		// 	case *complex128:
		// 		r := real(*v)
		// 		i := imag(*v)
		// 		return float64(math.Mod(r, 2) - math.Mod(i, 2))
	default:
		panic(fmt.Errorf("cannot cast: %#v", w))
	}
}

func sqrt(a interface{}) float64 {
	A := castToFloat64(a)
	return math.Sqrt(A)
}

func max(a, b int) int {
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

func epsilonf32() float32 {
	return 0.11920929E-06
}

func epsilonf64() float64 {
	return 0.22204460E-15
}

func conjg(c complex128) complex128 {
	return complex128(cmplx.Conj(complex128(c)))
}

func dconjg(c complex128) complex128 {
	return cmplx.Conj(c)
}

func dble(a interface{}) float64 {
	switch a.(type) {
	case int:
		return float64(a.(int))
	case int32:
		return float64(a.(int32))
	case int64:
		return float64(a.(int64))
	case float32:
		return float64(a.(float32))
	case complex64:
		return float64(real(a.(complex64)))
	case complex128:
		return float64(real(a.(complex128)))
	case float64:
		return a.(float64)
	}
	panic(fmt.Errorf("Cannot find type : %T", a))
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

func cabs(a complex128) float64 {
	return cmplx.Abs(a)
}

func powf32(a, b float32) float32 {
	return float32(math.Pow(float64(a), float64(b)))
}

func powf64(a, b float64) float64 {
	return math.Pow(a, b)
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

func sqrtf32(a float32) float32 {
	return float32(math.Sqrt(float64(a)))
}

func sqrtf64(a float64) float64 {
	return math.Sqrt(a)
}

func mod(a, b int) int {
	return a % b
}

func cmplxs(a interface{}) complex128 {
	A := castToFloat64(a)
	return complex(A, 0)
}
