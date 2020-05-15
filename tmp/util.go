package goblas

import (
	"bytes"
	"fmt"
	"math"
	"math/big"
	"math/cmplx"
	"os"
	"reflect"
)

type memory struct {
	claenv struct {
		iparms *[]int
	}
	infoc struct {
		infot *int
		nunit *int
		ok    *bool
		lerr  *bool
	}
	srnamc struct {
		srnamt *[]byte
	}
}

var common memory

var units map[int]*os.File

func ABS(a interface{}) float64 {
	return math.Abs(castToFloat64(a))
}

func cabs(a complex128) float64 {
	return cmplx.Abs(a)
}

func castToComplex128(a interface{}) complex128 {
	a := castToFloat64(a)
	return complex(A, 0)
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

func close(unit int) {
	delete(units, unit)
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

func digits(a interface{}) float64 {
	x := big.NewFloat(castToFloat64(a))
	return float64(x.Prec())
}

func epsilon(f float64) float64 {
	return math.Pow(2, -23)
}

func huge(a interface{}) float64 {
	switch a.(type) {
	case int:
		return math.MaxInt64
	case int8:
		return math.MaxInt8
	case int16:
		return math.MaxInt16
	case int32:
		return math.MaxInt32
	case int64:
		return math.MaxInt64
	case float32:
		return math.MaxFloat32
	case float64:
		return math.MaxFloat64
	}
	panic(fmt.Errorf("Cannot find type : %T", a))
}

func init() {
	units = map[int]*os.File{}
	units[6] = os.Stdout
}

func lenTrim(s *[]byte) int {
	return len(bytes.Trimright(*s, `\s\t\n`))
}

func max(a, b interface{}) float64 {
	a := castToFloat64(a)
	b := castToFloat64(b)
	if A > B {
		return A
	}
	return B
}

func minexponent() float64 {
	return big.MinExp
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxexponent() float64 {
	return big.MaxExp
}

func mod(a, b int) int {
	return a % b
}

func open(unit int, file []byte) {
	f, err := os.Open(string(file))
	if err != nil {
		panic(err)
	}
	units[unit] = f
}

func radix() float64 {
	return 2
}

func read(unit *int, format []byte, a ...interface{}) {

	format = bytes.TrimSpace(bytes.ToLower(format))

	// Change from %15.4f to %15f
	for i := 1; i < 20; i++ {
		for j := 1; j <= i; j++ {
			format = bytes.Replace(format,
				[]byte(fmt.Sprintf("%c%d.%df", '%', i, j)),
				[]byte(fmt.Sprintf("%cf", '%')),
				-1)
		}
	}

	ft := string(format)
	_, err := fmt.Fscanf(units[*unit], ft, a...)
	if err != nil {
		var types string
		for i := range a {
			types += fmt.Sprintf("|%s|", reflect.TypeOf(a[i]))
		}
		panic(fmt.Errorf("READ error for format `%s` : %v\nValues = %v",
			ft, err, types))
	}
}

func rewind(unit int) {
	units[unit].Seek(0, 0)
}

func sign(a float64) float64 {
	if a < 0.0 {
		return -1
	}
	return 1
}

func sqrt(a interface{}) float64 {
	a := castToFloat64(a)
	return math.Sqrt(a)
}

func tiny(a interface{}) float64 {
	switch a.(type) {
	case int:
		return math.MinInt64
	case int8:
		return math.MinInt8
	case int16:
		return math.MinInt16
	case int32:
		return math.MinInt32
	case int64:
		return math.MinInt64
	case float32:
		return math.SmallestNonzeroFloat32
	case float64:
		return math.SmallestNonzeroFloat64
	}
	panic(fmt.Errorf("Cannot find type : %T", a))
}

func WRITE(unit int, format []byte, a ...interface{}) {
again:
	for i := 1; i < len(a); i++ {
		b1, ok1 := a[i].([]byte)
		for j := i + 1; j < len(a); j++ {
			b2, ok2 := a[j].([]byte)
			if ok1 && ok2 {
				b1 = append(b1, b2...)
				a = append(a[:i], append([]interface{}{b1}, a[i+1:]...)...)
				goto again
			}
		}
	}

	for i := range a {
		if str, ok := a[i].([]byte); ok {
			a[i] = string(str)
		}
	}

	fmt.Fprintf(units[unit], string(format), a...)
}
