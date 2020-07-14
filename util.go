package golapack

type memory struct {
	combla struct {
		icase int
		n     int
		incx  int
		incy  int
		off   int
		pass  bool
	}
	begc struct {
		mi int
		ic int
		i  int
	}
}

var common memory

func maxlocf32(a []float32) (idx int) {
	if len(a) == 1 {
		return
	}
	amax := a[0]
	for i := 1; i < len(a); i++ {
		if a[i] > amax {
			idx = i
			amax = a[i]
		}
	}
	return
}
