package golapack

type memory struct {
	begc struct {
		mi int
		mj int
		ic int
		i  int
		j  int
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
