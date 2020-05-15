// Copyright (c) 2020 Matt Braunstein

package archive

/*
#cgo cfLAGS: -g -O2
#cgo ldfLAGS: -L/opt/OpenBLAS/lib -lopenblas
#include </opt/OpenBLAS/include/lapacke.h>
*/
// import "C"
import (
	"fmt"
	// "testing"
)

// func TestCaxpy(t *testing.T) {
// 	cx :=[]complex64{1, 1, 1, 2, 3, 4, 3, 5, 2, 4}
// 	cy :=[]complex64{-10, -3, 12, 14, 14, 12, 16, 16, 18, 16}
// 	exemplar :=[]complex64{-8 + 2i, -1 + 2i, 14 + 2i, 18 + 4i, 20 + 6i, 20 + 8i, 22 + 6i, 26 + 10i, 22 + 4i, 24 + 8i}
// 	var ca complex64
// 	var n, incX, incY int

// 	n = 3
// 	ca = 2 + 2i
// 	incX = 1
// 	incY = 1

// 	Caxpy(n, ca, cx, incX, cy, incY)
// 	// resultwork := Caxpywork(n, ca, cx, incX, cy, incY)

// 	for i, v := range cy {
// 		if v != exemplar[i] {
// 			t.Errorf("incorrect results: got %f, want %f", v, exemplar[i])
// 		}
// 	}
// }

func main() {
	cx :=[]complex64{1, 1, 1, 2, 3, 4, 3, 5, 2, 4}
	cy :=[]complex64{-10, -3, 12, 14, 14, 12, 16, 16, 18, 16}
	exemplar :=[]complex64{-8 + 2i, -1 + 2i, 14 + 2i, 18 + 4i, 20 + 6i, 20 + 8i, 22 + 6i, 26 + 10i, 22 + 4i, 24 + 8i}
	var ca complex64
	var n, incX, incY int

	n = 3
	ca = 2 + 2i
	incX = 1
	incY = 1

	Caxpy(n, ca, cx, incX, cy, incY)
	// resultwork := Caxpywork(n, ca, cx, incX, cy, incY)

	for i, v := range cy {
		fmt.Printf("incorrect results: got %f, want %f", v, exemplar[i])
		// if v != exemplar[i] {
		// 	t.Errorf("incorrect results: got %f, want %f", v, exemplar[i])
		// }
	}
}
