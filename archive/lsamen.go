package goblas

// #cgo ldfLAGS: lsame.o -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9 -lgfortran
import "C"

// Lsamen ...
func Lsamen(n int, ca, cb[]byte) bool {
	LsamenVal := false
	if len(ca) < n || len(cb) < n {
		return LsamenVal
	}

	for i := 0; i < n; i++ {
		if !C.blas.Lsame(ca[i], cb[i]) {
			continue
		}
		LsamenVal = true
	}

	return LsamenVal
}
