package goblas

// #cgo LDFLAGS: librefblas.a -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9/ -lgfortran
// #include <complex.h>
// double dznrm2_(int*, complex double*, int*);
// double dzasum_(int*, complex double*, int*);
// void zscal_(int*, complex double*, complex double*, int*);
// void zdscal_(int*, double*, complex double*, int*);
// int izamax_(int*, complex double*, int*);
// complex double zdotc_(int*, complex double*, int*, complex double*, int*);
// complex double zdotu_(int*, complex double*, int*, complex double*, int*);
// void zaxpy_(int*, complex double*, complex double*, int*, complex double*, int*);
// void zcopy_(int*, complex double*, int*, complex double*, int*);
// void zswap_(int*, complex double*, int*, complex double*, int*);
// void zgemv_(char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zgbmv_(char*, int*, int*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zhemv_(char*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zhbmv_(char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zhpmv_(char*, int*, complex double*, complex double*, complex double*, int*, complex double*, complex double*, int*);
// void ztrmv_(char*, char*, char*, int*, complex double*, int*, complex double*, int*);
// void ztbmv_(char*, char*, char*, int*, int*, complex double*, int*, complex double*, int*);
// void ztpmv_(char*, char*, char*, int*, complex double*, complex double*, int*);
// void ztrsv_(char*, char*, char*, int*, complex double*, int*, complex double*, int*);
// void ztbsv_(char*, char*, char*, int*, int*, complex double*, int*, complex double*, int*);
// void ztpsv_(char*, char*, char*, int*, complex double*, complex double*, int*);
// void zgerc_(int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, int*);
// void zgeru_(int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, int*);
// void zher_(char*, int*, double*, complex double*, int*, complex double*, int*);
// void zhpr_(char*, int*, double*, complex double*, int*, complex double*);
// void zher2_(char*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, int*);
// void zhpr2_(char*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*);
// void zgemm_(char*, char*, int*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zhemm_(char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void zsymm_(char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
// void ztrmm_(char*, char*, char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*);
// void ztrsm_(char*, char*, char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*);
// void zherk_(char*, char*, int*, int*, double*, complex double*, int*, double*, complex double*, int*);
// void zsyrk_(char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, complex double*, int*);
// void zher2k_(char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, double*, complex double*, int*);
// void zsyr2k_(char*, char*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*);
import "C"

func dznrm2Wrapper(n *int, sx *[]complex128, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float64(C.dznrm2_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
}

func dznrm2WrapperTest(n *int, sx *[]complex128, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float64(C.dznrm2_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func dzasumWrapper(n *int, sx *[]complex128, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float64(C.dzasum_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
}

func dzasumWrapperTest(n *int, sx *[]complex128, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float64(C.dzasum_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func zscalWrapper(n *int, sa *complex128, sx *[]complex128, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zscal_(&_n, (*C.complexdouble)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx)
}

func zscalWrapperTest(n *int, sa *complex128, sx *[]complex128, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zscal_(&_n, (*C.complexdouble)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func zdscalWrapper(n *int, sa *float64, sx *[]complex128, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zdscal_(&_n, (*C.double)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx)
}

func zdscalWrapperTest(n *int, sa *float64, sx *[]complex128, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zdscal_(&_n, (*C.double)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func izamaxWrapper(n *int, sx *[]complex128, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return int(C.izamax_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
}

func izamaxWrapperTest(n *int, sx *[]complex128, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := int(C.izamax_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func zdotcWrapper(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) complex128 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return complex128(C.zdotc_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy))
}

func zdotcWrapperTest(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) complex128 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := complex128(C.zdotc_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func zdotuWrapper(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) complex128 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return complex128(C.zdotu_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy))
}

func zdotuWrapperTest(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) complex128 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := complex128(C.zdotu_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func zaxpyWrapper(n *int, sa *complex128, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zaxpy_(&_n, (*C.complexdouble)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
}

func zaxpyWrapperTest(n *int, sa *complex128, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zaxpy_(&_n, (*C.complexdouble)(sa), (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zcopyWrapper(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zcopy_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
}

func zcopyWrapperTest(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zcopy_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zswapWrapper(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zswap_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
}

func zswapWrapperTest(n *int, sx *[]complex128, incx *int, sy *[]complex128, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zswap_(&_n, (*C.complexdouble)(&(*sx)[0]), &_incx, (*C.complexdouble)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zgemvWrapper(trans *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgemv_(&_trans, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
}

func zgemvWrapperTest(trans *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgemv_(&_trans, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zgbmvWrapper(trans *byte, m, n, kl, ku *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
}

func zgbmvWrapperTest(trans *byte, m, n, kl, ku *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*kl = int(_kl)
	*ku = int(_ku)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zhemvWrapper(uplo *byte, n *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhemv_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
}

func zhemvWrapperTest(uplo *byte, n *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhemv_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zhbmvWrapper(uplo *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhbmv_(&_uplo, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
}

func zhbmvWrapperTest(uplo *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhbmv_(&_uplo, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zhpmvWrapper(uplo *byte, n *int, alpha *complex128, a *[]complex128, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhpmv_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
}

func zhpmvWrapperTest(uplo *byte, n *int, alpha *complex128, a *[]complex128, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhpmv_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(beta), (*C.complexdouble)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ztrmvWrapper(uplo, trans, diag *byte, n *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztrmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztrmvWrapperTest(uplo, trans, diag *byte, n *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztrmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ztbmvWrapper(uplo, trans, diag *byte, n, k *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztbmvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ztpmvWrapper(uplo, trans, diag *byte, n *int, a, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ztpmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztpmvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ztpmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func ztrsvWrapper(uplo, trans, diag *byte, n *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztrsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztrsvWrapperTest(uplo, trans, diag *byte, n *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztrsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ztbsvWrapper(uplo, trans, diag *byte, n, k *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztbsvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]complex128, lda *int, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ztbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ztpsvWrapper(uplo, trans, diag *byte, n *int, a, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ztpsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx)
}

func ztpsvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]complex128, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ztpsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexdouble)(&(*a)[0]), (*C.complexdouble)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func zgercWrapper(m, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgerc_(&_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
}

func zgercWrapperTest(m, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgerc_(&_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zgeruWrapper(m, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgeru_(&_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
}

func zgeruWrapperTest(m, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zgeru_(&_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zherWrapper(uplo *byte, n *int, ralpha *float64, x *[]complex128, incx *int, a *[]complex128, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.zher_(&_uplo, &_n, (*C.double)(ralpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*a)[0]), &_lda)
}

func zherWrapperTest(uplo *byte, n *int, ralpha *float64, x *[]complex128, incx *int, a *[]complex128, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.zher_(&_uplo, &_n, (*C.double)(ralpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func zhprWrapper(uplo *byte, n *int, ralpha *float64, x *[]complex128, incx *int, a *[]complex128) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zhpr_(&_uplo, &_n, (*C.double)(ralpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*a)[0]))
}

func zhprWrapperTest(uplo *byte, n *int, ralpha *float64, x *[]complex128, incx *int, a *[]complex128) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.zhpr_(&_uplo, &_n, (*C.double)(ralpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
}

func zher2Wrapper(uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	_lda := C.int(*lda)
	C.zher2_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
}

func zher2WrapperTest(uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	_lda := C.int(*lda)
	C.zher2_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	*lda = int(_lda)
}

func zhpr2Wrapper(uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhpr2_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]))
}

func zhpr2WrapperTest(uplo *byte, n *int, alpha *complex128, x *[]complex128, incx *int, y *[]complex128, incy *int, a *[]complex128) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.zhpr2_(&_uplo, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*x)[0]), &_incx, (*C.complexdouble)(&(*y)[0]), &_incy, (*C.complexdouble)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func zgemmWrapper(transa, transb *byte, m, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zgemmWrapperTest(transa, transb *byte, m, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*transa = byte(_transa)
	*transb = byte(_transb)
	*m = int(_m)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func zhemmWrapper(side, uplo *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zhemm_(&_side, &_uplo, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zhemmWrapperTest(side, uplo *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zhemm_(&_side, &_uplo, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func zsymmWrapper(side, uplo *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zsymm_(&_side, &_uplo, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zsymmWrapperTest(side, uplo *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zsymm_(&_side, &_uplo, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func ztrmmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ztrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb)
}

func ztrmmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ztrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func ztrsmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ztrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb)
}

func ztrsmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ztrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func zherkWrapper(uplo, trans *byte, n, k *int, ralpha *float64, a *[]complex128, lda *int, rbeta *float64, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.zherk_(&_uplo, &_trans, &_n, &_k, (*C.double)(ralpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.double)(rbeta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zherkWrapperTest(uplo, trans *byte, n, k *int, ralpha *float64, a *[]complex128, lda *int, rbeta *float64, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.zherk_(&_uplo, &_trans, &_n, &_k, (*C.double)(ralpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.double)(rbeta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func zsyrkWrapper(uplo, transa *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, beta *complex128, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.zsyrk_(&_uplo, &_transa, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zsyrkWrapperTest(uplo, transa *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, beta *complex128, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.zsyrk_(&_uplo, &_transa, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func zher2kWrapper(uplo, transa *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, rbeta *float64, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zher2k_(&_uplo, &_transa, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.double)(rbeta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zher2kWrapperTest(uplo, transa *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, rbeta *float64, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zher2k_(&_uplo, &_transa, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.double)(rbeta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func zsyr2kWrapper(uplo, trans *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zsyr2k_(&_uplo, &_trans, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
}

func zsyr2kWrapperTest(uplo, trans *byte, n, k *int, alpha *complex128, a *[]complex128, lda *int, b *[]complex128, ldb *int, beta *complex128, c *[]complex128, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.zsyr2k_(&_uplo, &_trans, &_n, &_k, (*C.complexdouble)(alpha), (*C.complexdouble)(&(*a)[0]), &_lda, (*C.complexdouble)(&(*b)[0]), &_ldb, (*C.complexdouble)(beta), (*C.complexdouble)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}
