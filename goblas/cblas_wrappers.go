package goblas

// #cgo LDFLAGS: librefblas.a -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9/ -lgfortran
// #include <complex.h>
// float scnrm2_(int*, complex float*, int*);
// float scasum_(int*, complex float*, int*);
// void cscal_(int*, complex float*, complex float*, int*);
// void csscal_(int*, float*, complex float*, int*);
// int icamax_(int*, complex float*, int*);
// complex float cdotc_(int*, complex float*, int*, complex float*, int*);
// complex float cdotu_(int*, complex float*, int*, complex float*, int*);
// void caxpy_(int*, complex float*, complex float*, int*, complex float*, int*);
// void ccopy_(int*, complex float*, int*, complex float*, int*);
// void cswap_(int*, complex float*, int*, complex float*, int*);
// void cgemv_(char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void cgbmv_(char*, int*, int*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void chemv_(char*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void chbmv_(char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void chpmv_(char*, int*, complex float*, complex float*, complex float*, int*, complex float*, complex float*, int*);
// void ctrmv_(char*, char*, char*, int*, complex float*, int*, complex float*, int*);
// void ctbmv_(char*, char*, char*, int*, int*, complex float*, int*, complex float*, int*);
// void ctpmv_(char*, char*, char*, int*, complex float*, complex float*, int*);
// void ctrsv_(char*, char*, char*, int*, complex float*, int*, complex float*, int*);
// void ctbsv_(char*, char*, char*, int*, int*, complex float*, int*, complex float*, int*);
// void ctpsv_(char*, char*, char*, int*, complex float*, complex float*, int*);
// void cgerc_(int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, int*);
// void cgeru_(int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, int*);
// void cher_(char*, int*, float*, complex float*, int*, complex float*, int*);
// void chpr_(char*, int*, float*, complex float*, int*, complex float*);
// void cher2_(char*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, int*);
// void chpr2_(char*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*);
// void cgemm_(char*, char*, int*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void chemm_(char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void csymm_(char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
// void ctrmm_(char*, char*, char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*);
// void ctrsm_(char*, char*, char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*);
// void cherk_(char*, char*, int*, int*, float*, complex float*, int*, float*, complex float*, int*);
// void csyrk_(char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, complex float*, int*);
// void cher2k_(char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, float*, complex float*, int*);
// void csyr2k_(char*, char*, int*, int*, complex float*, complex float*, int*, complex float*, int*, complex float*, complex float*, int*);
import "C"

func scnrm2Wrapper(n *int, sx *[]complex64, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float32(C.scnrm2_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
}

func scnrm2WrapperTest(n *int, sx *[]complex64, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float32(C.scnrm2_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func scasumWrapper(n *int, sx *[]complex64, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float32(C.scasum_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
}

func scasumWrapperTest(n *int, sx *[]complex64, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float32(C.scasum_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func cscalWrapper(n *int, sa *complex64, sx *[]complex64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.cscal_(&_n, (*C.complexfloat)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx)
}

func cscalWrapperTest(n *int, sa *complex64, sx *[]complex64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.cscal_(&_n, (*C.complexfloat)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func csscalWrapper(n *int, sa *float32, sx *[]complex64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.csscal_(&_n, (*C.float)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx)
}

func csscalWrapperTest(n *int, sa *float32, sx *[]complex64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.csscal_(&_n, (*C.float)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func icamaxWrapper(n *int, sx *[]complex64, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return int(C.icamax_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
}

func icamaxWrapperTest(n *int, sx *[]complex64, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := int(C.icamax_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func cdotcWrapper(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) complex64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return complex64(C.cdotc_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy))
}

func cdotcWrapperTest(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) complex64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := complex64(C.cdotc_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func cdotuWrapper(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) complex64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return complex64(C.cdotu_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy))
}

func cdotuWrapperTest(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) complex64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := complex64(C.cdotu_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func caxpyWrapper(n *int, sa *complex64, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.caxpy_(&_n, (*C.complexfloat)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
}

func caxpyWrapperTest(n *int, sa *complex64, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.caxpy_(&_n, (*C.complexfloat)(sa), (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ccopyWrapper(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ccopy_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
}

func ccopyWrapperTest(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ccopy_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cswapWrapper(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cswap_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
}

func cswapWrapperTest(n *int, sx *[]complex64, incx *int, sy *[]complex64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cswap_(&_n, (*C.complexfloat)(&(*sx)[0]), &_incx, (*C.complexfloat)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cgemvWrapper(trans *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgemv_(&_trans, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
}

func cgemvWrapperTest(trans *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgemv_(&_trans, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cgbmvWrapper(trans *byte, m, n, kl, ku *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
}

func cgbmvWrapperTest(trans *byte, m, n, kl, ku *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*kl = int(_kl)
	*ku = int(_ku)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func chemvWrapper(uplo *byte, n *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chemv_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
}

func chemvWrapperTest(uplo *byte, n *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chemv_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func chbmvWrapper(uplo *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chbmv_(&_uplo, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
}

func chbmvWrapperTest(uplo *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chbmv_(&_uplo, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func chpmvWrapper(uplo *byte, n *int, alpha *complex64, a *[]complex64, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chpmv_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
}

func chpmvWrapperTest(uplo *byte, n *int, alpha *complex64, a *[]complex64, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chpmv_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(beta), (*C.complexfloat)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ctrmvWrapper(uplo, trans, diag *byte, n *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctrmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctrmvWrapperTest(uplo, trans, diag *byte, n *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctrmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ctbmvWrapper(uplo, trans, diag *byte, n, k *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctbmvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ctpmvWrapper(uplo, trans, diag *byte, n *int, a, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ctpmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctpmvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ctpmv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func ctrsvWrapper(uplo, trans, diag *byte, n *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctrsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctrsvWrapperTest(uplo, trans, diag *byte, n *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctrsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ctbsvWrapper(uplo, trans, diag *byte, n, k *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctbsvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]complex64, lda *int, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ctbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ctpsvWrapper(uplo, trans, diag *byte, n *int, a, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ctpsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx)
}

func ctpsvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]complex64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.ctpsv_(&_uplo, &_trans, &_diag, &_n, (*C.complexfloat)(&(*a)[0]), (*C.complexfloat)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func cgercWrapper(m, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgerc_(&_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
}

func cgercWrapperTest(m, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgerc_(&_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cgeruWrapper(m, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgeru_(&_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
}

func cgeruWrapperTest(m, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.cgeru_(&_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cherWrapper(uplo *byte, n *int, ralpha *float32, x *[]complex64, incx *int, a *[]complex64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.cher_(&_uplo, &_n, (*C.float)(ralpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*a)[0]), &_lda)
}

func cherWrapperTest(uplo *byte, n *int, ralpha *float32, x *[]complex64, incx *int, a *[]complex64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.cher_(&_uplo, &_n, (*C.float)(ralpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func chprWrapper(uplo *byte, n *int, ralpha *float32, x *[]complex64, incx *int, a *[]complex64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.chpr_(&_uplo, &_n, (*C.float)(ralpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*a)[0]))
}

func chprWrapperTest(uplo *byte, n *int, ralpha *float32, x *[]complex64, incx *int, a *[]complex64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.chpr_(&_uplo, &_n, (*C.float)(ralpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
}

func cher2Wrapper(uplo *byte, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	_lda := C.int(*lda)
	C.cher2_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
}

func cher2WrapperTest(uplo *byte, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	_lda := C.int(*lda)
	C.cher2_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	*lda = int(_lda)
}

func chpr2Wrapper(uplo *byte, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chpr2_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]))
}

func chpr2WrapperTest(uplo *byte, n *int, alpha *complex64, x *[]complex64, incx *int, y *[]complex64, incy *int, a *[]complex64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.chpr2_(&_uplo, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*x)[0]), &_incx, (*C.complexfloat)(&(*y)[0]), &_incy, (*C.complexfloat)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func cgemmWrapper(transa, transb *byte, m, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.cgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func cgemmWrapperTest(transa, transb *byte, m, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.cgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*transa = byte(_transa)
	*transb = byte(_transb)
	*m = int(_m)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func chemmWrapper(side, uplo *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.chemm_(&_side, &_uplo, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func chemmWrapperTest(side, uplo *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.chemm_(&_side, &_uplo, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func csymmWrapper(side, uplo *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.csymm_(&_side, &_uplo, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func csymmWrapperTest(side, uplo *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.csymm_(&_side, &_uplo, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func ctrmmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ctrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb)
}

func ctrmmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ctrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func ctrsmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ctrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb)
}

func ctrsmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.ctrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func cherkWrapper(uplo, trans *byte, n, k *int, ralpha *float32, a *[]complex64, lda *int, rbeta *float32, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.cherk_(&_uplo, &_trans, &_n, &_k, (*C.float)(ralpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.float)(rbeta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func cherkWrapperTest(uplo, trans *byte, n, k *int, ralpha *float32, a *[]complex64, lda *int, rbeta *float32, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.cherk_(&_uplo, &_trans, &_n, &_k, (*C.float)(ralpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.float)(rbeta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func csyrkWrapper(uplo, transa *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, beta *complex64, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.csyrk_(&_uplo, &_transa, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func csyrkWrapperTest(uplo, transa *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, beta *complex64, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.csyrk_(&_uplo, &_transa, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func cher2kWrapper(uplo, transa *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, rbeta *float32, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.cher2k_(&_uplo, &_transa, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.float)(rbeta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func cher2kWrapperTest(uplo, transa *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, rbeta *float32, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.cher2k_(&_uplo, &_transa, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.float)(rbeta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func csyr2kWrapper(uplo, trans *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.csyr2k_(&_uplo, &_trans, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
}

func csyr2kWrapperTest(uplo, trans *byte, n, k *int, alpha *complex64, a *[]complex64, lda *int, b *[]complex64, ldb *int, beta *complex64, c *[]complex64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.csyr2k_(&_uplo, &_trans, &_n, &_k, (*C.complexfloat)(alpha), (*C.complexfloat)(&(*a)[0]), &_lda, (*C.complexfloat)(&(*b)[0]), &_ldb, (*C.complexfloat)(beta), (*C.complexfloat)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}
