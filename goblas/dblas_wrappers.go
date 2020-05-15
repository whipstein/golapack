package goblas

// #cgo LDFLAGS: librefblas.a -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9/ -lgfortran
// void drotg_(double*, double*, double*, double*);
// void drotmg_(double*, double*, double*, double*, double*);
// double dnrm2_(int*, double*, int*);
// double dasum_(int*, double*, int*);
// void dscal_(int*, double*, double*, int*);
// int idamax_(int*, double*, int*);
// double ddot_(int*, double*, int*, double*, int*);
// void daxpy_(int*, double*, double*, int*, double*, int*);
// void dcopy_(int*, double*, int*, double*, int*);
// void dswap_(int*, double*, int*, double*, int*);
// void drotm_(int*, double*, int*, double*, int*, double*);
// double dsdot_(int*, float*, int*, float*, int*);
// void drot_(int*, double*, int*, double*, int*, double*, double*);
// void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dgbmv_(char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dsymv_(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dsbmv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dspmv_(char*, int*, double*, double*, double*, int*, double*, double*, int*);
// void dtrmv_(char*, char*, char*, int*, double*, int*, double*, int*);
// void dtbmv_(char*, char*, char*, int*, int*, double*, int*, double*, int*);
// void dtpmv_(char*, char*, char*, int*, double*, double*, int*);
// void dtrsv_(char*, char*, char*, int*, double*, int*, double*, int*);
// void dtbsv_(char*, char*, char*, int*, int*, double*, int*, double*, int*);
// void dtpsv_(char*, char*, char*, int*, double*, double*, int*);
// void dger_(int*, int*, double*, double*, int*, double*, int*, double*, int*);
// void dsyr_(char*, int*, double*, double*, int*, double*, int*);
// void dspr_(char*, int*, double*, double*, int*, double*);
// void dsyr2_(char*, int*, double*, double*, int*, double*, int*, double*, int*);
// void dspr2_(char*, int*, double*, double*, int*, double*, int*, double*);
// void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dsymm_(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
// void dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
// void dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
// void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);
// void dsyr2k_(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
import "C"

func drotgWrapper(a, b, c, s *float64) {
	C.drotg_((*C.double)(a), (*C.double)(b), (*C.double)(c), (*C.double)(s))
}

func drotmgWrapper(d1, d2, x, y *float64, sparam *[]float64) {
	C.drotmg_((*C.double)(d1), (*C.double)(d2), (*C.double)(x), (*C.double)(y), (*C.double)(&(*sparam)[0]))
}

func dnrm2Wrapper(n *int, sx *[]float64, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float64(C.dnrm2_(&_n, (*C.double)(&(*sx)[0]), &_incx))
}

func dnrm2WrapperTest(n *int, sx *[]float64, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float64(C.dnrm2_(&_n, (*C.double)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func dasumWrapper(n *int, sx *[]float64, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float64(C.dasum_(&_n, (*C.double)(&(*sx)[0]), &_incx))
}

func dasumWrapperTest(n *int, sx *[]float64, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float64(C.dasum_(&_n, (*C.double)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func dscalWrapper(n *int, sa *float64, sx *[]float64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dscal_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx)
}

func dscalWrapperTest(n *int, sa *float64, sx *[]float64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dscal_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func idamaxWrapper(n *int, sx *[]float64, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return int(C.idamax_(&_n, (*C.double)(&(*sx)[0]), &_incx))
}

func idamaxWrapperTest(n *int, sx *[]float64, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := int(C.idamax_(&_n, (*C.double)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func ddotWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float64(C.ddot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy))
}

func ddotWrapperTest(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := float64(C.ddot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func daxpyWrapper(n *int, sa *float64, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.daxpy_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func daxpyWrapperTest(n *int, sa *float64, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.daxpy_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dcopyWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dcopy_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func dcopyWrapperTest(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dcopy_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dswapWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dswap_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func dswapWrapperTest(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dswap_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func drotmWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, temp *[]float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drotm_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(&(*temp)[0]))
}

func drotmWrapperTest(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, temp *[]float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drotm_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(&(*temp)[0]))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dsdotWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float64(C.dsdot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
}

func dsdotWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := float64(C.dsdot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func drotWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, sc, ss *float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(sc), (*C.double)(ss))
}

func drotWrapperTest(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, sc, ss *float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(sc), (*C.double)(ss))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dgemvWrapper(trans *byte, m, n *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dgemv_(&_trans, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
}

func dgemvWrapperTest(trans *byte, m, n *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dgemv_(&_trans, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dgbmvWrapper(trans *byte, m, n, kl, ku *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
}

func dgbmvWrapperTest(trans *byte, m, n, kl, ku *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*kl = int(_kl)
	*ku = int(_ku)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dsymvWrapper(uplo *byte, n *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsymv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
}

func dsymvWrapperTest(uplo *byte, n *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsymv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dsbmvWrapper(uplo *byte, n, k *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsbmv_(&_uplo, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
}

func dsbmvWrapperTest(uplo *byte, n, k *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsbmv_(&_uplo, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dspmvWrapper(uplo *byte, n *int, alpha *float64, a *[]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspmv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
}

func dspmvWrapperTest(uplo *byte, n *int, alpha *float64, a *[]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspmv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dtrmvWrapper(uplo, trans, diag *byte, n *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtrmv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
}

func dtrmvWrapperTest(uplo, trans, diag *byte, n *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtrmv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func dtbmvWrapper(uplo, trans, diag *byte, n, k *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
}

func dtbmvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func dtpmvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpmv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
}

func dtpmvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpmv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func dtrsvWrapper(uplo, trans, diag *byte, n *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtrsv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
}

func dtrsvWrapperTest(uplo, trans, diag *byte, n *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtrsv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func dtbsvWrapper(uplo, trans, diag *byte, n, k *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
}

func dtbsvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]float64, lda *int, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dtbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func dtpsvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpsv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
}

func dtpsvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpsv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func dgerWrapper(m, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dger_(&_m, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
}

func dgerWrapperTest(m, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dger_(&_m, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dsyrWrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dsyr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]), &_lda)
}

func dsyrWrapperTest(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dsyr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func dsprWrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dspr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]))
}

func dsprWrapperTest(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dspr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
}

func dsyr2Wrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsyr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
}

func dsyr2WrapperTest(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsyr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dspr2Wrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]))
}

func dspr2WrapperTest(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func dgemmWrapper(transa, transb *byte, m, n, k *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
}

func dgemmWrapperTest(transa, transb *byte, m, n, k *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
	*transa = byte(_transa)
	*transb = byte(_transb)
	*m = int(_m)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func dsymmWrapper(side, uplo *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dsymm_(&_side, &_uplo, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
}

func dsymmWrapperTest(side, uplo *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dsymm_(&_side, &_uplo, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func dtrmmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.dtrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb)
}

func dtrmmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.dtrmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func dtrsmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.dtrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb)
}

func dtrsmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.dtrsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func dsyrkWrapper(uplo, transa *byte, n, k *int, alpha *float64, a *[]float64, lda *int, beta *float64, c *[]float64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.dsyrk_(&_uplo, &_transa, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
}

func dsyrkWrapperTest(uplo, transa *byte, n, k *int, alpha *float64, a *[]float64, lda *int, beta *float64, c *[]float64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.dsyrk_(&_uplo, &_transa, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func dsyr2kWrapper(uplo, trans *byte, n, k *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dsyr2k_(&_uplo, &_trans, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
}

func dsyr2kWrapperTest(uplo, trans *byte, n, k *int, alpha *float64, a *[]float64, lda *int, b *[]float64, ldb *int, beta *float64, c *[]float64, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.dsyr2k_(&_uplo, &_trans, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*b)[0]), &_ldb, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}
