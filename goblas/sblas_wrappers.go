package goblas

// #cgo LDFLAGS: librefblas.a -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9/ -lgfortran
// void srotg_(float*, float*, float*, float*);
// void srotmg_(float*, float*, float*, float*, float*);
// float snrm2_(int*, float*, int*);
// float sasum_(int*, float*, int*);
// void sscal_(int*, float*, float*, int*);
// int isamax_(int*, float*, int*);
// float sdot_(int*, float*, int*, float*, int*);
// void saxpy_(int*, float*, float*, int*, float*, int*);
// void scopy_(int*, float*, int*, float*, int*);
// void sswap_(int*, float*, int*, float*, int*);
// void srotm_(int*, float*, int*, float*, int*, float*);
// float sdsdot_(int*, float*, float*, int*, float*, int*);
// void srot_(int*, float*, int*, float*, int*, float*, float*);
// void sgemv_(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void sgbmv_(char*, int*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void ssymv_(char*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void ssbmv_(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void sspmv_(char*, int*, float*, float*, float*, int*, float*, float*, int*);
// void strmv_(char*, char*, char*, int*, float*, int*, float*, int*);
// void stbmv_(char*, char*, char*, int*, int*, float*, int*, float*, int*);
// void stpmv_(char*, char*, char*, int*, float*, float*, int*);
// void strsv_(char*, char*, char*, int*, float*, int*, float*, int*);
// void stbsv_(char*, char*, char*, int*, int*, float*, int*, float*, int*);
// void stpsv_(char*, char*, char*, int*, float*, float*, int*);
// void sger_(int*, int*, float*, float*, int*, float*, int*, float*, int*);
// void ssyr_(char*, int*, float*, float*, int*, float*, int*);
// void sspr_(char*, int*, float*, float*, int*, float*);
// void ssyr2_(char*, int*, float*, float*, int*, float*, int*, float*, int*);
// void sspr2_(char*, int*, float*, float*, int*, float*, int*, float*);
// void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void ssymm_(char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
// void strmm_(char*, char*, char*, char*, int*, int*, float*, float*, int*, float*, int*);
// void strsm_(char*, char*, char*, char*, int*, int*, float*, float*, int*, float*, int*);
// void ssyrk_(char*, char*, int*, int*, float*, float*, int*, float*, float*, int*);
// void ssyr2k_(char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
import "C"

func srotgWrapper(a, b, c, s *float32) {
	C.srotg_((*C.float)(a), (*C.float)(b), (*C.float)(c), (*C.float)(s))
}

func srotmgWrapper(d1, d2, x, y *float32, sparam *[]float32) {
	C.srotmg_((*C.float)(d1), (*C.float)(d2), (*C.float)(x), (*C.float)(y), (*C.float)(&(*sparam)[0]))
}

func snrm2Wrapper(n *int, sx *[]float32, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float32(C.snrm2_(&_n, (*C.float)(&(*sx)[0]), &_incx))
}

func snrm2WrapperTest(n *int, sx *[]float32, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float32(C.snrm2_(&_n, (*C.float)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func sasumWrapper(n *int, sx *[]float32, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float32(C.sasum_(&_n, (*C.float)(&(*sx)[0]), &_incx))
}

func sasumWrapperTest(n *int, sx *[]float32, incx *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := float32(C.sasum_(&_n, (*C.float)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func sscalWrapper(n *int, sa *float32, sx *[]float32, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.sscal_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx)
}

func sscalWrapperTest(n *int, sa *float32, sx *[]float32, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.sscal_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx)
	*n = int(_n)
	*incx = int(_incx)
}

func isamaxWrapper(n *int, sx *[]float32, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return int(C.isamax_(&_n, (*C.float)(&(*sx)[0]), &_incx))
}

func isamaxWrapperTest(n *int, sx *[]float32, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	res := int(C.isamax_(&_n, (*C.float)(&(*sx)[0]), &_incx))
	*n = int(_n)
	*incx = int(_incx)
	return res
}

func sdotWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float32(C.sdot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
}

func sdotWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := float32(C.sdot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func saxpyWrapper(n *int, sa *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.saxpy_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
}

func saxpyWrapperTest(n *int, sa *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.saxpy_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func scopyWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.scopy_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
}

func scopyWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.scopy_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sswapWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sswap_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
}

func sswapWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sswap_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func srotmWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, temp *[]float32) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.srotm_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy, (*C.float)(&(*temp)[0]))
}

func srotmWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, temp *[]float32) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.srotm_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy, (*C.float)(&(*temp)[0]))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sdsdotWrapper(n *int, sa *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float32(C.sdsdot_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
}

func sdsdotWrapperTest(n *int, sa *float32, sx *[]float32, incx *int, sy *[]float32, incy *int) float32 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	res := float32(C.sdsdot_(&_n, (*C.float)(sa), (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
	return res
}

func srotWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, sc, ss *float32) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.srot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy, (*C.float)(sc), (*C.float)(ss))
}

func srotWrapperTest(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int, sc, ss *float32) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.srot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy, (*C.float)(sc), (*C.float)(ss))
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sgemvWrapper(trans *byte, m, n *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sgemv_(&_trans, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
}

func sgemvWrapperTest(trans *byte, m, n *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sgemv_(&_trans, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sgbmvWrapper(trans *byte, m, n, kl, ku *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
}

func sgbmvWrapperTest(trans *byte, m, n, kl, ku *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_trans := C.char(*trans)
	_m := C.int(*m)
	_n := C.int(*n)
	_kl := C.int(*kl)
	_ku := C.int(*ku)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sgbmv_(&_trans, &_m, &_n, &_kl, &_ku, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
	*trans = byte(_trans)
	*m = int(_m)
	*n = int(_n)
	*kl = int(_kl)
	*ku = int(_ku)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ssymvWrapper(uplo *byte, n *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssymv_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
}

func ssymvWrapperTest(uplo *byte, n *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssymv_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ssbmvWrapper(uplo *byte, n, k *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssbmv_(&_uplo, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
}

func ssbmvWrapperTest(uplo *byte, n, k *int, alpha *float32, a *[]float32, lda *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssbmv_(&_uplo, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sspmvWrapper(uplo *byte, n *int, alpha *float32, a *[]float32, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sspmv_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
}

func sspmvWrapperTest(uplo *byte, n *int, alpha *float32, a *[]float32, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sspmv_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx, (*C.float)(beta), (*C.float)(&(*y)[0]), &_incy)
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func strmvWrapper(uplo, trans, diag *byte, n *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.strmv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
}

func strmvWrapperTest(uplo, trans, diag *byte, n *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.strmv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func stbmvWrapper(uplo, trans, diag *byte, n, k *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.stbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
}

func stbmvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.stbmv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func stpmvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.stpmv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx)
}

func stpmvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.stpmv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func strsvWrapper(uplo, trans, diag *byte, n *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.strsv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
}

func strsvWrapperTest(uplo, trans, diag *byte, n *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.strsv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func stbsvWrapper(uplo, trans, diag *byte, n, k *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.stbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
}

func stbsvWrapperTest(uplo, trans, diag *byte, n, k *int, a *[]float32, lda *int, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.stbsv_(&_uplo, &_trans, &_diag, &_n, &_k, (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*incx = int(_incx)
}

func stpsvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.stpsv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx)
}

func stpsvWrapperTest(uplo, trans, diag *byte, n *int, a, x *[]float32, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.stpsv_(&_uplo, &_trans, &_diag, &_n, (*C.float)(&(*a)[0]), (*C.float)(&(*x)[0]), &_incx)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*diag = byte(_diag)
	*n = int(_n)
	*incx = int(_incx)
}

func sgerWrapper(m, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sger_(&_m, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]), &_lda)
}

func sgerWrapperTest(m, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sger_(&_m, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]), &_lda)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func ssyrWrapper(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, a *[]float32, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ssyr_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*a)[0]), &_lda)
}

func ssyrWrapperTest(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, a *[]float32, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.ssyr_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
}

func ssprWrapper(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, a *[]float32) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.sspr_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*a)[0]))
}

func ssprWrapperTest(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, a *[]float32) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.sspr_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
}

func ssyr2Wrapper(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssyr2_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]), &_lda)
}

func ssyr2WrapperTest(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.ssyr2_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]), &_lda)
	*uplo = byte(_uplo)
	*n = int(_n)
	*lda = int(_lda)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sspr2Wrapper(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sspr2_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]))
}

func sspr2WrapperTest(uplo *byte, n *int, alpha *float32, x *[]float32, incx *int, y *[]float32, incy *int, a *[]float32) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.sspr2_(&_uplo, &_n, (*C.float)(alpha), (*C.float)(&(*x)[0]), &_incx, (*C.float)(&(*y)[0]), &_incy, (*C.float)(&(*a)[0]))
	*uplo = byte(_uplo)
	*n = int(_n)
	*incx = int(_incx)
	*incy = int(_incy)
}

func sgemmWrapper(transa, transb *byte, m, n, k *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.sgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
}

func sgemmWrapperTest(transa, transb *byte, m, n, k *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_transa := C.char(*transa)
	_transb := C.char(*transb)
	_m := C.int(*m)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.sgemm_(&_transa, &_transb, &_m, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
	*transa = byte(_transa)
	*transb = byte(_transb)
	*m = int(_m)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func ssymmWrapper(side, uplo *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.ssymm_(&_side, &_uplo, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
}

func ssymmWrapperTest(side, uplo *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.ssymm_(&_side, &_uplo, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}

func strmmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.strmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb)
}

func strmmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.strmm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func strsmWrapper(side, uplo, transa, diag *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.strsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb)
}

func strsmWrapperTest(side, uplo, transa, diag *byte, m, n *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int) {
	_side := C.char(*side)
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_diag := C.char(*diag)
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	C.strsm_(&_side, &_uplo, &_transa, &_diag, &_m, &_n, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb)
	*side = byte(_side)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*diag = byte(_diag)
	*m = int(_m)
	*n = int(_n)
	*lda = int(_lda)
	*ldb = int(_ldb)
}

func ssyrkWrapper(uplo, transa *byte, n, k *int, alpha *float32, a *[]float32, lda *int, beta *float32, c *[]float32, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.ssyrk_(&_uplo, &_transa, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
}

func ssyrkWrapperTest(uplo, transa *byte, n, k *int, alpha *float32, a *[]float32, lda *int, beta *float32, c *[]float32, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.ssyrk_(&_uplo, &_transa, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*transa = byte(_transa)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldc = int(_ldc)
}

func ssyr2kWrapper(uplo, trans *byte, n, k *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.ssyr2k_(&_uplo, &_trans, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
}

func ssyr2kWrapperTest(uplo, trans *byte, n, k *int, alpha *float32, a *[]float32, lda *int, b *[]float32, ldb *int, beta *float32, c *[]float32, ldc *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldb := C.int(*ldb)
	_ldc := C.int(*ldc)
	C.ssyr2k_(&_uplo, &_trans, &_n, &_k, (*C.float)(alpha), (*C.float)(&(*a)[0]), &_lda, (*C.float)(&(*b)[0]), &_ldb, (*C.float)(beta), (*C.float)(&(*c)[0]), &_ldc)
	*uplo = byte(_uplo)
	*trans = byte(_trans)
	*n = int(_n)
	*k = int(_k)
	*lda = int(_lda)
	*ldb = int(_ldb)
	*ldc = int(_ldc)
}
