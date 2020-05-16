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

func dasumWrapper(n *int, sx *[]float64, incx *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return float64(C.dasum_(&_n, (*C.double)(&(*sx)[0]), &_incx))
}

func dscalWrapper(n *int, sa *float64, sx *[]float64, incx *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dscal_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx)
}

func idamaxWrapper(n *int, sx *[]float64, incx *int) int {
	_n := C.int(*n)
	_incx := C.int(*incx)
	return int(C.idamax_(&_n, (*C.double)(&(*sx)[0]), &_incx))
}

func ddotWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float64(C.ddot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy))
}

func daxpyWrapper(n *int, sa *float64, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.daxpy_(&_n, (*C.double)(sa), (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func dcopyWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dcopy_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func dswapWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dswap_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy)
}

func drotmWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, temp *[]float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drotm_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(&(*temp)[0]))
}

func dsdotWrapper(n *int, sx *[]float32, incx *int, sy *[]float32, incy *int) float64 {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	return float64(C.dsdot_(&_n, (*C.float)(&(*sx)[0]), &_incx, (*C.float)(&(*sy)[0]), &_incy))
}

func drotWrapper(n *int, sx *[]float64, incx *int, sy *[]float64, incy *int, sc, ss *float64) {
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.drot_(&_n, (*C.double)(&(*sx)[0]), &_incx, (*C.double)(&(*sy)[0]), &_incy, (*C.double)(sc), (*C.double)(ss))
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

func dsymvWrapper(uplo *byte, n *int, alpha *float64, a *[]float64, lda *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsymv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
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

func dspmvWrapper(uplo *byte, n *int, alpha *float64, a *[]float64, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspmv_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx, (*C.double)(beta), (*C.double)(&(*y)[0]), &_incy)
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

func dtpmvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpmv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
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

func dtpsvWrapper(uplo, trans, diag *byte, n *int, a, x *[]float64, incx *int) {
	_uplo := C.char(*uplo)
	_trans := C.char(*trans)
	_diag := C.char(*diag)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dtpsv_(&_uplo, &_trans, &_diag, &_n, (*C.double)(&(*a)[0]), (*C.double)(&(*x)[0]), &_incx)
}

func dgerWrapper(m, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_m := C.int(*m)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dger_(&_m, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
}

func dsyrWrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	C.dsyr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]), &_lda)
}

func dsprWrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	C.dspr_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*a)[0]))
}

func dsyr2Wrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64, lda *int) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_lda := C.int(*lda)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dsyr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]), &_lda)
}

func dspr2Wrapper(uplo *byte, n *int, alpha *float64, x *[]float64, incx *int, y *[]float64, incy *int, a *[]float64) {
	_uplo := C.char(*uplo)
	_n := C.int(*n)
	_incx := C.int(*incx)
	_incy := C.int(*incy)
	C.dspr2_(&_uplo, &_n, (*C.double)(alpha), (*C.double)(&(*x)[0]), &_incx, (*C.double)(&(*y)[0]), &_incy, (*C.double)(&(*a)[0]))
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

func dsyrkWrapper(uplo, transa *byte, n, k *int, alpha *float64, a *[]float64, lda *int, beta *float64, c *[]float64, ldc *int) {
	_uplo := C.char(*uplo)
	_transa := C.char(*transa)
	_n := C.int(*n)
	_k := C.int(*k)
	_lda := C.int(*lda)
	_ldc := C.int(*ldc)
	C.dsyrk_(&_uplo, &_transa, &_n, &_k, (*C.double)(alpha), (*C.double)(&(*a)[0]), &_lda, (*C.double)(beta), (*C.double)(&(*c)[0]), &_ldc)
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
