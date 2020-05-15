package goblas

type memory struct {
	combla struct {
		icase *int
		n     *int
		incx  *int
		incy  *int
		mode  *int
		pass  *bool
	}
	infoc struct {
		infot *int
		nout  *int
		ok    *bool
		lerr  *bool
		noutc *int
	}
	srnamc struct {
		srnamt *[]byte
	}
}

var common memory
