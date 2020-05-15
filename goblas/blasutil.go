package goblas

import "C"

func ccompress(major byte, y *[][]complex64, nrows *int, ncols *int) *[]complex64 {
	var i, j int
	var _y []complex64
	_y = make([]complex64, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func ccompressC(major byte, y *[][]C.complexfloat, nrows *C.int, ncols *C.int) *[]C.complexfloat {
	var i, j C.int
	var _y []C.complexfloat
	_y = make([]C.complexfloat, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cexpand(major byte, y *[]complex64, nrows *int, ncols *int, inc *int) *[][]complex64 {
	var i, j int
	var _y [][]complex64
	if major == 'R' {
		_y = make([][]complex64, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]complex64, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex64, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]complex64, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cexpandC(major byte, y *[]C.complexfloat, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.complexfloat {
	var i, j C.int
	var _y [][]C.complexfloat
	if major == 'R' {
		_y = make([][]C.complexfloat, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.complexfloat, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexfloat, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.complexfloat, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cexpand2(major byte, y []complex64, nrows *int, ncols *int, inc *int) *[][]complex64 {
	var i, j int
	var _y [][]complex64

	if major == 'R' {
		_y = make([][]complex64, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]complex64, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex64, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]complex64, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cexpand2C(major byte, y []C.complexfloat, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.complexfloat {
	var i, j C.int
	var _y [][]C.complexfloat

	if major == 'R' {
		_y = make([][]C.complexfloat, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.complexfloat, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexfloat, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.complexfloat, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cslice(major byte, y *[]complex64, nrows *int, ncols *int, ioff *int) *[][]complex64 {
	var onei int
	var _x []complex64
	_y := new([][]complex64)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = cexpand(major, &_x, nrows, ncols, &onei)

	return _y
}

func csliceC(major byte, y *[]C.complexfloat, nrows *C.int, ncols *C.int, ioff *C.int) *[][]C.complexfloat {
	var onei C.int
	var _x []C.complexfloat
	_y := new([][]C.complexfloat)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = cexpandC(major, &_x, nrows, ncols, &onei)

	return _y
}

func cslice2(major byte, y *[][]complex64, nrows *int, ncols *int, irow *int, icol *int) *[][]complex64 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]complex64

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap((*y)[0])-irown)
	ncoln = min(*ncols, cap(*y)-icoln)

	if major == 'R' {
		_y = make([][]complex64, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]complex64, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex64, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]complex64, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cslice2C(major byte, y *[][]C.complexfloat, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.complexfloat {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.complexfloat

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap((*y)[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(*y))-icoln)

	if major == 'R' {
		_y = make([][]C.complexfloat, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.complexfloat, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexfloat, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.complexfloat, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func cslice3(major byte, y *[][]complex64, nrows *int, ncols *int, irow *int, icol *int) *[]complex64 {
	var i, j, ncolsn, nrowsn int
	var _x [][]complex64
	_y := new([]complex64)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = cap(*y) - ((*icol) - 1)
	nrowsn = cap((*y)[0]) - ((*irow) - 1)

	_x = make([][]complex64, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]complex64, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = ccompress(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func cslice3C(major byte, y *[][]C.complexfloat, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[]C.complexfloat {
	var i, j, ncolsn, nrowsn C.int
	var _x [][]C.complexfloat
	_y := new([]C.complexfloat)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = C.int(cap(*y)) - ((*icol) - 1)
	nrowsn = C.int(cap((*y)[0])) - ((*irow) - 1)

	_x = make([][]C.complexfloat, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]C.complexfloat, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = ccompressC(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func dcompress(major byte, y *[][]float64, nrows *int, ncols *int) *[]float64 {
	var i, j int
	var _y []float64
	_y = make([]float64, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dcompressC(major byte, y *[][]C.double, nrows *C.int, ncols *C.int) *[]C.double {
	var i, j C.int
	var _y []C.double
	_y = make([]C.double, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dexpand(major byte, y *[]float64, nrows *int, ncols *int, inc *int) *[][]float64 {
	var i, j int
	var _y [][]float64

	if major == 'R' {
		_y = make([][]float64, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]float64, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]float64, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]float64, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dexpandC(major byte, y *[]C.double, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.double {
	var i, j C.int
	var _y [][]C.double

	if major == 'R' {
		_y = make([][]C.double, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.double, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.double, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.double, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dexpand2(major byte, y []float64, nrows *int, ncols *int, inc *int) *[][]float64 {
	var i, j int
	var _y [][]float64

	if major == 'R' {
		_y = make([][]float64, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]float64, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]float64, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]float64, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dexpand2C(major byte, y []C.double, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.double {
	var i, j C.int
	var _y [][]C.double

	if major == 'R' {
		_y = make([][]C.double, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.double, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.double, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.double, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice(major byte, y *[]float64, nrows *int, ncols *int, ioff *int) *[][]float64 {
	var onei int
	var _x []float64
	_y := new([][]float64)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = dexpand(major, &_x, nrows, ncols, &onei)

	return _y
}

func dsliceC(major byte, y *[]C.double, nrows *C.int, ncols *C.int, ioff *C.int) *[][]C.double {
	var onei C.int
	var _x []C.double
	_y := new([][]C.double)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = dexpandC(major, &_x, nrows, ncols, &onei)

	return _y
}

func dslice2(major byte, y *[][]float64, nrows *int, ncols *int, irow *int, icol *int) *[][]float64 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]float64

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap((*y)[0])-irown)
	ncoln = min(*ncols, cap(*y)-icoln)

	if major == 'R' {
		_y = make([][]float64, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]float64, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]float64, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]float64, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice2C(major byte, y *[][]C.double, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.double {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.double

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap((*y)[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(*y))-icoln)

	if major == 'R' {
		_y = make([][]C.double, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.double, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.double, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.double, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice2d(major byte, y *[][]float64, nrows *int, ncols *int, irow *int, icol *int) *[][]float64 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]float64

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap((*y)[0])-irown)
	ncoln = min(*ncols, cap(*y)-icoln)

	if major == 'R' {
		_y = make([][]float64, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]float64, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]float64, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]float64, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice2dC(major byte, y *[][]C.double, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.double {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.double

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap((*y)[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(*y))-icoln)

	if major == 'R' {
		_y = make([][]C.double, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.double, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.double, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.double, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice2dLit(major byte, y [][]float64, nrows *int, ncols *int, irow *int, icol *int) *[][]float64 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]float64

	if *nrows == 0 || *ncols == 0 {
		return &y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap(y[0])-irown)
	ncoln = min(*ncols, cap(y)-icoln)

	if major == 'R' {
		_y = make([][]float64, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]float64, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = y[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]float64, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]float64, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = y[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice2dLitC(major byte, y [][]C.double, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.double {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.double

	if *nrows == 0 || *ncols == 0 {
		return &y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap(y[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(y))-icoln)

	if major == 'R' {
		_y = make([][]C.double, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.double, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = y[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.double, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.double, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = y[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func dslice3(major byte, y *[][]float64, nrows *int, ncols *int, irow *int, icol *int) *[]float64 {
	var i, j, ncolsn, nrowsn int
	var _x [][]float64
	_y := new([]float64)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = cap(*y) - ((*icol) - 1)
	nrowsn = cap((*y)[0]) - ((*irow) - 1)

	_x = make([][]float64, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]float64, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = dcompress(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func dslice3C(major byte, y *[][]C.double, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[]C.double {
	var i, j, ncolsn, nrowsn C.int
	var _x [][]C.double
	_y := new([]C.double)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = C.int(cap(*y)) - ((*icol) - 1)
	nrowsn = C.int(cap((*y)[0])) - ((*irow) - 1)

	_x = make([][]C.double, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]C.double, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = dcompressC(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func scompress(major byte, y *[][]float32, nrows *int, ncols *int) *[]float32 {
	var i, j int
	var _y []float32
	_y = make([]float32, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func scompressC(major byte, y *[][]C.float, nrows *C.int, ncols *C.int) *[]C.float {
	var i, j C.int
	var _y []C.float
	_y = make([]C.float, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sexpand(major byte, y *[]float32, nrows *int, ncols *int, inc *int) *[][]float32 {
	var i, j int
	var _y [][]float32
	if major == 'R' {
		_y = make([][]float32, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]float32, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]float32, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]float32, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sexpandC(major byte, y *[]C.float, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.float {
	var i, j C.int
	var _y [][]C.float
	if major == 'R' {
		_y = make([][]C.float, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.float, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.float, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.float, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sexpand2(major byte, y []float32, nrows *int, ncols *int, inc *int) *[][]float32 {
	var i, j int
	var _y [][]float32
	if major == 'R' {
		_y = make([][]float32, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]float32, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]float32, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]float32, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sexpand2C(major byte, y []C.float, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.float {
	var i, j C.int
	var _y [][]C.float
	if major == 'R' {
		_y = make([][]C.float, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.float, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.float, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.float, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sslice(major byte, y *[]float32, nrows *int, ncols *int, ioff *int) *[][]float32 {
	var onei int
	var _x []float32
	_y := new([][]float32)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = sexpand(major, &_x, nrows, ncols, &onei)

	return _y
}

func ssliceC(major byte, y *[]C.float, nrows *C.int, ncols *C.int, ioff *C.int) *[][]C.float {
	var onei C.int
	var _x []C.float
	_y := new([][]C.float)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = sexpandC(major, &_x, nrows, ncols, &onei)

	return _y
}

func sslice2(major byte, y *[][]float32, nrows *int, ncols *int, irow *int, icol *int) *[][]float32 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]float32

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap((*y)[0])-irown)
	ncoln = min(*ncols, cap(*y)-icoln)

	if major == 'R' {
		_y = make([][]float32, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]float32, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]float32, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]float32, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sslice2C(major byte, y *[][]C.float, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.float {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.float

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap((*y)[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(*y))-icoln)

	if major == 'R' {
		_y = make([][]C.float, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.float, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.float, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.float, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func sslice3(major byte, y *[][]float32, nrows *int, ncols *int, irow *int, icol *int) *[]float32 {
	var i, j, ncolsn, nrowsn int
	var _x [][]float32
	_y := new([]float32)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = cap(*y) - ((*icol) - 1)
	nrowsn = cap((*y)[0]) - ((*irow) - 1)

	_x = make([][]float32, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]float32, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = scompress(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func sslice3C(major byte, y *[][]C.float, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[]C.float {
	var i, j, ncolsn, nrowsn C.int
	var _x [][]C.float
	_y := new([]C.float)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = C.int(cap(*y)) - ((*icol) - 1)
	nrowsn = C.int(cap((*y)[0])) - ((*irow) - 1)

	_x = make([][]C.float, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]C.float, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = scompressC(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func zcompress(major byte, y *[][]complex128, nrows *int, ncols *int) *[]complex128 {
	var i, j int
	var _y []complex128
	_y = make([]complex128, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zcompressC(major byte, y *[][]C.complexdouble, nrows *C.int, ncols *C.int) *[]C.complexdouble {
	var i, j C.int
	var _y []C.complexdouble
	_y = make([]C.complexdouble, (*nrows)*(*ncols))
	if major == 'R' {
		for i = 0; i < *nrows; i++ {
			for j = 0; j < *ncols; j++ {
				_y[i*(*nrows)+j] = (*y)[i][j]
			}
		}
	} else if major == 'C' {
		for j = 0; j < *nrows; j++ {
			for i = 0; i < *ncols; i++ {
				_y[i+j*(*ncols)] = (*y)[i][j]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zexpand(major byte, y *[]complex128, nrows *int, ncols *int, inc *int) *[][]complex128 {
	var i, j int
	var _y [][]complex128
	if major == 'R' {
		_y = make([][]complex128, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]complex128, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex128, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]complex128, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zexpandC(major byte, y *[]C.complexdouble, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.complexdouble {
	var i, j C.int
	var _y [][]C.complexdouble
	if major == 'R' {
		_y = make([][]C.complexdouble, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.complexdouble, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = (*y)[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexdouble, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.complexdouble, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = (*y)[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zexpand2(major byte, y []complex128, nrows *int, ncols *int, inc *int) *[][]complex128 {
	var i, j int
	var _y [][]complex128

	if major == 'R' {
		_y = make([][]complex128, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]complex128, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex128, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]complex128, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zexpand2C(major byte, y []C.complexdouble, nrows *C.int, ncols *C.int, inc *C.int) *[][]C.complexdouble {
	var i, j C.int
	var _y [][]C.complexdouble

	if major == 'R' {
		_y = make([][]C.complexdouble, *nrows)
		for i = 0; i < *nrows; i++ {
			_y[i] = make([]C.complexdouble, *ncols)
			for j = 0; j < *ncols; j++ {
				_y[i][j] = y[i*(*nrows)+j*(*inc)]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexdouble, *ncols)
		for i = 0; i < *ncols; i++ {
			_y[i] = make([]C.complexdouble, *nrows)
			for j = 0; j < *nrows; j++ {
				_y[i][j] = y[i*(*inc)+j*(*ncols)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zslice(major byte, y *[]complex128, nrows *int, ncols *int, ioff *int) *[][]complex128 {
	var onei int
	var _x []complex128
	_y := new([][]complex128)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = zexpand(major, &_x, nrows, ncols, &onei)

	return _y
}

func zsliceC(major byte, y *[]C.complexdouble, nrows *C.int, ncols *C.int, ioff *C.int) *[][]C.complexdouble {
	var onei C.int
	var _x []C.complexdouble
	_y := new([][]C.complexdouble)

	onei = 1

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	_x = (*y)[(*ioff)-1:]
	_y = zexpandC(major, &_x, nrows, ncols, &onei)

	return _y
}

func zslice2(major byte, y *[][]complex128, nrows *int, ncols *int, irow *int, icol *int) *[][]complex128 {
	var i, j, irown, icoln, ncoln, nrown int
	var _y [][]complex128

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = min(*nrows, cap((*y)[0])-irown)
	ncoln = min(*ncols, cap(*y)-icoln)

	if major == 'R' {
		_y = make([][]complex128, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]complex128, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]complex128, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]complex128, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zslice2C(major byte, y *[][]C.complexdouble, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[][]C.complexdouble {
	var i, j, irown, icoln, ncoln, nrown C.int
	var _y [][]C.complexdouble

	if *nrows == 0 || *ncols == 0 {
		return y
	}

	irown = ((*irow) - 1)
	icoln = ((*icol) - 1)
	nrown = mincint(*nrows, C.int(cap((*y)[0]))-irown)
	ncoln = mincint(*ncols, C.int(cap(*y))-icoln)

	if major == 'R' {
		_y = make([][]C.complexdouble, (*nrows)-((*irow)-1))
		for i = 0; i < (*nrows)-((*irow)-1); i++ {
			_y[i] = make([]C.complexdouble, (*ncols)-((*icol)-1))
			for j = 0; j < (*ncols)-((*icol)-1); j++ {
				_y[i][j] = (*y)[i+(*irow)-1][j+(*icol)-1]
			}
		}
	} else if major == 'C' {
		_y = make([][]C.complexdouble, ncoln)
		for i = 0; i < ncoln; i++ {
			_y[i] = make([]C.complexdouble, nrown)
			for j = 0; j < nrown; j++ {
				_y[i][j] = (*y)[i+(icoln)][j+(irown)]
			}
		}
	} else {
		panic("Unrecognized major value")
	}

	return &_y
}

func zslice3(major byte, y *[][]complex128, nrows *int, ncols *int, irow *int, icol *int) *[]complex128 {
	var i, j, ncolsn, nrowsn int
	var _x [][]complex128
	_y := new([]complex128)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = cap(*y) - ((*icol) - 1)
	nrowsn = cap((*y)[0]) - ((*irow) - 1)

	_x = make([][]complex128, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]complex128, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = zcompress(major, &_x, &nrowsn, &ncolsn)

	return _y
}

func zslice3C(major byte, y *[][]C.complexdouble, nrows *C.int, ncols *C.int, irow *C.int, icol *C.int) *[]C.complexdouble {
	var i, j, ncolsn, nrowsn C.int
	var _x [][]C.complexdouble
	_y := new([]C.complexdouble)

	if *nrows == 0 || *ncols == 0 {
		return _y
	}

	ncolsn = C.int(cap(*y)) - ((*icol) - 1)
	nrowsn = C.int(cap((*y)[0])) - ((*irow) - 1)

	_x = make([][]C.complexdouble, ncolsn)
	for i = 0; i < ncolsn; i++ {
		_x[i] = make([]C.complexdouble, nrowsn)
		for j = 0; j < nrowsn; j++ {
			_x[i][j] = (*y)[i+(*icol)-1][j+(*irow)-1]
		}
	}
	_y = zcompressC(major, &_x, &nrowsn, &ncolsn)

	return _y
}
