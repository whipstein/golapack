package goblas

import "C"
import (
	"fmt"
	"testing"
)

type memory struct {
	combla struct {
		icase int
		nout  int
		ntra  int
		n     int
		incx  int
		incy  int
		mode  int
		pass  bool
	}
	infoc struct {
		infot int
		infox int
		nout  int
		ntra  int
		ndata int
		ok    bool
		lerr  bool
		noutc int
		test  bool
	}
	srnamc struct {
		srnamt string
		srnamx string
	}
	beg struct {
		mi int
		mj int
		ic int
		i  int
		j  int
	}
}

var common memory

func cbeg(reset *bool) complex64 {
	var cbegReturn complex64
	var i, j, ic, mi, mj *int
	//
	//  Generates complex numbers as pairs of random numbers uniformly
	//  distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	mj = &common.beg.mj
	ic = &common.beg.ic
	i = &common.beg.i
	j = &common.beg.j

	if *(reset) {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i or j is bounded between 1 and 999.
	//     If initial i or j = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i or j = 4 or 8, the period will be 25.
	//     If initial i or j = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i or j
	//     in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*j *= *mj
	*i -= 1000 * ((*i) / 1000)
	*j -= 1000 * ((*j) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	cbegReturn = cmplxc64(float32((*i)-500)/1001.0, float32((*j)-500)/1001.0)
	return cbegReturn
}

func cbegC(reset *bool) C.complexfloat {
	var cbegReturn C.complexfloat
	var i, j, ic, mi, mj *int
	//
	//  Generates complex numbers as pairs of random numbers uniformly
	//  distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	mj = &common.beg.mj
	ic = &common.beg.ic
	i = &common.beg.i
	j = &common.beg.j

	if *(reset) {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i or j is bounded between 1 and 999.
	//     If initial i or j = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i or j = 4 or 8, the period will be 25.
	//     If initial i or j = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i or j
	//     in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*j *= *mj
	*i -= 1000 * ((*i) / 1000)
	*j -= 1000 * ((*j) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	cbegReturn = cmplxcfloat(C.float((*i)-500)/1001.0, C.float((*j)-500)/1001.0)
	return cbegReturn
}

func chkerr(err error) {
	if err != nil {
		_ = fmt.Errorf("%v", err)
	}
}

func chkxer(t *testing.T) {
	var infot, infox *int
	var srnamt, srnamx *string
	//
	//  Tests whether XERBLA has detected an error when it should.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	infot = &common.infoc.infot
	infox = &common.infoc.infox
	srnamt = &common.srnamc.srnamt
	srnamx = &common.srnamc.srnamx

	if *srnamt != *srnamx {
		t.Errorf("Error Exit Test Failed: exit code: {%s} output, {%s} expected", *srnamx, *srnamt)
	}
	if *infot != *infox {
		t.Errorf("Error Exit Test Failed: %s: {%d} output, {%d} expected", *srnamt, *infox, *infot)
	}
	return
}

func cmakeL2(_type string, uplo byte, diag byte, m *int, n *int, a *[][]complex64, nmax *int, aa *[]complex64, lda *int, kl *int, ku *int, reset *bool, transl *complex64) {
	var zero, one, rogue complex64
	var rzero, rrogue float32
	var i, i1, i2, i3, ibeg, iend, ioff, j, jj, kk int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0 + (0.0)*1i)
	one = (1.0 + (0.0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0
	rrogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'H'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbeg(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = conjgc64((*a)[i-1][j-1])
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if sym {
			(*a)[j-1][j-1] = cmplxc64(realc64((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-2][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc64(realc64((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = max(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = min((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = kk + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc64(realc64((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
					if sym {
						(*aa)[ioff-1] = cmplxc64(realc64((*aa)[ioff-1]), rrogue)
					}
				}
			}
		}
	}
	return
}

func cmakeL2C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.complexfloat, nmax *C.int, aa *[]C.complexfloat, lda *C.int, kl *C.int, ku *C.int, reset *bool, transl *C.complexfloat) {
	var zero, one, rogue C.complexfloat
	var rzero, rrogue C.float
	var i, i1, i2, i3, ibeg, iend, ioff, j, jj, kk C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0 + (0.0)*1i)
	one = (1.0 + (0.0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0
	rrogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'H'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = cbegC(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = conjgcfloat((*a)[i-1][j-1])
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if sym {
			(*a)[j-1][j-1] = cmplxcfloat(realcfloat((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= mincint((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-2][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcfloat(realcfloat((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = maxcint(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = mincint((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = kk + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcfloat(realcfloat((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
					if sym {
						(*aa)[ioff-1] = cmplxcfloat(realcfloat((*aa)[ioff-1]), rrogue)
					}
				}
			}
		}
	}
	return
}

func cmakeL3(_type string, uplo byte, diag byte, m *int, n *int, a *[][]complex64, nmax *int, aa *[]complex64, lda *int, reset *bool, transl *complex64) {
	var zero, one, rogue complex64
	var rzero, rrogue float32
	var i, ibeg, iend, j, jj int
	var gen, her, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'HE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0 + (0.0)*1i)
	one = (1.0 + (0.0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0
	rrogue = -1.0e10

	gen = _type == "GE"
	her = _type == "HE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (her || sym || tri) && uplo == 'U'
	lower = (her || sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = cbeg(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if her {
						(*a)[j-1][i-1] = conjgc64((*a)[i-1][j-1])
					} else if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if her {
			(*a)[j-1][j-1] = cmplxc64(realc64((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if her {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc64(realc64((*aa)[jj-1]), rrogue)
			}
		}
	}
	return
}

func cmakeL3C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.complexfloat, nmax *C.int, aa *[]C.complexfloat, lda *C.int, reset *bool, transl *C.complexfloat) {
	var zero, one, rogue C.complexfloat
	var rzero, rrogue C.float
	var i, ibeg, iend, j, jj C.int
	var gen, her, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'HE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0 + (0.0)*1i)
	one = (1.0 + (0.0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0
	rrogue = -1.0e10

	gen = _type == "GE"
	her = _type == "HE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (her || sym || tri) && uplo == 'U'
	lower = (her || sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = cbegC(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if her {
						(*a)[j-1][i-1] = conjgcfloat((*a)[i-1][j-1])
					} else if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if her {
			(*a)[j-1][j-1] = cmplxcfloat(realcfloat((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if her {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcfloat(realcfloat((*aa)[jj-1]), rrogue)
			}
		}
	}
	return
}

func cmmch(transa byte, transb byte, m *int, n *int, kk *int, alpha *complex64, a *[][]complex64, lda *int, b *[][]complex64, ldb *int, beta *complex64, c *[][]complex64, ldc *int, ct *[]complex64, g *[]float32, cc *[][]complex64, LDCC *int, eps *float32, err *float32, fatal *bool, nout *int, mv bool) {
	var zero complex64
	var rzero, rone float32
	var erri float32
	var i, j, k int
	var ctrana, ctranb, trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0 + (0.0)*1i)
	rzero = 0.0
	rone = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	ctrana = transa == 'C'
	ctranb = transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = rzero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*(*b)[k-1][j-1]
					(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[i-1][k-1])*abssumf32((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			if ctrana {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + conjgc64((*a)[k-1][i-1])*(*b)[k-1][j-1]
						(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*(*b)[k-1][j-1]
						(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[k-1][j-1])
					}
				}
			}
		} else if !trana && tranb {
			if ctranb {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*conjgc64((*b)[j-1][k-1])
						(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[i-1][k-1])*abssumf32((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*(*b)[j-1][k-1]
						(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[i-1][k-1])*abssumf32((*b)[j-1][k-1])
					}
				}
			}
		} else if trana && tranb {
			if ctrana {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + conjgc64((*a)[k-1][i-1])*conjgc64((*b)[j-1][k-1])
							(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + conjgc64((*a)[k-1][i-1])*(*b)[j-1][k-1]
							(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[j-1][k-1])
						}
					}
				}
			} else {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*conjgc64((*b)[j-1][k-1])
							(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*(*b)[j-1][k-1]
							(*g)[i-1] = (*g)[i-1] + abssumf32((*a)[k-1][i-1])*abssumf32((*b)[j-1][k-1])
						}
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumf32(*alpha)*(*g)[i-1] + abssumf32(*beta)*abssumf32((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = rzero
		for i = 1; i <= *m; i++ {
			erri = abssumf32((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != rzero {
				erri /= (*g)[i-1]
			}
			*err = maxf32(*err, erri)
			if (*err)*sqrtf32(*eps) >= rone {
				goto Label230
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label250
	//
	//     Report fatal error.
	//
Label230:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label250:
	;
	return
}

func cmmchC(transa byte, transb byte, m *C.int, n *C.int, kk *C.int, alpha *C.complexfloat, a *[][]C.complexfloat, lda *C.int, b *[][]C.complexfloat, ldb *C.int, beta *C.complexfloat, c *[][]C.complexfloat, ldc *C.int, ct *[]C.complexfloat, g *[]C.float, cc *[][]C.complexfloat, LDCC *C.int, eps *C.float, err *C.float, fatal *bool, nout *int, mv bool) {
	var zero C.complexfloat
	var rzero, rone C.float
	// var cl C.complexfloat
	var erri C.float
	var i, j, k C.int
	var ctrana, ctranb, trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0 + (0.0)*1i)
	rzero = 0.0
	rone = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	ctrana = transa == 'C'
	ctranb = transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = rzero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*(*b)[k-1][j-1]
					(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[i-1][k-1])*abssumcfloat((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			if ctrana {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + conjgcfloat((*a)[k-1][i-1])*(*b)[k-1][j-1]
						(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*(*b)[k-1][j-1]
						(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[k-1][j-1])
					}
				}
			}
		} else if !trana && tranb {
			if ctranb {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*conjgcfloat((*b)[j-1][k-1])
						(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[i-1][k-1])*abssumcfloat((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] = (*ct)[i-1] + (*a)[i-1][k-1]*(*b)[j-1][k-1]
						(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[i-1][k-1])*abssumcfloat((*b)[j-1][k-1])
					}
				}
			}
		} else if trana && tranb {
			if ctrana {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + conjgcfloat((*a)[k-1][i-1])*conjgcfloat((*b)[j-1][k-1])
							(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + conjgcfloat((*a)[k-1][i-1])*(*b)[j-1][k-1]
							(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[j-1][k-1])
						}
					}
				}
			} else {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*conjgcfloat((*b)[j-1][k-1])
							(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] = (*ct)[i-1] + (*a)[k-1][i-1]*(*b)[j-1][k-1]
							(*g)[i-1] = (*g)[i-1] + abssumcfloat((*a)[k-1][i-1])*abssumcfloat((*b)[j-1][k-1])
						}
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumcfloat(*alpha)*(*g)[i-1] + abssumcfloat(*beta)*abssumcfloat((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = rzero
		for i = 1; i <= *m; i++ {
			erri = abssumcfloat((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != rzero {
				erri /= (*g)[i-1]
			}
			*err = maxcfloat(*err, erri)
			if (*err)*sqrtcfloat(*eps) >= rone {
				goto Label230
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label250
	//
	//     Report fatal error.
	//
Label230:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label250:
	;
	return
}

func cmvch(trans byte, m *int, n *int, alpha *complex64, a *[][]complex64, nmax *int, x *[]complex64, incx *int, beta *complex64, y *[]complex64, incy *int, yt *[]complex64, g *[]float32, yy *[]complex64, eps *float32, err *float32, fatal *bool, nout *int, mv bool) {
	var zero complex64
	var rzero, rone float32
	var erri float32
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl int
	var ctran, tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0 + (0.0)*1i)
	rzero = 0.0
	rone = 1.0

	tran = trans == 'T'
	ctran = trans == 'C'
	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = rzero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[j-1][i-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += conjgc64((*a)[j-1][i-1]) * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[j-1][i-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf32((*a)[i-1][j-1]) * abssumf32((*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumf32(*alpha)*(*g)[iy-1] + abssumf32(*beta)*abssumf32((*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = rzero
	for i = 1; i <= ml; i++ {
		erri = absc64((*yt)[i-1]-(*yy)[1+(i-1)*absint(*incy)-1]) / (*eps)
		if (*g)[i-1] != rzero {
			erri = erri / (*g)[i-1]
		}
		*err = maxf32(*err, erri)
		if (*err)*sqrtf32(*eps) >= rone {
			goto Label60
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label80
	//
	//     Report fatal error.
	//
Label60:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN half ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d  (%15.6f,%15.6f)(%15.6f,%15.6f)\n", i, (*yt)[i-1], (*yy)[1+(i-1)*absint(*incy)-1])
		} else {
			writeString(*nout, " %7d  (%15.6f,%15.6f)(%15.6f,%15.6f)\n", i, (*yy)[1+(i-1)*absint(*incy)-1], (*yt)[i-1])
		}
	}

Label80:
	;
	return
}

func cmvchC(trans byte, m *C.int, n *C.int, alpha *C.complexfloat, a *[][]C.complexfloat, nmax *C.int, x *[]C.complexfloat, incx *C.int, beta *C.complexfloat, y *[]C.complexfloat, incy *C.int, yt *[]C.complexfloat, g *[]C.float, yy *[]C.complexfloat, eps *C.float, err *C.float, fatal *bool, nout *int, mv bool) {
	var zero C.complexfloat
	var rzero, rone C.float
	var erri C.float
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl C.int
	var ctran, tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0 + (0.0)*1i)
	rzero = 0.0
	rone = 1.0

	tran = trans == 'T'
	ctran = trans == 'C'
	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = rzero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumcfloat((*a)[j-1][i-1]) * abssumcfloat((*x)[jx-1])
				jx += incxl
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += conjgcfloat((*a)[j-1][i-1]) * (*x)[jx-1]
				(*g)[iy-1] += abssumcfloat((*a)[j-1][i-1]) * abssumcfloat((*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumcfloat((*a)[i-1][j-1]) * abssumcfloat((*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumcfloat(*alpha)*(*g)[iy-1] + abssumcfloat(*beta)*abssumcfloat((*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = rzero
	for i = 1; i <= ml; i++ {
		erri = absccfloat((*yt)[i-1]-(*yy)[1+(i-1)*abscint(*incy)-1]) / (*eps)
		if (*g)[i-1] != rzero {
			erri = erri / (*g)[i-1]
		}
		*err = maxcfloat(*err, erri)
		if (*err)*sqrtcfloat(*eps) >= rone {
			goto Label60
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label80
	//
	//     Report fatal error.
	//
Label60:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d  (%15.6f)(%15.6f)\n", i, (*yt)[i-1], (*yy)[1+(i-1)*abscint(*incy)-1])
		} else {
			writeString(*nout, " %7d  (%15.6f)(%15.6f)\n", i, (*yy)[1+(i-1)*abscint(*incy)-1], (*yt)[i-1])
		}
	}

Label80:
	;
	return
}

func ctest(len int, ccomp *[]complex64, ctrue *[]complex64, csize *[]complex64, sfac *float32) {
	var i int
	scomp := make([]float32, 20)
	ssize := make([]float32, 20)
	strue := make([]float32, 20)
	//     **************************** CTEST *****************************
	//
	//     c.l. LAWSON, JPL, 1978 DEC 6
	//
	for i = 1; i <= len; i++ {
		scomp[2*i-0] = real((*ccomp)[i-1])
		scomp[2*i-1] = imag((*ccomp)[i-1])
		strue[2*i-0] = real((*ctrue)[i-1])
		strue[2*i-1] = imag((*ctrue)[i-1])
		ssize[2*i-0] = real((*csize)[i-1])
		ssize[2*i-1] = imag((*csize)[i-1])
	}

	stest(2*len, &scomp, &strue, &ssize, sfac)
	return
}

func ctestC(len C.int, ccomp *[]C.complexfloat, ctrue *[]C.complexfloat, csize *[]C.complexfloat, sfac *C.float) {
	var i C.int
	scomp := make([]C.float, 20)
	ssize := make([]C.float, 20)
	strue := make([]C.float, 20)
	//     **************************** CTEST *****************************
	//
	//     c.l. LAWSON, JPL, 1978 DEC 6
	//
	for i = 1; i <= len; i++ {
		scomp[2*i-0] = realcfloat((*ccomp)[i-1])
		scomp[2*i-1] = imagcfloat((*ccomp)[i-1])
		strue[2*i-0] = realcfloat((*ctrue)[i-1])
		strue[2*i-1] = imagcfloat((*ctrue)[i-1])
		ssize[2*i-0] = realcfloat((*csize)[i-1])
		ssize[2*i-1] = imagcfloat((*csize)[i-1])
	}

	stestC(2*len, &scomp, &strue, &ssize, sfac)
	return
}

func dbeg(reset *bool) float64 {
	var dbegReturn float64
	var i, ic, mi *int
	//
	//  Generates random numbers uniformly distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	ic = &common.beg.ic
	i = &common.beg.i

	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i is bounded between 1 and 999.
	//     If initial i = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i = 4 or 8, the period will be 25.
	//     If initial i = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*i -= 1000 * ((*i) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	dbegReturn = float64((*i)-500) / 1001.0
	return dbegReturn
}

func dbegC(reset *bool) C.double {
	var dbegReturn C.double
	var i, ic, mi *int
	//
	//  Generates random numbers uniformly distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	ic = &common.beg.ic
	i = &common.beg.i

	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i is bounded between 1 and 999.
	//     If initial i = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i = 4 or 8, the period will be 25.
	//     If initial i = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*i -= 1000 * ((*i) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	dbegReturn = C.double((*i)-500) / 1001.0
	return dbegReturn
}

func ddiff(da *float64, db *float64) float64 {
	//     ********************************* ddiff **************************
	//     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. l. LAWSON, JPL 1974 FEB 15
	//
	return (*da) - (*db)
}

func ddiffC(da *C.double, db *C.double) C.double {
	//     ********************************* ddiff **************************
	//     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. l. LAWSON, JPL 1974 FEB 15
	//
	return (*da) - (*db)
}

func dmakeL2(_type string, uplo byte, diag byte, m *int, n *int, a *[][]float64, nmax *int, aa *[]float64, lda *int, kl *int, ku *int, reset *bool, transl *float64) {
	var zero, one, rogue float64
	var i, i1, i2, i3, ibeg, iend, ioff, j, kk int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'S'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbeg(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-1-1][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = max(1, (*kl)+2-j)
				if unit {
					iend = (*kl)
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = min((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
				}
			}
		}
	}
	return
}

func dmakeL2C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.double, nmax *C.int, aa *[]C.double, lda *C.int, kl *C.int, ku *C.int, reset *bool, transl *C.double) {
	var zero, one, rogue C.double
	var i, i1, i2, i3, ibeg, iend, ioff, j, kk C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'S'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = dbegC(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= mincint((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-1-1][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = maxcint(1, (*kl)+2-j)
				if unit {
					iend = (*kl)
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = mincint((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
				}
			}
		}
	}
	return
}

func dmakeL3(_type string, uplo byte, diag byte, m *int, n *int, a *[][]float64, nmax *int, aa *[]float64, lda *int, reset *bool, transl *float64) {
	var zero, one, rogue float64
	var i, ibeg, iend, j int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type == "GE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = dbeg(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	}
	return
}

func dmakeL3C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.double, nmax *C.int, aa *[]C.double, lda *C.int, reset *bool, transl *C.double) {
	var zero, one, rogue C.double
	var i, ibeg, iend, j C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type == "GE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = dbegC(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	}
	return
}

func dmmch(transa byte, transb byte, m *int, n *int, kk *int, alpha *float64, a *[][]float64, lda *int, b *[][]float64, ldb *int, beta *float64, c *[][]float64, ldc *int, ct *[]float64, g *[]float64, cc *[][]float64, LDCC *int, eps *float64, err *float64, fatal *bool, nout *int, mv bool) {
	var zero, one, erri float64
	var i, j, k int
	var trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = zero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf64((*a)[i-1][k-1]) * absf64((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf64((*a)[k-1][i-1]) * absf64((*b)[k-1][j-1])
				}
			}
		} else if !trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf64((*a)[i-1][k-1]) * absf64((*b)[j-1][k-1])
				}
			}
		} else if trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf64((*a)[k-1][i-1]) * absf64((*b)[j-1][k-1])
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = absf64(*alpha)*(*g)[i-1] + absf64(*beta)*absf64((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		(*err) = zero
		for i = 1; i <= *m; i++ {
			erri = absf64((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != zero {
				erri /= (*g)[i-1]
			}
			*err = maxf64(*err, erri)
			if (*err)*sqrtf64(*eps) >= one {
				goto Label130
			}
		}

	}

	//     If the loop completes, all results are at least half accurate.
	goto Label150
	//
	//     Report fatal error.
	//
Label130:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label150:
	;
	return
}

func dmmchC(transa byte, transb byte, m *C.int, n *C.int, kk *C.int, alpha *C.double, a *[][]C.double, lda *C.int, b *[][]C.double, ldb *C.int, beta *C.double, c *[][]C.double, ldc *C.int, ct *[]C.double, g *[]C.double, cc *[][]C.double, LDCC *C.int, eps *C.double, err *C.double, fatal *bool, nout *int, mv bool) {
	var zero, one, erri C.double
	var i, j, k C.int
	var trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = zero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abscdouble((*a)[i-1][k-1]) * abscdouble((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abscdouble((*a)[k-1][i-1]) * abscdouble((*b)[k-1][j-1])
				}
			}
		} else if !trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += abscdouble((*a)[i-1][k-1]) * abscdouble((*b)[j-1][k-1])
				}
			}
		} else if trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += abscdouble((*a)[k-1][i-1]) * abscdouble((*b)[j-1][k-1])
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abscdouble(*alpha)*(*g)[i-1] + abscdouble(*beta)*abscdouble((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		(*err) = zero
		for i = 1; i <= *m; i++ {
			erri = abscdouble((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != zero {
				erri /= (*g)[i-1]
			}
			*err = maxcdouble(*err, erri)
			if (*err)*sqrtcdouble(*eps) >= one {
				goto Label130
			}
		}

	}

	//     If the loop completes, all results are at least half accurate.
	goto Label150
	//
	//     Report fatal error.
	//
Label130:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label150:
	;
	return
}

func dmvch(trans byte, m *int, n *int, alpha *float64, a *[][]float64, nmax *int, x *[]float64, incx *int, beta *float64, y *[]float64, incy *int, yt *[]float64, g *[]float64, yy *[]float64, eps *float64, err *float64, fatal *bool, nout *int, mv bool) {
	var zero, one, erri float64
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl int
	var tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0

	tran = trans == 'T' || trans == 'C'
	if tran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = zero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += absf64((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += absf64((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = absf64((*alpha))*(*g)[iy-1] + absf64((*beta)*(*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = zero
	for i = 1; i <= ml; i++ {
		erri = absf64((*yt)[i-1]-(*yy)[1+(i-1)*absint(*incy)-1]) / (*eps)
		if (*g)[i-1] != zero {
			erri /= (*g)[i-1]
		}
		*err = maxf64(*err, erri)
		if (*err)*sqrtf64(*eps) >= one {
			goto Label50
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label70
	//
	//     Report fatal error.
	//
Label50:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*absint(*incy)-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yy)[1+(i-1)*absint(*incy)-1], (*yt)[i-1])
		}
	}

Label70:
	;
	return
}

func dmvchC(trans byte, m *C.int, n *C.int, alpha *C.double, a *[][]C.double, nmax *C.int, x *[]C.double, incx *C.int, beta *C.double, y *[]C.double, incy *C.int, yt *[]C.double, g *[]C.double, yy *[]C.double, eps *C.double, err *C.double, fatal *bool, nout *int, mv bool) {
	var zero, one, erri C.double
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl C.int
	var tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0

	tran = trans == 'T' || trans == 'C'
	if tran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = zero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abscdouble((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abscdouble((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abscdouble((*alpha))*(*g)[iy-1] + abscdouble((*beta)*(*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = zero
	for i = 1; i <= ml; i++ {
		erri = abscdouble((*yt)[i-1]-(*yy)[1+(i-1)*abscint(*incy)-1]) / (*eps)
		if (*g)[i-1] != zero {
			erri /= (*g)[i-1]
		}
		*err = maxcdouble(*err, erri)
		if (*err)*sqrtcdouble(*eps) >= one {
			goto Label50
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label70
	//
	//     Report fatal error.
	//
Label50:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*abscint(*incy)-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yy)[1+(i-1)*abscint(*incy)-1], (*yt)[i-1])
		}
	}

Label70:
	;
	return
}

func dtest(len int, dcomp *[]float64, dtrue *[]float64, dsize *[]float64, dfac *float64) {
	var incx, incy, nout, icase *int
	var n *int
	var pass *bool
	var sd float64
	var i int
	//     ********************************* dtest **************************
	//
	//     THIS SUBR COMPARES ARRAYS  dcomp() AND dtrue() OF LENGTH _len TO
	//     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY dfac, ARE
	//     NEGLIGIBLE.
	//
	//     C. l. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	for i = 1; i <= len; i++ {
		sd = (*dcomp)[i-1] - (*dtrue)[i-1]
		if absf64((*dfac)*sd) <= absf64((*dsize)[i-1])*epsilonf64()+1e-18 {
			goto Label40
		}
		//
		//                             HERE    dcomp(i) IS NOT CLOSE TO dtrue(i).
		//
		if !*pass {
			goto Label20
		}
		//                             PRINT FAIL MESSAGE AND header.
		*pass = false
		writeString(*nout, "                                       FAIL\n")
		writeString(*nout, "\n CASE  N INCX INCY  I                             COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n \n")
	Label20:
		;
		writeString(*nout, " %4d%3d%5d%5d%3d%36.16e%36.16e%12.4e%12.4e\n", *icase, *n, *incx, *incy, i, (*dcomp)[i-1], (*dtrue)[i-1], sd, (*dsize)[i-1])
	Label40:
	}
	return
}

func dtestC(len C.int, dcomp *[]C.double, dtrue *[]C.double, dsize *[]C.double, dfac *C.double) {
	var incx, incy, nout, icase *int
	var n *int
	var pass *bool
	var sd C.double
	var i C.int
	//     ********************************* dtest **************************
	//
	//     THIS SUBR COMPARES ARRAYS  dcomp() AND dtrue() OF LENGTH _len TO
	//     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY dfac, ARE
	//     NEGLIGIBLE.
	//
	//     C. l. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	for i = 1; i <= len; i++ {
		sd = (*dcomp)[i-1] - (*dtrue)[i-1]
		if abscdouble((*dfac)*sd) <= abscdouble((*dsize)[i-1])*epsiloncdouble()+1e-18 {
			goto Label40
		}
		//
		//                             HERE    dcomp(i) IS NOT CLOSE TO dtrue(i).
		//
		if !*pass {
			goto Label20
		}
		//                             PRINT FAIL MESSAGE AND header.
		*pass = false
		writeString(*nout, "                                       FAIL\n")
		writeString(*nout, "\n CASE  N INCX INCY  I                             COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n \n")
	Label20:
		;
		writeString(*nout, " %4d%3d%5d%5d%3d%36.16e%36.16e%12.4e%12.4e\n", *icase, *n, *incx, *incy, i, (*dcomp)[i-1], (*dtrue)[i-1], sd, (*dsize)[i-1])
	Label40:
	}
	return
}

func dtest1(dcomp1 *float64, dtrue1 *float64, dsize *[]float64, dfac *float64) {
	dcomp := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	dtrue := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	//     ************************* dtest1 *****************************
	//
	//     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//
	//     C.l. LAWSON, JPL, 1978 DEC 6
	//
	(*dcomp)[0] = *dcomp1
	(*dtrue)[0] = *dtrue1
	dtest(1, dcomp, dtrue, dsize, dfac)

	return
}

func dtest1C(dcomp1 *C.double, dtrue1 *C.double, dsize *[]C.double, dfac *C.double) {
	dcomp := func() *[]C.double {
		arr := make([]C.double, 1)
		return &arr
	}()
	dtrue := func() *[]C.double {
		arr := make([]C.double, 1)
		return &arr
	}()
	//     ************************* dtest1 *****************************
	//
	//     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//
	//     C.l. LAWSON, JPL, 1978 DEC 6
	//
	(*dcomp)[0] = *dcomp1
	(*dtrue)[0] = *dtrue1
	dtestC(1, dcomp, dtrue, dsize, dfac)

	return
}

func itest1(icomp *int, itrue *int) {
	var nout, icase *int
	var incx, incy, n *int
	var pass *bool
	var id int
	//     ********************************* ITEST1 *************************
	//
	//     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
	//     EQUALITY.
	//     C. L. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	if *icomp == *itrue {
		goto Label40
	}
	//
	//                            HERE ICOMP IS NOT EQUAL TO ITRUE.
	//
	if !*pass {
		goto Label20
	}
	//                             PRINT FAIL MESSAGE AND HEADER.
	*pass = false
	writeString(*nout, "                                       FAIL\n")
	writeString(*nout, "\n CASE  N INCX INCY                                COMP                                TRUE     DIFFERENCE\n \n")
Label20:
	;
	id = (*icomp) - (*itrue)
	writeString(*nout, " %4d%3d%5d%5d%36d%36d%12d%12d\n", *icase, *n, *incx, *incy, *icomp, *itrue, id)
Label40:
	;
	return
}

func itest1C(icomp *C.int, itrue *C.int) {
	var nout, icase *int
	var incx, incy, n *int
	var pass *bool
	var id C.int
	//     ********************************* ITEST1 *************************
	//
	//     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
	//     EQUALITY.
	//     C. L. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	if *icomp == *itrue {
		goto Label40
	}
	//
	//                            HERE ICOMP IS NOT EQUAL TO ITRUE.
	//
	if !*pass {
		goto Label20
	}
	//                             PRINT FAIL MESSAGE AND HEADER.
	*pass = false
	writeString(*nout, "                                       FAIL\n")
	writeString(*nout, "\n CASE  N INCX INCY                                COMP                                TRUE     DIFFERENCE\n \n")
Label20:
	;
	id = (*icomp) - (*itrue)
	writeString(*nout, " %4d%3d%5d%5d%36d%36d%12d%12d\n", *icase, *n, *incx, *incy, *icomp, *itrue, id)
Label40:
	;
	return
}

func lce(ri *[]complex64, rj *[]complex64, lr *int) bool {
	var lceReturn bool
	var i int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lceReturn = true
	goto Label30
Label20:
	;
	lceReturn = false
Label30:
	;
	return lceReturn
}

func lceC(ri *[]C.complexfloat, rj *[]C.complexfloat, lr *C.int) bool {
	var lceReturn bool
	var i int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= int(*lr); i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lceReturn = true
	goto Label30
Label20:
	;
	lceReturn = false
Label30:
	;
	return lceReturn
}

func lceres(_type string, uplo byte, m *int, n *int, aa *[][]complex64, as *[][]complex64, lda *int) bool {
	var lceresReturn bool
	var i, ibeg, iend, j int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'HE' or 'HP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "HE" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lceresReturn = true
	goto Label80
Label70:
	;
	lceresReturn = false
Label80:
	;
	return lceresReturn
}

func lceresC(_type string, uplo byte, m *C.int, n *C.int, aa *[][]C.complexfloat, as *[][]C.complexfloat, lda *C.int) bool {
	var lceresReturn bool
	var i, ibeg, iend, j C.int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'HE' or 'HP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "HE" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lceresReturn = true
	goto Label80
Label70:
	;
	lceresReturn = false
Label80:
	;
	return lceresReturn
}

func lde(ri *[]float64, rj *[]float64, lr *int) bool {
	var ldeReturn bool
	var i int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	ldeReturn = true
	goto Label30
Label20:
	;
	ldeReturn = false
Label30:
	;
	return ldeReturn
}

func ldeC(ri *[]C.double, rj *[]C.double, lr *C.int) bool {
	var ldeReturn bool
	var i C.int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	ldeReturn = true
	goto Label30
Label20:
	;
	ldeReturn = false
Label30:
	;
	return ldeReturn
}

func lderes(_type string, uplo byte, m *int, n *int, aa *[][]float64, as *[][]float64, lda *int) bool {
	var lderesReturn bool
	var i, ibeg, iend, j int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'SY' or 'SP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "SY" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lderesReturn = true
	goto Label80
Label70:
	;
	lderesReturn = false
Label80:
	;
	return lderesReturn
}

func lderesC(_type string, uplo byte, m *C.int, n *C.int, aa *[][]C.double, as *[][]C.double, lda *C.int) bool {
	var lderesReturn bool
	var i, ibeg, iend, j C.int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'SY' or 'SP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "SY" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lderesReturn = true
	goto Label80
Label70:
	;
	lderesReturn = false
Label80:
	;
	return lderesReturn
}

func lse(ri *[]float32, rj *[]float32, lr *int) bool {
	var lseReturn bool
	var i int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lseReturn = true
	goto Label30
Label20:
	;
	lseReturn = false
Label30:
	;
	return lseReturn
}

func lseC(ri *[]C.float, rj *[]C.float, lr *C.int) bool {
	var lseReturn bool
	var i C.int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lseReturn = true
	goto Label30
Label20:
	;
	lseReturn = false
Label30:
	;
	return lseReturn
}

func lseres(_type string, uplo byte, m *int, n *int, aa *[][]float32, as *[][]float32, lda *int) bool {
	var lseresReturn bool
	var i, ibeg, iend, j int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'SY' or 'SP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "SY" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lseresReturn = true
	goto Label80
Label70:
	;
	lseresReturn = false
Label80:
	;
	return lseresReturn
}

func lseresC(_type string, uplo byte, m *C.int, n *C.int, aa *[][]C.float, as *[][]C.float, lda *C.int) bool {
	var lseresReturn bool
	var i, ibeg, iend, j C.int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'SY' or 'SP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "SY" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lseresReturn = true
	goto Label80
Label70:
	;
	lseresReturn = false
Label80:
	;
	return lseresReturn
}

func lze(ri *[]complex128, rj *[]complex128, lr *int) bool {
	var lzeReturn bool
	var i int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lzeReturn = true
	goto Label30
Label20:
	;
	lzeReturn = false
Label30:
	;
	return lzeReturn
}

func lzeC(ri *[]C.complexdouble, rj *[]C.complexdouble, lr *C.int) bool {
	var lzeReturn bool
	var i C.int
	//
	//  Tests if two arrays are identical.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	for i = 1; i <= *lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			goto Label20
		}
	}
	lzeReturn = true
	goto Label30
Label20:
	;
	lzeReturn = false
Label30:
	;
	return lzeReturn
}

func lzeres(_type string, uplo byte, m *int, n *int, aa *[][]complex128, as *[][]complex128, lda *int) bool {
	var lzeresReturn bool
	var i, ibeg, iend, j int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'HE' or 'HP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "HE" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lzeresReturn = true
	goto Label80
Label70:
	;
	lzeresReturn = false
Label80:
	;
	return lzeresReturn
}

func lzeresC(_type string, uplo byte, m *C.int, n *C.int, aa *[][]C.complexdouble, as *[][]C.complexdouble, lda *C.int) bool {
	var lzeresReturn bool
	var i, ibeg, iend, j C.int
	var upper bool
	//
	//  Tests if selected elements in two arrays are equal.
	//
	//  _type is 'GE', 'HE' or 'HP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = (*m) + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	} else if _type == "HE" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
			for i = iend + 1; i <= *lda; i++ {
				if (*aa)[i-1][j-1] != (*as)[i-1][j-1] {
					goto Label70
				}
			}
		}
	}

	lzeresReturn = true
	goto Label80
Label70:
	;
	lzeresReturn = false
Label80:
	;
	return lzeresReturn
}

func sbeg(reset *bool) float32 {
	var sbegReturn float32
	var i, ic, mi *int
	//
	//  Generates random numbers uniformly distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	ic = &common.beg.ic
	i = &common.beg.i

	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*i -= 1000 * ((*i) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	sbegReturn = float32((*i)-500) / 1001.0
	return sbegReturn
}

func sbegC(reset *bool) C.float {
	var sbegReturn C.float
	var i, ic, mi *int
	//
	//  Generates random numbers uniformly distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	ic = &common.beg.ic
	i = &common.beg.i

	if *reset {
		//        Initialize local variables.
		*mi = 891
		*i = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of I is bounded between 1 and 999.
	//     If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial I = 4 or 8, the period will be 25.
	//     If initial I = 5, the period will be 10.
	//     IC is used to break up the period by skipping 1 value of I in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*i -= 1000 * ((*i) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	sbegReturn = C.float((*i)-500) / 1001.0
	return sbegReturn
}

func sdiff(sa *float32, sb *float32) float32 {
	//     ********************************* SDIFF **************************
	//     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
	//
	return (*sa) - (*sb)
}

func smakeL2(_type string, uplo byte, diag byte, m *int, n *int, a *[][]float32, nmax *int, aa *[]float32, lda *int, kl *int, ku *int, reset *bool, transl *float32) {
	var zero, one, rogue float32
	var i, i1, i2, i3, ibeg, iend, ioff, j, kk int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'S'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbeg(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array AS in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-1-1][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = max(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = min((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff = ioff + 1
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
				}
			}
		}
	}
	return
}

func smakeL2C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.float, nmax *C.int, aa *[]C.float, lda *C.int, kl *C.int, ku *C.int, reset *bool, transl *C.float) {
	var zero, one, rogue C.float
	var i, i1, i2, i3, ibeg, iend, ioff, j, kk C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type[0] == 'G'
	sym = _type[0] == 'S'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = sbegC(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array AS in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= mincint((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-1-1][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = maxcint(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = mincint((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff = ioff + 1
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
				}
			}
		}
	}
	return
}

func smakeL3(_type string, uplo byte, diag byte, m *int, n *int, a *[][]float32, nmax *int, aa *[]float32, lda *int, reset *bool, transl *float32) {
	var zero, one, rogue float32
	var i, ibeg, iend, j int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type == "GE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = sbeg(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= (*lda); i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	}
	return
}

func smakeL3C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.float, nmax *C.int, aa *[]C.float, lda *C.int, reset *bool, transl *C.float) {
	var zero, one, rogue C.float
	var i, ibeg, iend, j C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0
	rogue = -1.0e10

	gen = _type == "GE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = sbegC(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= (*lda); i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	}
	return
}

func smmch(transa byte, transb byte, m *int, n *int, kk *int, alpha *float32, a *[][]float32, lda *int, b *[][]float32, ldb *int, beta *float32, c *[][]float32, ldc *int, ct *[]float32, g *[]float32, cc *[][]float32, ldcc *int, eps *float32, err *float32, fatal *bool, nout *int, mv bool) {
	var zero, one, erri float32
	var i, j, k int
	var trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and C.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = zero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf32((*a)[i-1][k-1]) * absf32((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += absf32((*a)[k-1][i-1]) * absf32((*b)[k-1][j-1])
				}
			}
		} else if !trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf32((*a)[i-1][k-1]) * absf32((*b)[j-1][k-1])
				}
			}
		} else if trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += absf32((*a)[k-1][i-1]) * absf32((*b)[j-1][k-1])
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = absf32(*alpha)*(*g)[i-1] + absf32(*beta)*absf32((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = zero
		for i = 1; i <= *m; i++ {
			erri = absf32((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != zero {
				erri /= (*g)[i-1]
			}
			*err = maxf32(*err, erri)
			if (*err)*sqrtf32(*eps) >= one {
				goto Label130
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label150
	//
	//     Report fatal error.
	//
Label130:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label150:
	;
	return
}

func smmchC(transa byte, transb byte, m *C.int, n *C.int, kk *C.int, alpha *C.float, a *[][]C.float, lda *C.int, b *[][]C.float, ldb *C.int, beta *C.float, c *[][]C.float, ldc *C.int, ct *[]C.float, g *[]C.float, cc *[][]C.float, ldcc *C.int, eps *C.float, err *C.float, fatal *bool, nout *int, mv bool) {
	var zero, one, erri C.float
	var i, j, k C.int
	var trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = 0.0
	one = 1.0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and C.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = zero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abscfloat((*a)[i-1][k-1]) * abscfloat((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abscfloat((*a)[k-1][i-1]) * abscfloat((*b)[k-1][j-1])
				}
			}
		} else if !trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
					(*g)[i-1] += abscfloat((*a)[i-1][k-1]) * abscfloat((*b)[j-1][k-1])
				}
			}
		} else if trana && tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
					(*g)[i-1] += abscfloat((*a)[k-1][i-1]) * abscfloat((*b)[j-1][k-1])
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abscfloat(*alpha)*(*g)[i-1] + abscfloat(*beta)*abscfloat((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = zero
		for i = 1; i <= *m; i++ {
			erri = abscfloat((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != zero {
				erri = erri / (*g)[i-1]
			}
			*err = maxcfloat(*err, erri)
			if (*err)*sqrtcfloat(*eps) >= one {
				goto Label130
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label150
	//
	//     Report fatal error.
	//
Label130:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label150:
	;
	return
}

func smvch(trans byte, m *int, n *int, alpha *float32, a *[][]float32, nmax *int, x *[]float32, incx *int, beta *float32, y *[]float32, incy *int, yt *[]float32, g *[]float32, yy *[]float32, eps *float32, err *float32, fatal *bool, nout *int, mv bool) {
	var zero, one, erri float32
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl int
	var tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0

	tran = trans == 'T' || trans == 'C'
	if tran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = zero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = absf32(*alpha)*(*g)[iy-1] + absf32((*beta)*(*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = zero
	for i = 1; i <= ml; i++ {
		erri = absf32((*yt)[i-1]-(*yy)[1+(i-1)*absint(*incy)-1]) / (*eps)
		if (*g)[i-1] != zero {
			erri /= (*g)[i-1]
		}
		*err = maxf32(*err, erri)
		if (*err)*sqrtf32(*eps) >= one {
			goto Label50
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label70
	//
	//     Report fatal error.
	//
Label50:
	;
	*fatal = true
	writeString(*nout, " ******* fatal ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*absint(*incy)-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yy)[1+(i-1)*absint(*incy)-1], (*yt)[i-1])
		}
	}

Label70:
	;
	return
}

func smvchC(trans byte, m *C.int, n *C.int, alpha *C.float, a *[][]C.float, nmax *C.int, x *[]C.float, incx *C.int, beta *C.float, y *[]C.float, incy *C.int, yt *[]C.float, g *[]C.float, yy *[]C.float, eps *C.float, err *C.float, fatal *bool, nout *int, mv bool) {
	var zero, one, erri C.float
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl C.int
	var tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = 0.0
	one = 1.0

	tran = trans == 'T' || trans == 'C'
	if tran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = zero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abscfloat((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abscfloat((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abscfloat(*alpha)*(*g)[iy-1] + abscfloat((*beta)*(*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = zero
	for i = 1; i <= ml; i++ {
		erri = abscfloat((*yt)[i-1]-(*yy)[1+(i-1)*abscint(*incy)-1]) / (*eps)
		if (*g)[i-1] != zero {
			erri /= (*g)[i-1]
		}
		*err = maxcfloat(*err, erri)
		if (*err)*sqrtcfloat(*eps) >= one {
			goto Label50
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label70
	//
	//     Report fatal error.
	//
Label50:
	;
	*fatal = true
	writeString(*nout, " ******* fatal ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*abscint(*incy)-1])
		} else {
			writeString(*nout, " %7d%18.6f%18.6f\n", i, (*yy)[1+(i-1)*abscint(*incy)-1], (*yt)[i-1])
		}
	}

Label70:
	;
	return
}

func stest(len int, scomp *[]float32, strue *[]float32, ssize *[]float32, sfac *float32) {
	var nout *int
	var icase *int
	var incx, incy, n *int
	var pass *bool
	var sd float32
	var i int
	//     ********************************* STEST **************************
	//
	//     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
	//     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY sfac, ARE
	//     NEGLIGIBLE.
	//
	//     C. L. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	for i = 1; i <= len; i++ {
		sd = (*scomp)[i-1] - (*strue)[i-1]
		if absf32((*sfac)*sd) <= absf32((*ssize)[i-1])*epsilonf32()+1e-9 {
			goto Label40
		}
		//
		//                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
		//
		if !(*pass) {
			goto Label20
		}
		//                             PRINT FAIL MESSAGE AND HEADER.
		*pass = false
		writeString(*nout, "                                       FAIL\n")
		writeString(*nout, "\n CASE  N incx incy  I                             COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n \n")
	Label20:
		;
		writeString(*nout, " %4d%3d%5d%5d%3d%36.8e%36.8e%12.4e%12.4e\n", *icase, *n, *incx, *incy, i, (*scomp)[i-1], (*strue)[i-1], sd, (*ssize)[i-1])
	Label40:
	}
	return
}

func stestC(len C.int, scomp *[]C.float, strue *[]C.float, ssize *[]C.float, sfac *C.float) {
	var nout, icase *int
	var incx, incy, n *int
	var pass *bool
	var sd C.float
	var i C.int
	//     ********************************* STEST **************************
	//
	//     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
	//     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY sfac, ARE
	//     NEGLIGIBLE.
	//
	//     C. L. LAWSON, JPL, 1974 DEC 10
	//
	nout = &common.combla.nout
	icase = &common.combla.icase
	n = &common.combla.n
	incx = &common.combla.incx
	incy = &common.combla.incy
	pass = &common.combla.pass

	for i = 1; i <= len; i++ {
		sd = (*scomp)[i-1] - (*strue)[i-1]
		if abscfloat((*sfac)*sd) <= abscfloat((*ssize)[i-1])*epsiloncfloat()+1e-9 {
			goto Label40
		}
		//
		//                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
		//
		if !(*pass) {
			goto Label20
		}
		//                             PRINT FAIL MESSAGE AND HEADER.
		*pass = false
		writeString(*nout, "                                       FAIL\n")
		writeString(*nout, "\n CASE  N incx incy  I                             COMP(I)                             TRUE(I)  DIFFERENCE     SIZE(I)\n \n")
	Label20:
		;
		writeString(*nout, " %4d%3d%5d%5d%3d%36.8e%36.8e%12.4e%12.4e\n", *icase, *n, *incx, *incy, i, (*scomp)[i-1], (*strue)[i-1], sd, (*ssize)[i-1])
	Label40:
	}
	return
}

func stest1(scomp1 *float32, strue1 *float32, ssize *[]float32, sfac *float32) {
	scomp := make([]float32, 1)
	strue := make([]float32, 1)
	//     ************************* STEST1 *****************************
	//
	//     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//
	//     C.L. LAWSON, JPL, 1978 DEC 6
	//
	scomp[0] = *scomp1
	strue[0] = *strue1
	stest(1, &scomp, &strue, ssize, sfac)

	return
}

func stest1C(scomp1 *C.float, strue1 *C.float, ssize *[]C.float, sfac *C.float) {
	scomp := make([]C.float, 1)
	strue := make([]C.float, 1)
	//     ************************* STEST1 *****************************
	//
	//     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//
	//     C.L. LAWSON, JPL, 1978 DEC 6
	//
	scomp[0] = *scomp1
	strue[0] = *strue1
	stestC(1, &scomp, &strue, ssize, sfac)

	return
}

func xerblaTest(srname *string, info *int) {
	var infox *int
	var lerr *bool
	var srnamx *string
	//
	//  This is a special version of XERBLA to be used only as part of
	//  the test program for testing error exits from the Level 2 & 3 BLAS
	//  routines.
	//
	//  XERBLA  is an error handler for the Level 3 BLAS routines.
	//
	//  It is called by the Level 2 & 3 BLAS routines if an input parameter is
	//  invalid.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	infox = &common.infoc.infox
	lerr = &common.infoc.lerr
	srnamx = &common.srnamc.srnamx

	*infox = *info
	*srnamx = *srname
	*lerr = true
	return
}

func zbeg(reset *bool) complex128 {
	var zbegReturn complex128
	var i, j, ic, mi, mj *int
	//
	//  Generates complex numbers as pairs of random numbers uniformly
	//  distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	mj = &common.beg.mj
	ic = &common.beg.ic
	i = &common.beg.i
	j = &common.beg.j

	if *(reset) {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i or j is bounded between 1 and 999.
	//     If initial i or j = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i or j = 4 or 8, the period will be 25.
	//     If initial i or j = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i or j
	//     in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*j *= *mj
	*i -= 1000 * ((*i) / 1000)
	*j -= 1000 * ((*j) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	zbegReturn = cmplxc128(float64((*i)-500)/1001.0, float64((*j)-500)/1001.0)
	return zbegReturn
}

func zbegC(reset *bool) C.complexdouble {
	var zbegReturn C.complexdouble
	var i, j, ic, mi, mj *int
	//
	//  Generates complex numbers as pairs of random numbers uniformly
	//  distributed between -0.5 and 0.5.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	mi = &common.beg.mi
	mj = &common.beg.mj
	ic = &common.beg.ic
	i = &common.beg.i
	j = &common.beg.j

	if *(reset) {
		//        Initialize local variables.
		*mi = 891
		*mj = 457
		*i = 7
		*j = 7
		*ic = 0
		*reset = false
	}
	//
	//     The sequence of values of i or j is bounded between 1 and 999.
	//     If initial i or j = 1,2,3,6,7 or 9, the period will be 50.
	//     If initial i or j = 4 or 8, the period will be 25.
	//     If initial i or j = 5, the period will be 10.
	//     ic is used to break up the period by skipping 1 value of i or j
	//     in 6.
	//
	*ic++
Label10:
	;
	*i *= *mi
	*j *= *mj
	*i -= 1000 * ((*i) / 1000)
	*j -= 1000 * ((*j) / 1000)
	if *ic >= 5 {
		*ic = 0
		goto Label10
	}
	zbegReturn = cmplxcdouble(C.double((*i)-500)/1001.0, C.double((*j)-500)/1001.0)
	return zbegReturn
}

func zmakeL2(_type string, uplo byte, diag byte, m *int, n *int, a *[][]complex128, nmax *int, aa *[]complex128, lda *int, kl *int, ku *int, reset *bool, transl *complex128) {
	var zero, one, rogue complex128
	var rzero, rrogue float64
	var i, i1, i2, i3, ibeg, iend, ioff, j, jj, kk int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0e0 + (0.0e0)*1i)
	one = (1.0e0 + (0.0e0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0e0
	rrogue = -1.0e10

	gen = _type[0] == 'g'
	sym = _type[0] == 'H'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbeg(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = conjgc128((*a)[i-1][j-1])
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if sym {
			(*a)[j-1][j-1] = cmplxc128(realc128((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= min((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-0][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc128(realc128((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = max(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = min((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = kk + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc128(realc128((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
					if sym {
						(*aa)[ioff-1] = cmplxc128(realc128((*aa)[ioff-1]), rrogue)
					}
				}
			}
		}
	}
	return
}

func zmakeL2C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.complexdouble, nmax *C.int, aa *[]C.complexdouble, lda *C.int, kl *C.int, ku *C.int, reset *bool, transl *C.complexdouble) {
	var zero, one, rogue C.complexdouble
	var rzero, rrogue C.double
	var i, i1, i2, i3, ibeg, iend, ioff, j, jj, kk C.int
	var gen, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a within the bandwidth
	//  defined by kl and ku.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0e0 + (0.0e0)*1i)
	one = (1.0e0 + (0.0e0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0e0
	rrogue = -1.0e10

	gen = _type[0] == 'g'
	sym = _type[0] == 'H'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= *ku) || (i >= j && i-j <= *kl) {
					(*a)[i-1][j-1] = zbegC(reset) + (*transl)
				} else {
					(*a)[i-1][j-1] = zero
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = conjgcdouble((*a)[i-1][j-1])
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if sym {
			(*a)[j-1][j-1] = cmplxcdouble(realcdouble((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "GB" {
		for j = 1; j <= *n; j++ {
			for i1 = 1; i1 <= (*ku)+1-j; i1++ {
				(*aa)[i1+(j-1)*(*lda)-1] = rogue
			}
			for i2 = i1; i2 <= mincint((*kl)+(*ku)+1, (*ku)+1+(*m)-j); i2++ {
				(*aa)[i2+(j-1)*(*lda)-1] = (*a)[i2+j-(*ku)-0][j-1]
			}
			for i3 = i2; i3 <= *lda; i3++ {
				(*aa)[i3+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcdouble(realcdouble((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HB" || _type == "TB" {
		for j = 1; j <= *n; j++ {
			if upper {
				kk = (*kl) + 1
				ibeg = maxcint(1, (*kl)+2-j)
				if unit {
					iend = *kl
				} else {
					iend = (*kl) + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = mincint((*kl)+1, 1+(*m)-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i+j-kk-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if sym {
				jj = kk + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcdouble(realcdouble((*aa)[jj-1]), rrogue)
			}
		}
	} else if _type == "HP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = *n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
					if sym {
						(*aa)[ioff-1] = cmplxcdouble(realcdouble((*aa)[ioff-1]), rrogue)
					}
				}
			}
		}
	}
	return
}

func zmakeL3(_type string, uplo byte, diag byte, m *int, n *int, a *[][]complex128, nmax *int, aa *[]complex128, lda *int, reset *bool, transl *complex128) {
	var zero, one, rogue complex128
	var rzero, rrogue float64
	var i, ibeg, iend, j, jj int
	var gen, her, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'HE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	//     .. Parameters ..
	zero = (0.0e0 + (0.0e0)*1i)
	one = (1.0e0 + (0.0e0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0e0
	rrogue = -1.0e10

	gen = _type == "GE"
	her = _type == "HE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (her || sym || tri) && uplo == 'U'
	lower = (her || sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = zbeg(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if her {
						(*a)[j-1][i-1] = conjgc128((*a)[i-1][j-1])
					} else if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if her {
			(*a)[j-1][j-1] = cmplxc128(realc128((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if her {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxc128(realc128((*aa)[jj-1]), rrogue)
			}
		}
	}
	return
}

func zmakeL3C(_type string, uplo byte, diag byte, m *C.int, n *C.int, a *[][]C.complexdouble, nmax *C.int, aa *[]C.complexdouble, lda *C.int, reset *bool, transl *C.complexdouble) {
	var zero, one, rogue C.complexdouble
	var rzero, rrogue C.double
	var i, ibeg, iend, j, jj C.int
	var gen, her, lower, sym, tri, unit, upper bool
	//
	//  Generates values for an m by n matrix a.
	//  Stores the values in the array aa in the data structure required
	//  by the routine, with unwanted elements set to rogue value.
	//
	//  _type is 'GE', 'HE', 'SY' or 'TR'.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	//     .. Parameters ..
	zero = (0.0e0 + (0.0e0)*1i)
	one = (1.0e0 + (0.0e0)*1i)
	rogue = (-1.0e10 + (1.0e10)*1i)
	rzero = 0.0e0
	rrogue = -1.0e10

	gen = _type == "GE"
	her = _type == "HE"
	sym = _type == "SY"
	tri = _type == "TR"
	upper = (her || sym || tri) && uplo == 'U'
	lower = (her || sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//     Generate data in array a.
	//
	for j = 1; j <= *n; j++ {
		for i = 1; i <= *m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				(*a)[i-1][j-1] = zbegC(reset) + (*transl)
				if i != j {
					//                 Set some elements to zero
					if *n > 3 && j == (*n)/2 {
						(*a)[i-1][j-1] = zero
					}
					if her {
						(*a)[j-1][i-1] = conjgcdouble((*a)[i-1][j-1])
					} else if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = zero
					}
				}
			}
		}
		if her {
			(*a)[j-1][j-1] = cmplxcdouble(realcdouble((*a)[j-1][j-1]), rzero)
		}
		if tri {
			(*a)[j-1][j-1] = (*a)[j-1][j-1] + one
		}
		if unit {
			(*a)[j-1][j-1] = one
		}
	}
	//
	//     Store elements in array as in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= *n; j++ {
			for i = 1; i <= *m; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = (*m) + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
		}
	} else if _type == "HE" || _type == "SY" || _type == "TR" {
		for j = 1; j <= *n; j++ {
			if upper {
				ibeg = 1
				if unit {
					iend = j - 1
				} else {
					iend = j
				}
			} else {
				if unit {
					ibeg = j + 1
				} else {
					ibeg = j
				}
				iend = *n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = (*a)[i-1][j-1]
			}
			for i = iend + 1; i <= *lda; i++ {
				(*aa)[i+(j-1)*(*lda)-1] = rogue
			}
			if her {
				jj = j + (j-1)*(*lda)
				(*aa)[jj-1] = cmplxcdouble(realcdouble((*aa)[jj-1]), rrogue)
			}
		}
	}
	return
}

func zmmch(transa byte, transb byte, m *int, n *int, kk *int, alpha *complex128, a *[][]complex128, lda *int, b *[][]complex128, ldb *int, beta *complex128, c *[][]complex128, ldc *int, ct *[]complex128, g *[]float64, cc *[][]complex128, LDCC *int, eps *float64, err *float64, fatal *bool, nout *int, mv bool) {
	var zero complex128
	var rzero, rone float64
	var erri float64
	var i, j, k int
	var ctrana, ctranb, trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0e0 + (0.0e0)*1i)
	rzero = 0.0e0
	rone = 1.0e0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	ctrana = transa == 'C'
	ctranb = transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = rzero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			if ctrana {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += conjgc128((*a)[k-1][i-1]) * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
						(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[k-1][j-1])
					}
				}
			}
		} else if !trana && tranb {
			if ctranb {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * conjgc128((*b)[j-1][k-1])
						(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumf64((*a)[i-1][k-1]) * abssumf64((*b)[j-1][k-1])
					}
				}
			}
		} else if trana && tranb {
			if ctrana {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += conjgc128((*a)[k-1][i-1]) * conjgc128((*b)[j-1][k-1])
							(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += conjgc128((*a)[k-1][i-1]) * (*b)[j-1][k-1]
							(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
						}
					}
				}
			} else {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += (*a)[k-1][i-1] * conjgc128((*b)[j-1][k-1])
							(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
							(*g)[i-1] += abssumf64((*a)[k-1][i-1]) * abssumf64((*b)[j-1][k-1])
						}
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumf64(*alpha)*(*g)[i-1] + abssumf64(*beta)*abssumf64((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = rzero
		for i = 1; i <= *m; i++ {
			erri = abssumf64((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != rzero {
				erri /= (*g)[i-1]
			}
			*err = maxf64(*err, erri)
			if (*err)*sqrtf64(*eps) >= rone {
				goto Label230
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label250
	//
	//     Report fatal error.
	//
Label230:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label250:
	;
	return
}

func zmmchC(transa byte, transb byte, m *C.int, n *C.int, kk *C.int, alpha *C.complexdouble, a *[][]C.complexdouble, lda *C.int, b *[][]C.complexdouble, ldb *C.int, beta *C.complexdouble, c *[][]C.complexdouble, ldc *C.int, ct *[]C.complexdouble, g *[]C.double, cc *[][]C.complexdouble, LDCC *C.int, eps *C.double, err *C.double, fatal *bool, nout *int, mv bool) {
	var zero C.complexdouble
	var rzero, rone C.double
	var erri C.double
	var i, j, k C.int
	var ctrana, ctranb, trana, tranb bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 3 Blas.
	//
	//  -- Written on 8-February-1989.
	//     Jack Dongarra, Argonne National Laboratory.
	//     Iain Duff, AERE Harwell.
	//     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//     Sven Hammarling, Numerical Algorithms Group Ltd.
	//
	zero = (0.0e0 + (0.0e0)*1i)
	rzero = 0.0e0
	rone = 1.0e0

	trana = transa == 'T' || transa == 'C'
	tranb = transb == 'T' || transb == 'C'
	ctrana = transa == 'C'
	ctranb = transb == 'C'
	//
	//     Compute expected result, one column at a time, in ct using data
	//     in a, b and c.
	//     Compute gauges in g.
	//
	for j = 1; j <= *n; j++ {

		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = zero
			(*g)[i-1] = rzero
		}
		if !trana && !tranb {
			for k = 1; k <= *kk; k++ {
				for i = 1; i <= *m; i++ {
					(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[k-1][j-1]
					(*g)[i-1] += abssumcdouble((*a)[i-1][k-1]) * abssumcdouble((*b)[k-1][j-1])
				}
			}
		} else if trana && !tranb {
			if ctrana {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += conjgcdouble((*a)[k-1][i-1]) * (*b)[k-1][j-1]
						(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[k-1][j-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[k-1][j-1]
						(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[k-1][j-1])
					}
				}
			}
		} else if !trana && tranb {
			if ctranb {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * conjgcdouble((*b)[j-1][k-1])
						(*g)[i-1] += abssumcdouble((*a)[i-1][k-1]) * abssumcdouble((*b)[j-1][k-1])
					}
				}
			} else {
				for k = 1; k <= *kk; k++ {
					for i = 1; i <= *m; i++ {
						(*ct)[i-1] += (*a)[i-1][k-1] * (*b)[j-1][k-1]
						(*g)[i-1] += abssumcdouble((*a)[i-1][k-1]) * abssumcdouble((*b)[j-1][k-1])
					}
				}
			}
		} else if trana && tranb {
			if ctrana {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += conjgcdouble((*a)[k-1][i-1]) * conjgcdouble((*b)[j-1][k-1])
							(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += conjgcdouble((*a)[k-1][i-1]) * (*b)[j-1][k-1]
							(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[j-1][k-1])
						}
					}
				}
			} else {
				if ctranb {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += (*a)[k-1][i-1] * conjgcdouble((*b)[j-1][k-1])
							(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[j-1][k-1])
						}
					}
				} else {
					for k = 1; k <= *kk; k++ {
						for i = 1; i <= *m; i++ {
							(*ct)[i-1] += (*a)[k-1][i-1] * (*b)[j-1][k-1]
							(*g)[i-1] += abssumcdouble((*a)[k-1][i-1]) * abssumcdouble((*b)[j-1][k-1])
						}
					}
				}
			}
		}
		for i = 1; i <= *m; i++ {
			(*ct)[i-1] = (*alpha)*(*ct)[i-1] + (*beta)*(*c)[i-1][j-1]
			(*g)[i-1] = abssumcdouble(*alpha)*(*g)[i-1] + abssumcdouble(*beta)*abssumcdouble((*c)[i-1][j-1])
		}
		//
		//        Compute the error ratio for this result.
		//
		*err = rzero
		for i = 1; i <= *m; i++ {
			erri = abssumcdouble((*ct)[i-1]-(*cc)[i-1][j-1]) / (*eps)
			if (*g)[i-1] != rzero {
				erri /= (*g)[i-1]
			}
			*err = maxcdouble(*err, erri)
			if (*err)*sqrtcdouble(*eps) >= rone {
				goto Label230
			}
		}

	}
	//
	//     If the loop completes, all results are at least half accurate.
	goto Label250
	//
	//     Report fatal error.
	//
Label230:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= *m; i++ {
		if mv {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*ct)[i-1], (*cc)[i-1][j-1])
		} else {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*cc)[i-1][j-1], (*ct)[i-1])
		}
	}
	if *n > 1 {
		writeString(*nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n", j)
	}

Label250:
	;
	return
}

func zmvch(trans byte, m *int, n *int, alpha *complex128, a *[][]complex128, nmax *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int, yt *[]complex128, g *[]float64, yy *[]complex128, eps *float64, err *float64, fatal *bool, nout *int, mv bool) {
	var zero complex128
	var rzero, rone float64
	var erri float64
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl int
	var ctran, tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0e0 + (0.0e0)*1i)
	rzero = 0.0e0
	rone = 1.0e0

	tran = trans == 'T'
	ctran = trans == 'C'
	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = rzero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[j-1][i-1]) * abssumf64((*x)[jx-1])
				jx = jx + incxl
				//Label10:
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += conjgc128((*a)[j-1][i-1]) * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[j-1][i-1]) * abssumf64((*x)[jx-1])
				jx = jx + incxl
				//Label20:
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumf64((*a)[i-1][j-1]) * abssumf64((*x)[jx-1])
				jx = jx + incxl
				//Label30:
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumf64(*alpha)*(*g)[iy-1] + abssumf64(*beta)*abssumf64((*y)[iy-1])
		iy = iy + incyl
		//Label40:
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = rzero
	for i = 1; i <= ml; i++ {
		erri = absc128((*yt)[i-1]-(*yy)[1+(i-1)*absint(*incy)-1]) / (*eps)
		if (*g)[i-1] != rzero {
			erri = erri / (*g)[i-1]
		}
		*err = maxf64(*err, erri)
		if (*err)*sqrtf64(*eps) >= rone {
			goto Label60
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label80
	//
	//     Report fatal error.
	//
Label60:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*absint(*incy)-1])
		} else {
			writeString(*nout, " %7d  %15.6f%15.6f\n", i, (*yy)[1+(i-1)*absint(*incy)-1], (*yt)[i-1])
		}
	}

Label80:
	;
	return
}

func zmvchC(trans byte, m *C.int, n *C.int, alpha *C.complexdouble, a *[][]C.complexdouble, nmax *C.int, x *[]C.complexdouble, incx *C.int, beta *C.complexdouble, y *[]C.complexdouble, incy *C.int, yt *[]C.complexdouble, g *[]C.double, yy *[]C.complexdouble, eps *C.double, err *C.double, fatal *bool, nout *int, mv bool) {
	var zero C.complexdouble
	var rzero, rone C.double
	var erri C.double
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl C.int
	var ctran, tran bool
	//
	//  Checks the results of the computational tests.
	//
	//  Auxiliary routine for test program for Level 2 Blas.
	//
	//  -- Written on 10-August-1987.
	//     Richard Hanson, Sandia National Labs.
	//     Jeremy Du Croz, NAG Central Office.
	//
	zero = (0.0 + (0.0)*1i)
	rzero = 0.0
	rone = 1.0

	tran = trans == 'T'
	ctran = trans == 'C'
	if tran || ctran {
		ml = *n
		nl = *m
	} else {
		ml = *m
		nl = *n
	}
	if *incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if *incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//     Compute expected result in yt using data in a, x and y.
	//     Compute gauges in g.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = zero
		(*g)[iy-1] = rzero
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumcdouble((*a)[j-1][i-1]) * abssumcdouble((*x)[jx-1])
				jx += incxl
			}
		} else if ctran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += conjgcdouble((*a)[j-1][i-1]) * (*x)[jx-1]
				(*g)[iy-1] += abssumcdouble((*a)[j-1][i-1]) * abssumcdouble((*x)[jx-1])
				jx += incxl
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += abssumcdouble((*a)[i-1][j-1]) * abssumcdouble((*x)[jx-1])
				jx += incxl
			}
		}
		(*yt)[iy-1] = (*alpha)*(*yt)[iy-1] + (*beta)*(*y)[iy-1]
		(*g)[iy-1] = abssumcdouble(*alpha)*(*g)[iy-1] + abssumcdouble(*beta)*abssumcdouble((*y)[iy-1])
		iy += incyl
	}
	//
	//     Compute the error ratio for this result.
	//
	*err = rzero
	for i = 1; i <= ml; i++ {
		erri = absccdouble((*yt)[i-1]-(*yy)[1+(i-1)*abscint(*incy)-1]) / (*eps)
		if (*g)[i-1] != rzero {
			erri = erri / (*g)[i-1]
		}
		*err = maxcdouble(*err, erri)
		if (*err)*sqrtcdouble(*eps) >= rone {
			goto Label60
		}
	}
	//     If the loop completes, all results are at least half accurate.
	goto Label80
	//
	//     Report fatal error.
	//
Label60:
	;
	*fatal = true
	writeString(*nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(*nout, " %7d  (%15.6f)(%15.6f)\n", i, (*yt)[i-1], (*yy)[1+(i-1)*abscint(*incy)-1])
		} else {
			writeString(*nout, " %7d  (%15.6f)(%15.6f)\n", i, (*yy)[1+(i-1)*abscint(*incy)-1], (*yt)[i-1])
		}
	}

Label80:
	;
	return
}

func ztest(len int, ccomp *[]complex128, ctrue *[]complex128, csize *[]complex128, sfac *float64) {
	var i int
	scomp := make([]float64, 20)
	ssize := make([]float64, 20)
	strue := make([]float64, 20)
	//     **************************** CTEST *****************************
	//
	//     c.l. LAWSON, JPL, 1978 DEC 6
	//
	for i = 1; i <= len; i++ {
		scomp[2*i-0] = real((*ccomp)[i-1])
		scomp[2*i-1] = imag((*ccomp)[i-1])
		strue[2*i-0] = real((*ctrue)[i-1])
		strue[2*i-1] = imag((*ctrue)[i-1])
		ssize[2*i-0] = real((*csize)[i-1])
		ssize[2*i-1] = imag((*csize)[i-1])
	}

	dtest(2*len, &scomp, &strue, &ssize, sfac)
	return
}

func ztestC(len C.int, ccomp *[]C.complexdouble, ctrue *[]C.complexdouble, csize *[]C.complexdouble, sfac *C.double) {
	var i C.int
	scomp := make([]C.double, 20)
	ssize := make([]C.double, 20)
	strue := make([]C.double, 20)
	//     **************************** CTEST *****************************
	//
	//     c.l. LAWSON, JPL, 1978 DEC 6
	//
	for i = 1; i <= len; i++ {
		scomp[2*i-0] = realcdouble((*ccomp)[i-1])
		scomp[2*i-1] = imagcdouble((*ccomp)[i-1])
		strue[2*i-0] = realcdouble((*ctrue)[i-1])
		strue[2*i-1] = imagcdouble((*ctrue)[i-1])
		ssize[2*i-0] = realcdouble((*csize)[i-1])
		ssize[2*i-1] = imagcdouble((*csize)[i-1])
	}

	dtestC(2*len, &scomp, &strue, &ssize, sfac)
	return
}
