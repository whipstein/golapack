package goblas

import 

// \brief \b Dblat1
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       PROGRAM Dblat1
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Test program for the DOUBLE PRECISION Level 1 BLAS.
//
//    Based upon the original BLAS test routine together with:
//    F06EAF Example Program Text
// \endverbatim
//
//  Authors:
//  ========
//
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG Ltd.
//
// \date April 2012
//
// \ingroup double_blas_testing
//
//  =====================================================================
func main() {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sfac := new(float64)
	ic := new(int)
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*
	//*  -- Reference BLAS test routine (version 3.8.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     April 2012
	//*
	//*  =====================================================================
	//*
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. External Subroutines ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*sfac) = 9.765625e-4
	//*     .. Executable Statements ..
	WRITE((*nout), *func() *[]byte{y := []byte(" real blas test program results\n \n"); return &y}())
	for (*ic) = 1; (*ic) <= 13; (*ic)++ {
		(*icase) = (*ic)
		Header()
		//*
		//*        .. Initialize  pass,  incx,  and incy for a new case. ..
		//*        .. the value 9999 for incx or incy will appear in the ..
		//*        .. detailed  output, if any, for cases  that do not involve ..
		//*        .. these parameters ..
		//*
		(*pass) = true
		(*incx) = 9999
		(*incy) = 9999
		if (*icase) == 3 || (*icase) == 11 {
			CHECK0(sfac)
		} else if (*icase) == 7 || (*icase) == 8 || (*icase) == 9 || (*icase) == 10 {
			Check1(sfac)
		} else if (*icase) == 1 || (*icase) == 2 || (*icase) == 5 || (*icase) == 6 || (*icase) == 12 || (*icase) == 13 {
			Check2(sfac)
		} else if (*icase) == 4 {
			CHECK3(sfac)
		}
		//*        -- Print
		if *pass {
			WRITE((*nout), *func() *[]byte{y := []byte("                                    ----- pass -----\n"); return &y}())
		}
	//Label20:
	}
	panic("")
	//*
}

func Header() {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	l := func() *[][]byte {
		arr := make([][]byte, 6)
		for u := 0; u < 6; u++ {
			arr[u] = make([]byte, 13)
		}
		return &arr
	}()
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalars in Common ..
	//*     .. Local Arrays ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*l)[0][0], (*l)[1][0], (*l)[2][0], (*l)[3][0], (*l)[4][0], (*l)[5][0] = ' ', 'd', 'd', 'o', 't', ' '
	(*l)[0][1], (*l)[1][1], (*l)[2][1], (*l)[3][1], (*l)[4][1], (*l)[5][1] = 'd', 'a', 'x', 'p', 'y', ' '
	(*l)[0][2], (*l)[1][2], (*l)[2][2], (*l)[3][2], (*l)[4][2], (*l)[5][2] = 'd', 'r', 'o', 't', 'g', ' '
	(*l)[0][3], (*l)[1][3], (*l)[2][3], (*l)[3][3], (*l)[4][3], (*l)[5][3] = ' ', 'd', 'r', 'o', 't', ' '
	(*l)[0][4], (*l)[1][4], (*l)[2][4], (*l)[3][4], (*l)[4][4], (*l)[5][4] = 'd', 'c', 'o', 'p', 'y', ' '
	(*l)[0][5], (*l)[1][5], (*l)[2][5], (*l)[3][5], (*l)[4][5], (*l)[5][5] = 'd', 's', 'w', 'a', 'p', ' '
	(*l)[0][6], (*l)[1][6], (*l)[2][6], (*l)[3][6], (*l)[4][6], (*l)[5][6] = 'd', 'n', 'r', 'm', '2', ' '
	(*l)[0][7], (*l)[1][7], (*l)[2][7], (*l)[3][7], (*l)[4][7], (*l)[5][7] = 'd', 'a', 's', 'u', 'm', ' '
	(*l)[0][8], (*l)[1][8], (*l)[2][8], (*l)[3][8], (*l)[4][8], (*l)[5][8] = 'd', 's', 'c', 'a', 'l', ' '
	(*l)[0][9], (*l)[1][9], (*l)[2][9], (*l)[3][9], (*l)[4][9], (*l)[5][9] = 'i', 'd', 'a', 'm', 'a', 'x'
	(*l)[0][10], (*l)[1][10], (*l)[2][10], (*l)[3][10], (*l)[4][10], (*l)[5][10] = 'd', 'r', 'o', 't', 'm', 'g'
	(*l)[0][11], (*l)[1][11], (*l)[2][11], (*l)[3][11], (*l)[4][11], (*l)[5][11] = 'd', 'r', 'o', 't', 'm', ' '
	(*l)[0][12], (*l)[1][12], (*l)[2][12], (*l)[3][12], (*l)[4][12], (*l)[5][12] = 'd', 's', 'd', 'o', 't', ' '
	//*     .. Executable Statements ..
	WRITE((*nout), *func() *[]byte{y := []byte("\n test of subprogram number%3d            %6s\n"); return &y}(), (*icase), (*l)[(*icase)-1])
	return
	//*
}

func Check0(sfac *float64) {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sa := new(float64)
	sb := new(float64)
	sc := new(float64)
	ss := new(float64)
	d12 := new(float64)
	i := new(int)
	k := new(int)
	da1 := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	datrue := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	db1 := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	dbtrue := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	dc1 := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	ds1 := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	dab := func() *[][]float64 {
		arr := make([][]float64, 4)
		for u := 0; u < 4; u++ {
			arr[u] = make([]float64, 9)
		}
		return &arr
	}()
	dtemp := func() *[]float64 {
		arr := make([]float64, 9)
		return &arr
	}()
	dtrue := func() *[][]float64 {
		arr := make([][]float64, 9)
		for u := 0; u < 9; u++ {
			arr[u] = make([]float64, 9)
		}
		return &arr
	}()
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Subroutines ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*da1)[0], (*da1)[1], (*da1)[2], (*da1)[3], (*da1)[4], (*da1)[5], (*da1)[6], (*da1)[7] = 0.3, 0.4, -0.3, -0.4, -0.3, 0.0, 0.0, 1.0
	(*db1)[0], (*db1)[1], (*db1)[2], (*db1)[3], (*db1)[4], (*db1)[5], (*db1)[6], (*db1)[7] = 0.4, 0.3, 0.4, 0.3, -0.4, 0.0, 1.0, 0.0
	(*dc1)[0], (*dc1)[1], (*dc1)[2], (*dc1)[3], (*dc1)[4], (*dc1)[5], (*dc1)[6], (*dc1)[7] = 0.6, 0.8, -0.6, 0.8, 0.6, 1.0, 0.0, 1.0
	(*ds1)[0], (*ds1)[1], (*ds1)[2], (*ds1)[3], (*ds1)[4], (*ds1)[5], (*ds1)[6], (*ds1)[7] = 0.8, 0.6, 0.8, -0.6, 0.8, 0.0, 1.0, 0.0
	(*datrue)[0], (*datrue)[1], (*datrue)[2], (*datrue)[3], (*datrue)[4], (*datrue)[5], (*datrue)[6], (*datrue)[7] = 0.5, 0.5, 0.5, -0.5, -0.5, 0.0, 1.0, 1.0
	(*dbtrue)[0], (*dbtrue)[1], (*dbtrue)[2], (*dbtrue)[3], (*dbtrue)[4], (*dbtrue)[5], (*dbtrue)[6], (*dbtrue)[7] = 0.0, 0.6, 0.0, -0.6, 0.0, 0.0, 1.0, 0.0
	//*     INPUT FOR MODIFIED GIVENS
	(*dab)[0][0], (*dab)[1][0], (*dab)[2][0], (*dab)[3][0], (*dab)[0][1], (*dab)[1][1], (*dab)[2][1], (*dab)[3][1], (*dab)[0][2], (*dab)[1][2], (*dab)[2][2], (*dab)[3][2], (*dab)[0][3], (*dab)[1][3], (*dab)[2][3], (*dab)[3][3], (*dab)[0][4], (*dab)[1][4], (*dab)[2][4], (*dab)[3][4], (*dab)[0][5], (*dab)[1][5], (*dab)[2][5], (*dab)[3][5], (*dab)[0][6], (*dab)[1][6], (*dab)[2][6], (*dab)[3][6], (*dab)[0][7], (*dab)[1][7], (*dab)[2][7], (*dab)[3][7], (*dab)[0][8], (*dab)[1][8], (*dab)[2][8], (*dab)[3][8] = .1, .3, 1.2, .2, .7, .2, .6, 4.2, 0., 0., 0., 0., 4., -1., 2., 4., 6.e-10, 2.e-2, 1.e5, 10., 4.e10, 2.e-2, 1.e-5, 10., 2.e-10, 4.e-2, 1.e5, 10., 2.e10, 4.e-2, 1.e-5, 10., 4., -2., 8., 4.
	//*    TRUE RESULTS FOR MODIFIED GIVENS
	(*dtrue)[0][0], (*dtrue)[1][0], (*dtrue)[2][0], (*dtrue)[3][0], (*dtrue)[4][0], (*dtrue)[5][0], (*dtrue)[6][0], (*dtrue)[7][0], (*dtrue)[8][0], (*dtrue)[0][1], (*dtrue)[1][1], (*dtrue)[2][1], (*dtrue)[3][1], (*dtrue)[4][1], (*dtrue)[5][1], (*dtrue)[6][1], (*dtrue)[7][1], (*dtrue)[8][1], (*dtrue)[0][2], (*dtrue)[1][2], (*dtrue)[2][2], (*dtrue)[3][2], (*dtrue)[4][2], (*dtrue)[5][2], (*dtrue)[6][2], (*dtrue)[7][2], (*dtrue)[8][2], (*dtrue)[0][3], (*dtrue)[1][3], (*dtrue)[2][3], (*dtrue)[3][3], (*dtrue)[4][3], (*dtrue)[5][3], (*dtrue)[6][3], (*dtrue)[7][3], (*dtrue)[8][3], (*dtrue)[0][4], (*dtrue)[1][4], (*dtrue)[2][4], (*dtrue)[3][4], (*dtrue)[4][4], (*dtrue)[5][4], (*dtrue)[6][4], (*dtrue)[7][4], (*dtrue)[8][4], (*dtrue)[0][5], (*dtrue)[1][5], (*dtrue)[2][5], (*dtrue)[3][5], (*dtrue)[4][5], (*dtrue)[5][5], (*dtrue)[6][5], (*dtrue)[7][5], (*dtrue)[8][5], (*dtrue)[0][6], (*dtrue)[1][6], (*dtrue)[2][6], (*dtrue)[3][6], (*dtrue)[4][6], (*dtrue)[5][6], (*dtrue)[6][6], (*dtrue)[7][6], (*dtrue)[8][6], (*dtrue)[0][7], (*dtrue)[1][7], (*dtrue)[2][7], (*dtrue)[3][7], (*dtrue)[4][7], (*dtrue)[5][7], (*dtrue)[6][7], (*dtrue)[7][7], (*dtrue)[8][7], (*dtrue)[0][8], (*dtrue)[1][8], (*dtrue)[2][8], (*dtrue)[3][8], (*dtrue)[4][8], (*dtrue)[5][8], (*dtrue)[6][8], (*dtrue)[7][8], (*dtrue)[8][8] = 0., 0., 1.3, .2, 0., 0., 0., .5, 0., 0., 0., 4.5, 4.2, 1., .5, 0., 0., 0., 0., 0., 0., 0., -2., 0., 0., 0., 0., 0., 0., 0., 4., -1., 0., 0., 0., 0., 0., 15.e-3, 0., 10., -1., 0., -1.e-4, 0., 1., 0., 0., 6144.e-5, 10., -1., 4096., -1.e6, 0., 1., 0., 0., 15., 10., -1., 5.e-5, 0., 1., 0., 0., 0., 15., 10., -1.(*d0), 5.e5, -4096., 1., 4096.e-6, 0., 0., 7., 4., 0., 0., -.5, -.25, 0.
	//*                   4096 = 2 ** 12
	(*d12) = 4096.
	(*dtrue)[0][0] = 12. / 130.
	(*dtrue)[1][0] = 36. / 130.
	(*dtrue)[6][0] = -1. / 6.
	(*dtrue)[0][1] = 14. / 75.
	(*dtrue)[1][1] = 49. / 75.
	(*dtrue)[8][1] = 1. / 7.
	(*dtrue)[0][4] = 45.e-11 * ((*d12) * (*d12))
	(*dtrue)[2][4] = 4.e5 / (3. * (*d12))
	(*dtrue)[5][4] = 1. / (*d12)
	(*dtrue)[7][4] = 1.e4 / (3. * (*d12))
	(*dtrue)[0][5] = 4.e10 / (1.5 * (*d12) * (*d12))
	(*dtrue)[1][5] = 2.e-2 / 1.5
	(*dtrue)[7][5] = 5.e-7 * (*d12)
	(*dtrue)[0][6] = 4. / 150.
	(*dtrue)[1][6] = (2.e-10 / 1.5) * ((*d12) * (*d12))
	(*dtrue)[6][6] = -(*dtrue)[5][4]
	(*dtrue)[8][6] = 1.e4 / (*d12)
	(*dtrue)[0][7] = (*dtrue)[0][6]
	(*dtrue)[1][7] = 2.e10 / (1.5 * (*d12) * (*d12))
	(*dtrue)[0][8] = 32. / 7.
	(*dtrue)[1][8] = -16. / 7.
	//*     .. Executable Statements ..
	//*
	//*     Compute true values which cannot be prestored
	//*     in decimal notation
	//*
	(*dbtrue)[0] = 1.0 / 0.6
	(*dbtrue)[2] = -1.0 / 0.6
	(*dbtrue)[4] = 1.0 / 0.6
	//*
	for (*k) = 1; (*k) <= 8; (*k)++ {
		//*        .. Set N=K for identification in output if any ..
		(*n) = (*k)
		if (*icase) == 3 {
			//*           .. Drotg ..
			if (*k) > 8 {
				goto Label40
			}
			(*sa) = (*da1)[(*k)-1]
			(*sb) = (*db1)[(*k)-1]
			Drotg(sa, sb, sc, ss)
			Stest1(sa, &((*datrue)[(*k)-1]), &((*datrue)[(*k)-1]), sfac)
			Stest1(sb, &((*dbtrue)[(*k)-1]), &((*dbtrue)[(*k)-1]), sfac)
			Stest1(sc, &((*dc1)[(*k)-1]), &((*dc1)[(*k)-1]), sfac)
			Stest1(ss, &((*ds1)[(*k)-1]), &((*ds1)[(*k)-1]), sfac)
		} else if (*icase) == 11 {
			//*           .. Drotmg ..
			for (*i) = 1; (*i) <= 4; (*i)++ {
				(*dtemp)[(*i)-1] = (*dab)[(*i)-1][(*k)-1]
				(*dtemp)[(*i)+3] = 0.0
			}
			(*dtemp)[8] = 0.0
			Drotmg(&((*dtemp)[0]), &((*dtemp)[1]), &((*dtemp)[2]), &((*dtemp)[3]), &((*dtemp)[4]))
			Stest(func() *int{y := 9; return &y}(), dtemp, &((*dtrue)[0][(*k)-1]), &((*dtrue)[0][(*k)-1]), sfac)
		} else {
			WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}(), *func() *[]byte{y := []byte(" shouldn't be here in check0"); return &y}())
			panic("")
		}
	//Label20:
	}
Label40:
	;
	return
}

func Check1(sfac *float64) {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	i := new(int)
	len := new(int)
	np1 := new(int)
	dtrue1 := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	dtrue3 := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	dtrue5 := func() *[][][]float64 {
		arr := make([][][]float64, 8)
		for u := 0; u < 8; u++ {
			arr[u] = make([][]float64, 5)
			for w := 0; w < 5; w++ {
				arr[u][w] = make([]float64, 2)
			}
		}
		return &arr
	}()
	dv := func() *[][][]float64 {
		arr := make([][][]float64, 8)
		for u := 0; u < 8; u++ {
			arr[u] = make([][]float64, 5)
			for w := 0; w < 5; w++ {
				arr[u][w] = make([]float64, 2)
			}
		}
		return &arr
	}()
	sa := func() *[]float64 {
		arr := make([]float64, 10)
		return &arr
	}()
	stemp := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	strue := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	sx := func() *[]float64 {
		arr := make([]float64, 8)
		return &arr
	}()
	itrue2 := func() *[]int {
		arr := make([]int, 5)
		return &arr
	}()
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*sa)[0], (*sa)[1], (*sa)[2], (*sa)[3], (*sa)[4], (*sa)[5], (*sa)[6], (*sa)[7], (*sa)[8], (*sa)[9] = 0.3, -1.0, 0.0, 1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3
	(*dv)[0][0][0], (*dv)[1][0][0], (*dv)[2][0][0], (*dv)[3][0][0], (*dv)[4][0][0], (*dv)[5][0][0], (*dv)[6][0][0], (*dv)[7][0][0], (*dv)[0][1][0], (*dv)[1][1][0], (*dv)[2][1][0], (*dv)[3][1][0], (*dv)[4][1][0], (*dv)[5][1][0], (*dv)[6][1][0], (*dv)[7][1][0], (*dv)[0][2][0], (*dv)[1][2][0], (*dv)[2][2][0], (*dv)[3][2][0], (*dv)[4][2][0], (*dv)[5][2][0], (*dv)[6][2][0], (*dv)[7][2][0], (*dv)[0][3][0], (*dv)[1][3][0], (*dv)[2][3][0], (*dv)[3][3][0], (*dv)[4][3][0], (*dv)[5][3][0], (*dv)[6][3][0], (*dv)[7][3][0], (*dv)[0][4][0], (*dv)[1][4][0], (*dv)[2][4][0], (*dv)[3][4][0], (*dv)[4][4][0], (*dv)[5][4][0], (*dv)[6][4][0], (*dv)[7][4][0], (*dv)[0][0][1], (*dv)[1][0][1], (*dv)[2][0][1], (*dv)[3][0][1], (*dv)[4][0][1], (*dv)[5][0][1], (*dv)[6][0][1], (*dv)[7][0][1], (*dv)[0][1][1], (*dv)[1][1][1], (*dv)[2][1][1], (*dv)[3][1][1], (*dv)[4][1][1], (*dv)[5][1][1], (*dv)[6][1][1], (*dv)[7][1][1], (*dv)[0][2][1], (*dv)[1][2][1], (*dv)[2][2][1], (*dv)[3][2][1], (*dv)[4][2][1], (*dv)[5][2][1], (*dv)[6][2][1], (*dv)[7][2][1], (*dv)[0][3][1], (*dv)[1][3][1], (*dv)[2][3][1], (*dv)[3][3][1], (*dv)[4][3][1], (*dv)[5][3][1], (*dv)[6][3][1], (*dv)[7][3][1], (*dv)[0][4][1], (*dv)[1][4][1], (*dv)[2][4][1], (*dv)[3][4][1], (*dv)[4][4][1], (*dv)[5][4][1], (*dv)[6][4][1], (*dv)[7][4][1] = 0.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.3, -0.4, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.2, -0.6, 0.3, 5.0, 5.0, 5.0, 5.0, 5.0, 0.1, -0.3, 0.5, -0.1, 6.0, 6.0, 6.0, 6.0, 0.1, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.3, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.3, 2.0, -0.4, 2.0, 2.0, 2.0, 2.0, 2.0, 0.2, 3.0, -0.6, 5.0, 0.3, 2.0, 2.0, 2.0, 0.1, 4.0, -0.3, 6.0, -0.5, 7.0, -0.1, 3.0
	(*dtrue1)[0], (*dtrue1)[1], (*dtrue1)[2], (*dtrue1)[3], (*dtrue1)[4] = 0.0, 0.3, 0.5, 0.7, 0.6
	(*dtrue3)[0], (*dtrue3)[1], (*dtrue3)[2], (*dtrue3)[3], (*dtrue3)[4] = 0.0, 0.3, 0.7, 1.1, 1.0
	(*dtrue5)[0][0][0], (*dtrue5)[1][0][0], (*dtrue5)[2][0][0], (*dtrue5)[3][0][0], (*dtrue5)[4][0][0], (*dtrue5)[5][0][0], (*dtrue5)[6][0][0], (*dtrue5)[7][0][0], (*dtrue5)[0][1][0], (*dtrue5)[1][1][0], (*dtrue5)[2][1][0], (*dtrue5)[3][1][0], (*dtrue5)[4][1][0], (*dtrue5)[5][1][0], (*dtrue5)[6][1][0], (*dtrue5)[7][1][0], (*dtrue5)[0][2][0], (*dtrue5)[1][2][0], (*dtrue5)[2][2][0], (*dtrue5)[3][2][0], (*dtrue5)[4][2][0], (*dtrue5)[5][2][0], (*dtrue5)[6][2][0], (*dtrue5)[7][2][0], (*dtrue5)[0][3][0], (*dtrue5)[1][3][0], (*dtrue5)[2][3][0], (*dtrue5)[3][3][0], (*dtrue5)[4][3][0], (*dtrue5)[5][3][0], (*dtrue5)[6][3][0], (*dtrue5)[7][3][0], (*dtrue5)[0][4][0], (*dtrue5)[1][4][0], (*dtrue5)[2][4][0], (*dtrue5)[3][4][0], (*dtrue5)[4][4][0], (*dtrue5)[5][4][0], (*dtrue5)[6][4][0], (*dtrue5)[7][4][0], (*dtrue5)[0][0][1], (*dtrue5)[1][0][1], (*dtrue5)[2][0][1], (*dtrue5)[3][0][1], (*dtrue5)[4][0][1], (*dtrue5)[5][0][1], (*dtrue5)[6][0][1], (*dtrue5)[7][0][1], (*dtrue5)[0][1][1], (*dtrue5)[1][1][1], (*dtrue5)[2][1][1], (*dtrue5)[3][1][1], (*dtrue5)[4][1][1], (*dtrue5)[5][1][1], (*dtrue5)[6][1][1], (*dtrue5)[7][1][1], (*dtrue5)[0][2][1], (*dtrue5)[1][2][1], (*dtrue5)[2][2][1], (*dtrue5)[3][2][1], (*dtrue5)[4][2][1], (*dtrue5)[5][2][1], (*dtrue5)[6][2][1], (*dtrue5)[7][2][1], (*dtrue5)[0][3][1], (*dtrue5)[1][3][1], (*dtrue5)[2][3][1], (*dtrue5)[3][3][1], (*dtrue5)[4][3][1], (*dtrue5)[5][3][1], (*dtrue5)[6][3][1], (*dtrue5)[7][3][1], (*dtrue5)[0][4][1], (*dtrue5)[1][4][1], (*dtrue5)[2][4][1], (*dtrue5)[3][4][1], (*dtrue5)[4][4][1], (*dtrue5)[5][4][1], (*dtrue5)[6][4][1], (*dtrue5)[7][4][1] = 0.10, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, -0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.20, -0.60, 0.30, 5.0, 5.0, 5.0, 5.0, 5.0, 0.03, -0.09, 0.15, -0.03, 6.0, 6.0, 6.0, 6.0, 0.10, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.09, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 0.09, 2.0, -0.12, 2.0, 2.0, 2.0, 2.0, 2.0, 0.06, 3.0, -0.18, 5.0, 0.09, 2.0, 2.0, 2.0, 0.03, 4.0, -0.09, 6.0, -0.15, 7.0, -0.03, 3.0
	(*itrue2)[0], (*itrue2)[1], (*itrue2)[2], (*itrue2)[3], (*itrue2)[4] = 0, 1, 2, 2, 3
	//*     .. Executable Statements ..
	for (*incx) = 1; (*incx) <= 2; (*incx)++ {
		for (*np1) = 1; (*np1) <= 5; (*np1)++ {
			(*n) = (*np1) - 1
			(*len) = 2 * MAX((*n), 1)
			//*           .. Set vector arguments ..
			for (*i) = 1; (*i) <= (*len); (*i)++ {
				(*sx)[(*i)-1] = (*dv)[(*i)-1][(*np1)-1][(*incx)-1]
			//Label20:
			}
			//*
			if (*icase) == 7 {
				//*              .. Dnrm2 ..
				(*stemp)[0] = (*dtrue1)[(*np1)-1]
				Stest1(Dnrm2(n, sx, incx), &((*stemp)[0]), stemp, sfac)
			} else if (*icase) == 8 {
				//*              .. Dasum ..
				(*stemp)[0] = (*dtrue3)[(*np1)-1]
				Stest1(Dasum(n, sx, incx), &((*stemp)[0]), stemp, sfac)
			} else if (*icase) == 9 {
				//*              .. Dscal ..
				Dscal(n, &((*sa)[((*incx)-1)*5+(*np1)-1]), sx, incx)
				for (*i) = 1; (*i) <= (*len); (*i)++ {
					(*strue)[(*i)-1] = (*dtrue5)[(*i)-1][(*np1)-1][(*incx)-1]
				//Label40:
				}
				Stest(len, sx, strue, strue, sfac)
			} else if (*icase) == 10 {
				//*              .. Idamax ..
				Itest1(Idamax(n, sx, incx), &((*itrue2)[(*np1)-1]))
			} else {
				WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}(), *func() *[]byte{y := []byte(" shouldn't be here in check1"); return &y}())
				panic("")
			}
		//Label60:
		}
	//Label80:
	}
	return
}

func Check2(sfac *float64) {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sa := new(float64)
	i := new(int)
	j := new(int)
	ki := new(int)
	kn := new(int)
	kni := new(int)
	kpar := new(int)
	ksize := new(int)
	lenx := new(int)
	leny := new(int)
	mx := new(int)
	my := new(int)
	dt10x := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt10y := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt7 := func() *[][]float64 {
		arr := make([][]float64, 4)
		for u := 0; u < 4; u++ {
			arr[u] = make([]float64, 4)
		}
		return &arr
	}()
	dt8 := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dx1 := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	dy1 := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	ssize1 := func() *[]float64 {
		arr := make([]float64, 4)
		return &arr
	}()
	ssize2 := func() *[][]float64 {
		arr := make([][]float64, 14)
		for u := 0; u < 14; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	ssize := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	stx := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sty := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sx := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sy := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	dpar := func() *[][]float64 {
		arr := make([][]float64, 5)
		for u := 0; u < 5; u++ {
			arr[u] = make([]float64, 4)
		}
		return &arr
	}()
	dt19x := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 16)
			}
		}
		return &arr
	}()
	dt19xa := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19xb := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19xc := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19xd := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19y := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 16)
			}
		}
		return &arr
	}()
	dt19ya := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19yb := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19yc := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt19yd := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dtemp := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	incxs := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	incys := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	lens := func() *[][]int {
		arr := make([][]int, 4)
		for u := 0; u < 4; u++ {
			arr[u] = make([]int, 2)
		}
		return &arr
	}()
	ns := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*sa) = 0.3
	(*incxs)[0], (*incxs)[1], (*incxs)[2], (*incxs)[3] = 1, 2, -2, -1
	(*incys)[0], (*incys)[1], (*incys)[2], (*incys)[3] = 1, -2, 1, -2
	(*lens)[0][0], (*lens)[1][0], (*lens)[2][0], (*lens)[3][0], (*lens)[0][1], (*lens)[1][1], (*lens)[2][1], (*lens)[3][1] = 1, 1, 2, 4, 1, 1, 3, 7
	(*ns)[0], (*ns)[1], (*ns)[2], (*ns)[3] = 0, 1, 2, 4
	(*dx1)[0], (*dx1)[1], (*dx1)[2], (*dx1)[3], (*dx1)[4], (*dx1)[5], (*dx1)[6] = 0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4
	(*dy1)[0], (*dy1)[1], (*dy1)[2], (*dy1)[3], (*dy1)[4], (*dy1)[5], (*dy1)[6] = 0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8
	(*dt7)[0][0], (*dt7)[1][0], (*dt7)[2][0], (*dt7)[3][0], (*dt7)[0][1], (*dt7)[1][1], (*dt7)[2][1], (*dt7)[3][1], (*dt7)[0][2], (*dt7)[1][2], (*dt7)[2][2], (*dt7)[3][2], (*dt7)[0][3], (*dt7)[1][3], (*dt7)[2][3], (*dt7)[3][3] = 0.0, 0.30, 0.21, 0.62, 0.0, 0.30, -0.07, 0.85, 0.0, 0.30, -0.79, -0.74, 0.0, 0.30, 0.33, 1.27
	(*dt8)[0][0][0], (*dt8)[1][0][0], (*dt8)[2][0][0], (*dt8)[3][0][0], (*dt8)[4][0][0], (*dt8)[5][0][0], (*dt8)[6][0][0], (*dt8)[0][1][0], (*dt8)[1][1][0], (*dt8)[2][1][0], (*dt8)[3][1][0], (*dt8)[4][1][0], (*dt8)[5][1][0], (*dt8)[6][1][0], (*dt8)[0][2][0], (*dt8)[1][2][0], (*dt8)[2][2][0], (*dt8)[3][2][0], (*dt8)[4][2][0], (*dt8)[5][2][0], (*dt8)[6][2][0], (*dt8)[0][3][0], (*dt8)[1][3][0], (*dt8)[2][3][0], (*dt8)[3][3][0], (*dt8)[4][3][0], (*dt8)[5][3][0], (*dt8)[6][3][0], (*dt8)[0][0][1], (*dt8)[1][0][1], (*dt8)[2][0][1], (*dt8)[3][0][1], (*dt8)[4][0][1], (*dt8)[5][0][1], (*dt8)[6][0][1], (*dt8)[0][1][1], (*dt8)[1][1][1], (*dt8)[2][1][1], (*dt8)[3][1][1], (*dt8)[4][1][1], (*dt8)[5][1][1], (*dt8)[6][1][1], (*dt8)[0][2][1], (*dt8)[1][2][1], (*dt8)[2][2][1], (*dt8)[3][2][1], (*dt8)[4][2][1], (*dt8)[5][2][1], (*dt8)[6][2][1], (*dt8)[0][3][1], (*dt8)[1][3][1], (*dt8)[2][3][1], (*dt8)[3][3][1], (*dt8)[4][3][1], (*dt8)[5][3][1], (*dt8)[6][3][1], (*dt8)[0][0][2], (*dt8)[1][0][2], (*dt8)[2][0][2], (*dt8)[3][0][2], (*dt8)[4][0][2], (*dt8)[5][0][2], (*dt8)[6][0][2], (*dt8)[0][1][2], (*dt8)[1][1][2], (*dt8)[2][1][2], (*dt8)[3][1][2], (*dt8)[4][1][2], (*dt8)[5][1][2], (*dt8)[6][1][2], (*dt8)[0][2][2], (*dt8)[1][2][2], (*dt8)[2][2][2], (*dt8)[3][2][2], (*dt8)[4][2][2], (*dt8)[5][2][2], (*dt8)[6][2][2], (*dt8)[0][3][2], (*dt8)[1][3][2], (*dt8)[2][3][2], (*dt8)[3][3][2], (*dt8)[4][3][2], (*dt8)[5][3][2], (*dt8)[6][3][2], (*dt8)[0][0][3], (*dt8)[1][0][3], (*dt8)[2][0][3], (*dt8)[3][0][3], (*dt8)[4][0][3], (*dt8)[5][0][3], (*dt8)[6][0][3], (*dt8)[0][1][3], (*dt8)[1][1][3], (*dt8)[2][1][3], (*dt8)[3][1][3], (*dt8)[4][1][3], (*dt8)[5][1][3], (*dt8)[6][1][3], (*dt8)[0][2][3], (*dt8)[1][2][3], (*dt8)[2][2][3], (*dt8)[3][2][3], (*dt8)[4][2][3], (*dt8)[5][2][3], (*dt8)[6][2][3], (*dt8)[0][3][3], (*dt8)[1][3][3], (*dt8)[2][3][3], (*dt8)[3][3][3], (*dt8)[4][3][3], (*dt8)[5][3][3], (*dt8)[6][3][3] = 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.87, 0.15, 0.94, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.9, 0.48, 0.0, 0.0, 0.0, 0.0, 0.38, -0.9, 0.57, 0.7, -0.75, 0.2, 0.98, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35, -0.72, 0.0, 0.0, 0.0, 0.0, 0.0, 0.38, -0.63, 0.15, 0.88, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.0, 0.0, 0.0, 0.0, 0.68, -0.9, 0.33, 0.7, -0.75, 0.2, 1.04
	(*dt10x)[0][0][0], (*dt10x)[1][0][0], (*dt10x)[2][0][0], (*dt10x)[3][0][0], (*dt10x)[4][0][0], (*dt10x)[5][0][0], (*dt10x)[6][0][0], (*dt10x)[0][1][0], (*dt10x)[1][1][0], (*dt10x)[2][1][0], (*dt10x)[3][1][0], (*dt10x)[4][1][0], (*dt10x)[5][1][0], (*dt10x)[6][1][0], (*dt10x)[0][2][0], (*dt10x)[1][2][0], (*dt10x)[2][2][0], (*dt10x)[3][2][0], (*dt10x)[4][2][0], (*dt10x)[5][2][0], (*dt10x)[6][2][0], (*dt10x)[0][3][0], (*dt10x)[1][3][0], (*dt10x)[2][3][0], (*dt10x)[3][3][0], (*dt10x)[4][3][0], (*dt10x)[5][3][0], (*dt10x)[6][3][0], (*dt10x)[0][0][1], (*dt10x)[1][0][1], (*dt10x)[2][0][1], (*dt10x)[3][0][1], (*dt10x)[4][0][1], (*dt10x)[5][0][1], (*dt10x)[6][0][1], (*dt10x)[0][1][1], (*dt10x)[1][1][1], (*dt10x)[2][1][1], (*dt10x)[3][1][1], (*dt10x)[4][1][1], (*dt10x)[5][1][1], (*dt10x)[6][1][1], (*dt10x)[0][2][1], (*dt10x)[1][2][1], (*dt10x)[2][2][1], (*dt10x)[3][2][1], (*dt10x)[4][2][1], (*dt10x)[5][2][1], (*dt10x)[6][2][1], (*dt10x)[0][3][1], (*dt10x)[1][3][1], (*dt10x)[2][3][1], (*dt10x)[3][3][1], (*dt10x)[4][3][1], (*dt10x)[5][3][1], (*dt10x)[6][3][1], (*dt10x)[0][0][2], (*dt10x)[1][0][2], (*dt10x)[2][0][2], (*dt10x)[3][0][2], (*dt10x)[4][0][2], (*dt10x)[5][0][2], (*dt10x)[6][0][2], (*dt10x)[0][1][2], (*dt10x)[1][1][2], (*dt10x)[2][1][2], (*dt10x)[3][1][2], (*dt10x)[4][1][2], (*dt10x)[5][1][2], (*dt10x)[6][1][2], (*dt10x)[0][2][2], (*dt10x)[1][2][2], (*dt10x)[2][2][2], (*dt10x)[3][2][2], (*dt10x)[4][2][2], (*dt10x)[5][2][2], (*dt10x)[6][2][2], (*dt10x)[0][3][2], (*dt10x)[1][3][2], (*dt10x)[2][3][2], (*dt10x)[3][3][2], (*dt10x)[4][3][2], (*dt10x)[5][3][2], (*dt10x)[6][3][2], (*dt10x)[0][0][3], (*dt10x)[1][0][3], (*dt10x)[2][0][3], (*dt10x)[3][0][3], (*dt10x)[4][0][3], (*dt10x)[5][0][3], (*dt10x)[6][0][3], (*dt10x)[0][1][3], (*dt10x)[1][1][3], (*dt10x)[2][1][3], (*dt10x)[3][1][3], (*dt10x)[4][1][3], (*dt10x)[5][1][3], (*dt10x)[6][1][3], (*dt10x)[0][2][3], (*dt10x)[1][2][3], (*dt10x)[2][2][3], (*dt10x)[3][2][3], (*dt10x)[4][2][3], (*dt10x)[5][2][3], (*dt10x)[6][2][3], (*dt10x)[0][3][3], (*dt10x)[1][3][3], (*dt10x)[2][3][3], (*dt10x)[3][3][3], (*dt10x)[4][3][3], (*dt10x)[5][3][3], (*dt10x)[6][3][3] = 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.8, 0.1, -0.6, 0.8, 0.3, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.9, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.7, 0.1, 0.3, 0.8, -0.9, -0.3, 0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.3, -0.6, 0.8, 0.0, 0.0, 0.0
	(*dt10y)[0][0][0], (*dt10y)[1][0][0], (*dt10y)[2][0][0], (*dt10y)[3][0][0], (*dt10y)[4][0][0], (*dt10y)[5][0][0], (*dt10y)[6][0][0], (*dt10y)[0][1][0], (*dt10y)[1][1][0], (*dt10y)[2][1][0], (*dt10y)[3][1][0], (*dt10y)[4][1][0], (*dt10y)[5][1][0], (*dt10y)[6][1][0], (*dt10y)[0][2][0], (*dt10y)[1][2][0], (*dt10y)[2][2][0], (*dt10y)[3][2][0], (*dt10y)[4][2][0], (*dt10y)[5][2][0], (*dt10y)[6][2][0], (*dt10y)[0][3][0], (*dt10y)[1][3][0], (*dt10y)[2][3][0], (*dt10y)[3][3][0], (*dt10y)[4][3][0], (*dt10y)[5][3][0], (*dt10y)[6][3][0], (*dt10y)[0][0][1], (*dt10y)[1][0][1], (*dt10y)[2][0][1], (*dt10y)[3][0][1], (*dt10y)[4][0][1], (*dt10y)[5][0][1], (*dt10y)[6][0][1], (*dt10y)[0][1][1], (*dt10y)[1][1][1], (*dt10y)[2][1][1], (*dt10y)[3][1][1], (*dt10y)[4][1][1], (*dt10y)[5][1][1], (*dt10y)[6][1][1], (*dt10y)[0][2][1], (*dt10y)[1][2][1], (*dt10y)[2][2][1], (*dt10y)[3][2][1], (*dt10y)[4][2][1], (*dt10y)[5][2][1], (*dt10y)[6][2][1], (*dt10y)[0][3][1], (*dt10y)[1][3][1], (*dt10y)[2][3][1], (*dt10y)[3][3][1], (*dt10y)[4][3][1], (*dt10y)[5][3][1], (*dt10y)[6][3][1], (*dt10y)[0][0][2], (*dt10y)[1][0][2], (*dt10y)[2][0][2], (*dt10y)[3][0][2], (*dt10y)[4][0][2], (*dt10y)[5][0][2], (*dt10y)[6][0][2], (*dt10y)[0][1][2], (*dt10y)[1][1][2], (*dt10y)[2][1][2], (*dt10y)[3][1][2], (*dt10y)[4][1][2], (*dt10y)[5][1][2], (*dt10y)[6][1][2], (*dt10y)[0][2][2], (*dt10y)[1][2][2], (*dt10y)[2][2][2], (*dt10y)[3][2][2], (*dt10y)[4][2][2], (*dt10y)[5][2][2], (*dt10y)[6][2][2], (*dt10y)[0][3][2], (*dt10y)[1][3][2], (*dt10y)[2][3][2], (*dt10y)[3][3][2], (*dt10y)[4][3][2], (*dt10y)[5][3][2], (*dt10y)[6][3][2], (*dt10y)[0][0][3], (*dt10y)[1][0][3], (*dt10y)[2][0][3], (*dt10y)[3][0][3], (*dt10y)[4][0][3], (*dt10y)[5][0][3], (*dt10y)[6][0][3], (*dt10y)[0][1][3], (*dt10y)[1][1][3], (*dt10y)[2][1][3], (*dt10y)[3][1][3], (*dt10y)[4][1][3], (*dt10y)[5][1][3], (*dt10y)[6][1][3], (*dt10y)[0][2][3], (*dt10y)[1][2][3], (*dt10y)[2][2][3], (*dt10y)[3][2][3], (*dt10y)[4][2][3], (*dt10y)[5][2][3], (*dt10y)[6][2][3], (*dt10y)[0][3][3], (*dt10y)[1][3][3], (*dt10y)[2][3][3], (*dt10y)[3][3][3], (*dt10y)[4][3][3], (*dt10y)[5][3][3], (*dt10y)[6][3][3] = 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.9, 0.6, 0.0, 0.0, 0.0, 0.0, -0.4, -0.9, 0.9, 0.7, -0.5, 0.2, 0.6, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.4, 0.9, -0.5, 0.6, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.6, -0.9, 0.1, 0.7, -0.5, 0.2, 0.8
	(*ssize1)[0], (*ssize1)[1], (*ssize1)[2], (*ssize1)[3] = 0.0, 0.3, 1.6, 3.2
	(*ssize2)[0][0], (*ssize2)[1][0], (*ssize2)[2][0], (*ssize2)[3][0], (*ssize2)[4][0], (*ssize2)[5][0], (*ssize2)[6][0], (*ssize2)[7][0], (*ssize2)[8][0], (*ssize2)[9][0], (*ssize2)[10][0], (*ssize2)[11][0], (*ssize2)[12][0], (*ssize2)[13][0], (*ssize2)[0][1], (*ssize2)[1][1], (*ssize2)[2][1], (*ssize2)[3][1], (*ssize2)[4][1], (*ssize2)[5][1], (*ssize2)[6][1], (*ssize2)[7][1], (*ssize2)[8][1], (*ssize2)[9][1], (*ssize2)[10][1], (*ssize2)[11][1], (*ssize2)[12][1], (*ssize2)[13][1] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17
	//*
	//*                         FOR Drotm
	//*
	(*dpar)[0][0], (*dpar)[1][0], (*dpar)[2][0], (*dpar)[3][0], (*dpar)[4][0], (*dpar)[0][1], (*dpar)[1][1], (*dpar)[2][1], (*dpar)[3][1], (*dpar)[4][1], (*dpar)[0][2], (*dpar)[1][2], (*dpar)[2][2], (*dpar)[3][2], (*dpar)[4][2], (*dpar)[0][3], (*dpar)[1][3], (*dpar)[2][3], (*dpar)[3][3], (*dpar)[4][3] = -2., 0., 0., 0., 0., -1., 2., -3., -4., 5., 0., 0., 2., -3., 0., 1., 5., 2., 0., -4.
	//*                        TRUE X RESULTS F0R ROTATIONS Drotm
	(*dt19xa)[0][0][0], (*dt19xa)[1][0][0], (*dt19xa)[2][0][0], (*dt19xa)[3][0][0], (*dt19xa)[4][0][0], (*dt19xa)[5][0][0], (*dt19xa)[6][0][0], (*dt19xa)[0][1][0], (*dt19xa)[1][1][0], (*dt19xa)[2][1][0], (*dt19xa)[3][1][0], (*dt19xa)[4][1][0], (*dt19xa)[5][1][0], (*dt19xa)[6][1][0], (*dt19xa)[0][2][0], (*dt19xa)[1][2][0], (*dt19xa)[2][2][0], (*dt19xa)[3][2][0], (*dt19xa)[4][2][0], (*dt19xa)[5][2][0], (*dt19xa)[6][2][0], (*dt19xa)[0][3][0], (*dt19xa)[1][3][0], (*dt19xa)[2][3][0], (*dt19xa)[3][3][0], (*dt19xa)[4][3][0], (*dt19xa)[5][3][0], (*dt19xa)[6][3][0], (*dt19xa)[0][0][1], (*dt19xa)[1][0][1], (*dt19xa)[2][0][1], (*dt19xa)[3][0][1], (*dt19xa)[4][0][1], (*dt19xa)[5][0][1], (*dt19xa)[6][0][1], (*dt19xa)[0][1][1], (*dt19xa)[1][1][1], (*dt19xa)[2][1][1], (*dt19xa)[3][1][1], (*dt19xa)[4][1][1], (*dt19xa)[5][1][1], (*dt19xa)[6][1][1], (*dt19xa)[0][2][1], (*dt19xa)[1][2][1], (*dt19xa)[2][2][1], (*dt19xa)[3][2][1], (*dt19xa)[4][2][1], (*dt19xa)[5][2][1], (*dt19xa)[6][2][1], (*dt19xa)[0][3][1], (*dt19xa)[1][3][1], (*dt19xa)[2][3][1], (*dt19xa)[3][3][1], (*dt19xa)[4][3][1], (*dt19xa)[5][3][1], (*dt19xa)[6][3][1], (*dt19xa)[0][0][2], (*dt19xa)[1][0][2], (*dt19xa)[2][0][2], (*dt19xa)[3][0][2], (*dt19xa)[4][0][2], (*dt19xa)[5][0][2], (*dt19xa)[6][0][2], (*dt19xa)[0][1][2], (*dt19xa)[1][1][2], (*dt19xa)[2][1][2], (*dt19xa)[3][1][2], (*dt19xa)[4][1][2], (*dt19xa)[5][1][2], (*dt19xa)[6][1][2], (*dt19xa)[0][2][2], (*dt19xa)[1][2][2], (*dt19xa)[2][2][2], (*dt19xa)[3][2][2], (*dt19xa)[4][2][2], (*dt19xa)[5][2][2], (*dt19xa)[6][2][2], (*dt19xa)[0][3][2], (*dt19xa)[1][3][2], (*dt19xa)[2][3][2], (*dt19xa)[3][3][2], (*dt19xa)[4][3][2], (*dt19xa)[5][3][2], (*dt19xa)[6][3][2], (*dt19xa)[0][0][3], (*dt19xa)[1][0][3], (*dt19xa)[2][0][3], (*dt19xa)[3][0][3], (*dt19xa)[4][0][3], (*dt19xa)[5][0][3], (*dt19xa)[6][0][3], (*dt19xa)[0][1][3], (*dt19xa)[1][1][3], (*dt19xa)[2][1][3], (*dt19xa)[3][1][3], (*dt19xa)[4][1][3], (*dt19xa)[5][1][3], (*dt19xa)[6][1][3], (*dt19xa)[0][2][3], (*dt19xa)[1][2][3], (*dt19xa)[2][2][3], (*dt19xa)[3][2][3], (*dt19xa)[4][2][3], (*dt19xa)[5][2][3], (*dt19xa)[6][2][3], (*dt19xa)[0][3][3], (*dt19xa)[1][3][3], (*dt19xa)[2][3][3], (*dt19xa)[3][3][3], (*dt19xa)[4][3][3], (*dt19xa)[5][3][3], (*dt19xa)[6][3][3] = .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., -.8, 0., 0., 0., 0., 0., 0., -.9, 0., 0., 0., 0., 0., 0., 3.5, 0., 0., 0., 0., 0., 0., .6, .1, 0., 0., 0., 0., 0., -.8, 3.8, 0., 0., 0., 0., 0., -.9, 2.8, 0., 0., 0., 0., 0., 3.5, -.4, 0., 0., 0., 0., 0., .6, .1, -.5, .8, 0., 0., 0., -.8, 3.8, -2.2, -1.2, 0., 0., 0., -.9, 2.8, -1.4, -1.3, 0., 0., 0., 3.5, -.4, -2.2, 4.7, 0., 0., 0.
	//*
	(*dt19xb)[0][0][0], (*dt19xb)[1][0][0], (*dt19xb)[2][0][0], (*dt19xb)[3][0][0], (*dt19xb)[4][0][0], (*dt19xb)[5][0][0], (*dt19xb)[6][0][0], (*dt19xb)[0][1][0], (*dt19xb)[1][1][0], (*dt19xb)[2][1][0], (*dt19xb)[3][1][0], (*dt19xb)[4][1][0], (*dt19xb)[5][1][0], (*dt19xb)[6][1][0], (*dt19xb)[0][2][0], (*dt19xb)[1][2][0], (*dt19xb)[2][2][0], (*dt19xb)[3][2][0], (*dt19xb)[4][2][0], (*dt19xb)[5][2][0], (*dt19xb)[6][2][0], (*dt19xb)[0][3][0], (*dt19xb)[1][3][0], (*dt19xb)[2][3][0], (*dt19xb)[3][3][0], (*dt19xb)[4][3][0], (*dt19xb)[5][3][0], (*dt19xb)[6][3][0], (*dt19xb)[0][0][1], (*dt19xb)[1][0][1], (*dt19xb)[2][0][1], (*dt19xb)[3][0][1], (*dt19xb)[4][0][1], (*dt19xb)[5][0][1], (*dt19xb)[6][0][1], (*dt19xb)[0][1][1], (*dt19xb)[1][1][1], (*dt19xb)[2][1][1], (*dt19xb)[3][1][1], (*dt19xb)[4][1][1], (*dt19xb)[5][1][1], (*dt19xb)[6][1][1], (*dt19xb)[0][2][1], (*dt19xb)[1][2][1], (*dt19xb)[2][2][1], (*dt19xb)[3][2][1], (*dt19xb)[4][2][1], (*dt19xb)[5][2][1], (*dt19xb)[6][2][1], (*dt19xb)[0][3][1], (*dt19xb)[1][3][1], (*dt19xb)[2][3][1], (*dt19xb)[3][3][1], (*dt19xb)[4][3][1], (*dt19xb)[5][3][1], (*dt19xb)[6][3][1], (*dt19xb)[0][0][2], (*dt19xb)[1][0][2], (*dt19xb)[2][0][2], (*dt19xb)[3][0][2], (*dt19xb)[4][0][2], (*dt19xb)[5][0][2], (*dt19xb)[6][0][2], (*dt19xb)[0][1][2], (*dt19xb)[1][1][2], (*dt19xb)[2][1][2], (*dt19xb)[3][1][2], (*dt19xb)[4][1][2], (*dt19xb)[5][1][2], (*dt19xb)[6][1][2], (*dt19xb)[0][2][2], (*dt19xb)[1][2][2], (*dt19xb)[2][2][2], (*dt19xb)[3][2][2], (*dt19xb)[4][2][2], (*dt19xb)[5][2][2], (*dt19xb)[6][2][2], (*dt19xb)[0][3][2], (*dt19xb)[1][3][2], (*dt19xb)[2][3][2], (*dt19xb)[3][3][2], (*dt19xb)[4][3][2], (*dt19xb)[5][3][2], (*dt19xb)[6][3][2], (*dt19xb)[0][0][3], (*dt19xb)[1][0][3], (*dt19xb)[2][0][3], (*dt19xb)[3][0][3], (*dt19xb)[4][0][3], (*dt19xb)[5][0][3], (*dt19xb)[6][0][3], (*dt19xb)[0][1][3], (*dt19xb)[1][1][3], (*dt19xb)[2][1][3], (*dt19xb)[3][1][3], (*dt19xb)[4][1][3], (*dt19xb)[5][1][3], (*dt19xb)[6][1][3], (*dt19xb)[0][2][3], (*dt19xb)[1][2][3], (*dt19xb)[2][2][3], (*dt19xb)[3][2][3], (*dt19xb)[4][2][3], (*dt19xb)[5][2][3], (*dt19xb)[6][2][3], (*dt19xb)[0][3][3], (*dt19xb)[1][3][3], (*dt19xb)[2][3][3], (*dt19xb)[3][3][3], (*dt19xb)[4][3][3], (*dt19xb)[5][3][3], (*dt19xb)[6][3][3] = .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., -.8, 0., 0., 0., 0., 0., 0., -.9, 0., 0., 0., 0., 0., 0., 3.5, 0., 0., 0., 0., 0., 0., .6, .1, -.5, 0., 0., 0., 0., 0., .1, -3.0, 0., 0., 0., 0., -.3, .1, -2.0, 0., 0., 0., 0., 3.3, .1, -2.0, 0., 0., 0., 0., .6, .1, -.5, .8, .9, -.3, -.4, -2.0, .1, 1.4, .8, .6, -.3, -2.8, -1.8, .1, 1.3, .8, 0., -.3, -1.9, 3.8, .1, -3.1, .8, 4.8, -.3, -1.5
	//*
	(*dt19xc)[0][0][0], (*dt19xc)[1][0][0], (*dt19xc)[2][0][0], (*dt19xc)[3][0][0], (*dt19xc)[4][0][0], (*dt19xc)[5][0][0], (*dt19xc)[6][0][0], (*dt19xc)[0][1][0], (*dt19xc)[1][1][0], (*dt19xc)[2][1][0], (*dt19xc)[3][1][0], (*dt19xc)[4][1][0], (*dt19xc)[5][1][0], (*dt19xc)[6][1][0], (*dt19xc)[0][2][0], (*dt19xc)[1][2][0], (*dt19xc)[2][2][0], (*dt19xc)[3][2][0], (*dt19xc)[4][2][0], (*dt19xc)[5][2][0], (*dt19xc)[6][2][0], (*dt19xc)[0][3][0], (*dt19xc)[1][3][0], (*dt19xc)[2][3][0], (*dt19xc)[3][3][0], (*dt19xc)[4][3][0], (*dt19xc)[5][3][0], (*dt19xc)[6][3][0], (*dt19xc)[0][0][1], (*dt19xc)[1][0][1], (*dt19xc)[2][0][1], (*dt19xc)[3][0][1], (*dt19xc)[4][0][1], (*dt19xc)[5][0][1], (*dt19xc)[6][0][1], (*dt19xc)[0][1][1], (*dt19xc)[1][1][1], (*dt19xc)[2][1][1], (*dt19xc)[3][1][1], (*dt19xc)[4][1][1], (*dt19xc)[5][1][1], (*dt19xc)[6][1][1], (*dt19xc)[0][2][1], (*dt19xc)[1][2][1], (*dt19xc)[2][2][1], (*dt19xc)[3][2][1], (*dt19xc)[4][2][1], (*dt19xc)[5][2][1], (*dt19xc)[6][2][1], (*dt19xc)[0][3][1], (*dt19xc)[1][3][1], (*dt19xc)[2][3][1], (*dt19xc)[3][3][1], (*dt19xc)[4][3][1], (*dt19xc)[5][3][1], (*dt19xc)[6][3][1], (*dt19xc)[0][0][2], (*dt19xc)[1][0][2], (*dt19xc)[2][0][2], (*dt19xc)[3][0][2], (*dt19xc)[4][0][2], (*dt19xc)[5][0][2], (*dt19xc)[6][0][2], (*dt19xc)[0][1][2], (*dt19xc)[1][1][2], (*dt19xc)[2][1][2], (*dt19xc)[3][1][2], (*dt19xc)[4][1][2], (*dt19xc)[5][1][2], (*dt19xc)[6][1][2], (*dt19xc)[0][2][2], (*dt19xc)[1][2][2], (*dt19xc)[2][2][2], (*dt19xc)[3][2][2], (*dt19xc)[4][2][2], (*dt19xc)[5][2][2], (*dt19xc)[6][2][2], (*dt19xc)[0][3][2], (*dt19xc)[1][3][2], (*dt19xc)[2][3][2], (*dt19xc)[3][3][2], (*dt19xc)[4][3][2], (*dt19xc)[5][3][2], (*dt19xc)[6][3][2], (*dt19xc)[0][0][3], (*dt19xc)[1][0][3], (*dt19xc)[2][0][3], (*dt19xc)[3][0][3], (*dt19xc)[4][0][3], (*dt19xc)[5][0][3], (*dt19xc)[6][0][3], (*dt19xc)[0][1][3], (*dt19xc)[1][1][3], (*dt19xc)[2][1][3], (*dt19xc)[3][1][3], (*dt19xc)[4][1][3], (*dt19xc)[5][1][3], (*dt19xc)[6][1][3], (*dt19xc)[0][2][3], (*dt19xc)[1][2][3], (*dt19xc)[2][2][3], (*dt19xc)[3][2][3], (*dt19xc)[4][2][3], (*dt19xc)[5][2][3], (*dt19xc)[6][2][3], (*dt19xc)[0][3][3], (*dt19xc)[1][3][3], (*dt19xc)[2][3][3], (*dt19xc)[3][3][3], (*dt19xc)[4][3][3], (*dt19xc)[5][3][3], (*dt19xc)[6][3][3] = .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., -.8, 0., 0., 0., 0., 0., 0., -.9, 0., 0., 0., 0., 0., 0., 3.5, 0., 0., 0., 0., 0., 0., .6, .1, -.5, 0., 0., 0., 0., 4.8, .1, -3.0, 0., 0., 0., 0., 3.3, .1, -2.0, 0., 0., 0., 0., 2.1, .1, -2.0, 0., 0., 0., 0., .6, .1, -.5, .8, .9, -.3, -.4, -1.6, .1, -2.2, .8, 5.4, -.3, -2.8, -1.5, .1, -1.4, .8, 3.6, -.3, -1.9, 3.7, .1, -2.2, .8, 3.6, -.3, -1.5
	//*
	(*dt19xd)[0][0][0], (*dt19xd)[1][0][0], (*dt19xd)[2][0][0], (*dt19xd)[3][0][0], (*dt19xd)[4][0][0], (*dt19xd)[5][0][0], (*dt19xd)[6][0][0], (*dt19xd)[0][1][0], (*dt19xd)[1][1][0], (*dt19xd)[2][1][0], (*dt19xd)[3][1][0], (*dt19xd)[4][1][0], (*dt19xd)[5][1][0], (*dt19xd)[6][1][0], (*dt19xd)[0][2][0], (*dt19xd)[1][2][0], (*dt19xd)[2][2][0], (*dt19xd)[3][2][0], (*dt19xd)[4][2][0], (*dt19xd)[5][2][0], (*dt19xd)[6][2][0], (*dt19xd)[0][3][0], (*dt19xd)[1][3][0], (*dt19xd)[2][3][0], (*dt19xd)[3][3][0], (*dt19xd)[4][3][0], (*dt19xd)[5][3][0], (*dt19xd)[6][3][0], (*dt19xd)[0][0][1], (*dt19xd)[1][0][1], (*dt19xd)[2][0][1], (*dt19xd)[3][0][1], (*dt19xd)[4][0][1], (*dt19xd)[5][0][1], (*dt19xd)[6][0][1], (*dt19xd)[0][1][1], (*dt19xd)[1][1][1], (*dt19xd)[2][1][1], (*dt19xd)[3][1][1], (*dt19xd)[4][1][1], (*dt19xd)[5][1][1], (*dt19xd)[6][1][1], (*dt19xd)[0][2][1], (*dt19xd)[1][2][1], (*dt19xd)[2][2][1], (*dt19xd)[3][2][1], (*dt19xd)[4][2][1], (*dt19xd)[5][2][1], (*dt19xd)[6][2][1], (*dt19xd)[0][3][1], (*dt19xd)[1][3][1], (*dt19xd)[2][3][1], (*dt19xd)[3][3][1], (*dt19xd)[4][3][1], (*dt19xd)[5][3][1], (*dt19xd)[6][3][1], (*dt19xd)[0][0][2], (*dt19xd)[1][0][2], (*dt19xd)[2][0][2], (*dt19xd)[3][0][2], (*dt19xd)[4][0][2], (*dt19xd)[5][0][2], (*dt19xd)[6][0][2], (*dt19xd)[0][1][2], (*dt19xd)[1][1][2], (*dt19xd)[2][1][2], (*dt19xd)[3][1][2], (*dt19xd)[4][1][2], (*dt19xd)[5][1][2], (*dt19xd)[6][1][2], (*dt19xd)[0][2][2], (*dt19xd)[1][2][2], (*dt19xd)[2][2][2], (*dt19xd)[3][2][2], (*dt19xd)[4][2][2], (*dt19xd)[5][2][2], (*dt19xd)[6][2][2], (*dt19xd)[0][3][2], (*dt19xd)[1][3][2], (*dt19xd)[2][3][2], (*dt19xd)[3][3][2], (*dt19xd)[4][3][2], (*dt19xd)[5][3][2], (*dt19xd)[6][3][2], (*dt19xd)[0][0][3], (*dt19xd)[1][0][3], (*dt19xd)[2][0][3], (*dt19xd)[3][0][3], (*dt19xd)[4][0][3], (*dt19xd)[5][0][3], (*dt19xd)[6][0][3], (*dt19xd)[0][1][3], (*dt19xd)[1][1][3], (*dt19xd)[2][1][3], (*dt19xd)[3][1][3], (*dt19xd)[4][1][3], (*dt19xd)[5][1][3], (*dt19xd)[6][1][3], (*dt19xd)[0][2][3], (*dt19xd)[1][2][3], (*dt19xd)[2][2][3], (*dt19xd)[3][2][3], (*dt19xd)[4][2][3], (*dt19xd)[5][2][3], (*dt19xd)[6][2][3], (*dt19xd)[0][3][3], (*dt19xd)[1][3][3], (*dt19xd)[2][3][3], (*dt19xd)[3][3][3], (*dt19xd)[4][3][3], (*dt19xd)[5][3][3], (*dt19xd)[6][3][3] = .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., .6, 0., 0., 0., 0., 0., 0., -.8, 0., 0., 0., 0., 0., 0., -.9, 0., 0., 0., 0., 0., 0., 3.5, 0., 0., 0., 0., 0., 0., .6, .1, 0., 0., 0., 0., 0., -.8, -1.0, 0., 0., 0., 0., 0., -.9, -.8, 0., 0., 0., 0., 0., 3.5, .8, 0., 0., 0., 0., 0., .6, .1, -.5, .8, 0., 0., 0., -.8, -1.0, 1.4, -1.6, 0., 0., 0., -.9, -.8, 1.3, -1.6, 0., 0., 0., 3.5, .8, -3.1, 4.8, 0., 0., 0.
	//*                        TRUE Y RESULTS FOR ROTATIONS Drotm
	(*dt19ya)[0][0][0], (*dt19ya)[1][0][0], (*dt19ya)[2][0][0], (*dt19ya)[3][0][0], (*dt19ya)[4][0][0], (*dt19ya)[5][0][0], (*dt19ya)[6][0][0], (*dt19ya)[0][1][0], (*dt19ya)[1][1][0], (*dt19ya)[2][1][0], (*dt19ya)[3][1][0], (*dt19ya)[4][1][0], (*dt19ya)[5][1][0], (*dt19ya)[6][1][0], (*dt19ya)[0][2][0], (*dt19ya)[1][2][0], (*dt19ya)[2][2][0], (*dt19ya)[3][2][0], (*dt19ya)[4][2][0], (*dt19ya)[5][2][0], (*dt19ya)[6][2][0], (*dt19ya)[0][3][0], (*dt19ya)[1][3][0], (*dt19ya)[2][3][0], (*dt19ya)[3][3][0], (*dt19ya)[4][3][0], (*dt19ya)[5][3][0], (*dt19ya)[6][3][0], (*dt19ya)[0][0][1], (*dt19ya)[1][0][1], (*dt19ya)[2][0][1], (*dt19ya)[3][0][1], (*dt19ya)[4][0][1], (*dt19ya)[5][0][1], (*dt19ya)[6][0][1], (*dt19ya)[0][1][1], (*dt19ya)[1][1][1], (*dt19ya)[2][1][1], (*dt19ya)[3][1][1], (*dt19ya)[4][1][1], (*dt19ya)[5][1][1], (*dt19ya)[6][1][1], (*dt19ya)[0][2][1], (*dt19ya)[1][2][1], (*dt19ya)[2][2][1], (*dt19ya)[3][2][1], (*dt19ya)[4][2][1], (*dt19ya)[5][2][1], (*dt19ya)[6][2][1], (*dt19ya)[0][3][1], (*dt19ya)[1][3][1], (*dt19ya)[2][3][1], (*dt19ya)[3][3][1], (*dt19ya)[4][3][1], (*dt19ya)[5][3][1], (*dt19ya)[6][3][1], (*dt19ya)[0][0][2], (*dt19ya)[1][0][2], (*dt19ya)[2][0][2], (*dt19ya)[3][0][2], (*dt19ya)[4][0][2], (*dt19ya)[5][0][2], (*dt19ya)[6][0][2], (*dt19ya)[0][1][2], (*dt19ya)[1][1][2], (*dt19ya)[2][1][2], (*dt19ya)[3][1][2], (*dt19ya)[4][1][2], (*dt19ya)[5][1][2], (*dt19ya)[6][1][2], (*dt19ya)[0][2][2], (*dt19ya)[1][2][2], (*dt19ya)[2][2][2], (*dt19ya)[3][2][2], (*dt19ya)[4][2][2], (*dt19ya)[5][2][2], (*dt19ya)[6][2][2], (*dt19ya)[0][3][2], (*dt19ya)[1][3][2], (*dt19ya)[2][3][2], (*dt19ya)[3][3][2], (*dt19ya)[4][3][2], (*dt19ya)[5][3][2], (*dt19ya)[6][3][2], (*dt19ya)[0][0][3], (*dt19ya)[1][0][3], (*dt19ya)[2][0][3], (*dt19ya)[3][0][3], (*dt19ya)[4][0][3], (*dt19ya)[5][0][3], (*dt19ya)[6][0][3], (*dt19ya)[0][1][3], (*dt19ya)[1][1][3], (*dt19ya)[2][1][3], (*dt19ya)[3][1][3], (*dt19ya)[4][1][3], (*dt19ya)[5][1][3], (*dt19ya)[6][1][3], (*dt19ya)[0][2][3], (*dt19ya)[1][2][3], (*dt19ya)[2][2][3], (*dt19ya)[3][2][3], (*dt19ya)[4][2][3], (*dt19ya)[5][2][3], (*dt19ya)[6][2][3], (*dt19ya)[0][3][3], (*dt19ya)[1][3][3], (*dt19ya)[2][3][3], (*dt19ya)[3][3][3], (*dt19ya)[4][3][3], (*dt19ya)[5][3][3], (*dt19ya)[6][3][3] = .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .7, 0., 0., 0., 0., 0., 0., 1.7, 0., 0., 0., 0., 0., 0., -2.6, 0., 0., 0., 0., 0., 0., .5, -.9, 0., 0., 0., 0., 0., .7, -4.8, 0., 0., 0., 0., 0., 1.7, -.7, 0., 0., 0., 0., 0., -2.6, 3.5, 0., 0., 0., 0., 0., .5, -.9, .3, .7, 0., 0., 0., .7, -4.8, 3.0, 1.1, 0., 0., 0., 1.7, -.7, -.7, 2.3, 0., 0., 0., -2.6, 3.5, -.7, -3.6, 0., 0., 0.
	//*
	(*dt19yb)[0][0][0], (*dt19yb)[1][0][0], (*dt19yb)[2][0][0], (*dt19yb)[3][0][0], (*dt19yb)[4][0][0], (*dt19yb)[5][0][0], (*dt19yb)[6][0][0], (*dt19yb)[0][1][0], (*dt19yb)[1][1][0], (*dt19yb)[2][1][0], (*dt19yb)[3][1][0], (*dt19yb)[4][1][0], (*dt19yb)[5][1][0], (*dt19yb)[6][1][0], (*dt19yb)[0][2][0], (*dt19yb)[1][2][0], (*dt19yb)[2][2][0], (*dt19yb)[3][2][0], (*dt19yb)[4][2][0], (*dt19yb)[5][2][0], (*dt19yb)[6][2][0], (*dt19yb)[0][3][0], (*dt19yb)[1][3][0], (*dt19yb)[2][3][0], (*dt19yb)[3][3][0], (*dt19yb)[4][3][0], (*dt19yb)[5][3][0], (*dt19yb)[6][3][0], (*dt19yb)[0][0][1], (*dt19yb)[1][0][1], (*dt19yb)[2][0][1], (*dt19yb)[3][0][1], (*dt19yb)[4][0][1], (*dt19yb)[5][0][1], (*dt19yb)[6][0][1], (*dt19yb)[0][1][1], (*dt19yb)[1][1][1], (*dt19yb)[2][1][1], (*dt19yb)[3][1][1], (*dt19yb)[4][1][1], (*dt19yb)[5][1][1], (*dt19yb)[6][1][1], (*dt19yb)[0][2][1], (*dt19yb)[1][2][1], (*dt19yb)[2][2][1], (*dt19yb)[3][2][1], (*dt19yb)[4][2][1], (*dt19yb)[5][2][1], (*dt19yb)[6][2][1], (*dt19yb)[0][3][1], (*dt19yb)[1][3][1], (*dt19yb)[2][3][1], (*dt19yb)[3][3][1], (*dt19yb)[4][3][1], (*dt19yb)[5][3][1], (*dt19yb)[6][3][1], (*dt19yb)[0][0][2], (*dt19yb)[1][0][2], (*dt19yb)[2][0][2], (*dt19yb)[3][0][2], (*dt19yb)[4][0][2], (*dt19yb)[5][0][2], (*dt19yb)[6][0][2], (*dt19yb)[0][1][2], (*dt19yb)[1][1][2], (*dt19yb)[2][1][2], (*dt19yb)[3][1][2], (*dt19yb)[4][1][2], (*dt19yb)[5][1][2], (*dt19yb)[6][1][2], (*dt19yb)[0][2][2], (*dt19yb)[1][2][2], (*dt19yb)[2][2][2], (*dt19yb)[3][2][2], (*dt19yb)[4][2][2], (*dt19yb)[5][2][2], (*dt19yb)[6][2][2], (*dt19yb)[0][3][2], (*dt19yb)[1][3][2], (*dt19yb)[2][3][2], (*dt19yb)[3][3][2], (*dt19yb)[4][3][2], (*dt19yb)[5][3][2], (*dt19yb)[6][3][2], (*dt19yb)[0][0][3], (*dt19yb)[1][0][3], (*dt19yb)[2][0][3], (*dt19yb)[3][0][3], (*dt19yb)[4][0][3], (*dt19yb)[5][0][3], (*dt19yb)[6][0][3], (*dt19yb)[0][1][3], (*dt19yb)[1][1][3], (*dt19yb)[2][1][3], (*dt19yb)[3][1][3], (*dt19yb)[4][1][3], (*dt19yb)[5][1][3], (*dt19yb)[6][1][3], (*dt19yb)[0][2][3], (*dt19yb)[1][2][3], (*dt19yb)[2][2][3], (*dt19yb)[3][2][3], (*dt19yb)[4][2][3], (*dt19yb)[5][2][3], (*dt19yb)[6][2][3], (*dt19yb)[0][3][3], (*dt19yb)[1][3][3], (*dt19yb)[2][3][3], (*dt19yb)[3][3][3], (*dt19yb)[4][3][3], (*dt19yb)[5][3][3], (*dt19yb)[6][3][3] = .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .7, 0., 0., 0., 0., 0., 0., 1.7, 0., 0., 0., 0., 0., 0., -2.6, 0., 0., 0., 0., 0., 0., .5, -.9, .3, 0., 0., 0., 0., 4.0, -.9, -.3, 0., 0., 0., 0., -.5, -.9, 1.5, 0., 0., 0., 0., -1.5, -.9, -1.8, 0., 0., 0., 0., .5, -.9, .3, .7, -.6, .2, .8, 3.7, -.9, -1.2, .7, -1.5, .2, 2.2, -.3, -.9, 2.1, .7, -1.6, .2, 2.0, -1.6, -.9, -2.1, .7, 2.9, .2, -3.8
	//*
	(*dt19yc)[0][0][0], (*dt19yc)[1][0][0], (*dt19yc)[2][0][0], (*dt19yc)[3][0][0], (*dt19yc)[4][0][0], (*dt19yc)[5][0][0], (*dt19yc)[6][0][0], (*dt19yc)[0][1][0], (*dt19yc)[1][1][0], (*dt19yc)[2][1][0], (*dt19yc)[3][1][0], (*dt19yc)[4][1][0], (*dt19yc)[5][1][0], (*dt19yc)[6][1][0], (*dt19yc)[0][2][0], (*dt19yc)[1][2][0], (*dt19yc)[2][2][0], (*dt19yc)[3][2][0], (*dt19yc)[4][2][0], (*dt19yc)[5][2][0], (*dt19yc)[6][2][0], (*dt19yc)[0][3][0], (*dt19yc)[1][3][0], (*dt19yc)[2][3][0], (*dt19yc)[3][3][0], (*dt19yc)[4][3][0], (*dt19yc)[5][3][0], (*dt19yc)[6][3][0], (*dt19yc)[0][0][1], (*dt19yc)[1][0][1], (*dt19yc)[2][0][1], (*dt19yc)[3][0][1], (*dt19yc)[4][0][1], (*dt19yc)[5][0][1], (*dt19yc)[6][0][1], (*dt19yc)[0][1][1], (*dt19yc)[1][1][1], (*dt19yc)[2][1][1], (*dt19yc)[3][1][1], (*dt19yc)[4][1][1], (*dt19yc)[5][1][1], (*dt19yc)[6][1][1], (*dt19yc)[0][2][1], (*dt19yc)[1][2][1], (*dt19yc)[2][2][1], (*dt19yc)[3][2][1], (*dt19yc)[4][2][1], (*dt19yc)[5][2][1], (*dt19yc)[6][2][1], (*dt19yc)[0][3][1], (*dt19yc)[1][3][1], (*dt19yc)[2][3][1], (*dt19yc)[3][3][1], (*dt19yc)[4][3][1], (*dt19yc)[5][3][1], (*dt19yc)[6][3][1], (*dt19yc)[0][0][2], (*dt19yc)[1][0][2], (*dt19yc)[2][0][2], (*dt19yc)[3][0][2], (*dt19yc)[4][0][2], (*dt19yc)[5][0][2], (*dt19yc)[6][0][2], (*dt19yc)[0][1][2], (*dt19yc)[1][1][2], (*dt19yc)[2][1][2], (*dt19yc)[3][1][2], (*dt19yc)[4][1][2], (*dt19yc)[5][1][2], (*dt19yc)[6][1][2], (*dt19yc)[0][2][2], (*dt19yc)[1][2][2], (*dt19yc)[2][2][2], (*dt19yc)[3][2][2], (*dt19yc)[4][2][2], (*dt19yc)[5][2][2], (*dt19yc)[6][2][2], (*dt19yc)[0][3][2], (*dt19yc)[1][3][2], (*dt19yc)[2][3][2], (*dt19yc)[3][3][2], (*dt19yc)[4][3][2], (*dt19yc)[5][3][2], (*dt19yc)[6][3][2], (*dt19yc)[0][0][3], (*dt19yc)[1][0][3], (*dt19yc)[2][0][3], (*dt19yc)[3][0][3], (*dt19yc)[4][0][3], (*dt19yc)[5][0][3], (*dt19yc)[6][0][3], (*dt19yc)[0][1][3], (*dt19yc)[1][1][3], (*dt19yc)[2][1][3], (*dt19yc)[3][1][3], (*dt19yc)[4][1][3], (*dt19yc)[5][1][3], (*dt19yc)[6][1][3], (*dt19yc)[0][2][3], (*dt19yc)[1][2][3], (*dt19yc)[2][2][3], (*dt19yc)[3][2][3], (*dt19yc)[4][2][3], (*dt19yc)[5][2][3], (*dt19yc)[6][2][3], (*dt19yc)[0][3][3], (*dt19yc)[1][3][3], (*dt19yc)[2][3][3], (*dt19yc)[3][3][3], (*dt19yc)[4][3][3], (*dt19yc)[5][3][3], (*dt19yc)[6][3][3] = .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .7, 0., 0., 0., 0., 0., 0., 1.7, 0., 0., 0., 0., 0., 0., -2.6, 0., 0., 0., 0., 0., 0., .5, -.9, 0., 0., 0., 0., 0., 4.0, -6.3, 0., 0., 0., 0., 0., -.5, .3, 0., 0., 0., 0., 0., -1.5, 3.0, 0., 0., 0., 0., 0., .5, -.9, .3, .7, 0., 0., 0., 3.7, -7.2, 3.0, 1.7, 0., 0., 0., -.3, .9, -.7, 1.9, 0., 0., 0., -1.6, 2.7, -.7, -3.4, 0., 0., 0.
	//*
	(*dt19yd)[0][0][0], (*dt19yd)[1][0][0], (*dt19yd)[2][0][0], (*dt19yd)[3][0][0], (*dt19yd)[4][0][0], (*dt19yd)[5][0][0], (*dt19yd)[6][0][0], (*dt19yd)[0][1][0], (*dt19yd)[1][1][0], (*dt19yd)[2][1][0], (*dt19yd)[3][1][0], (*dt19yd)[4][1][0], (*dt19yd)[5][1][0], (*dt19yd)[6][1][0], (*dt19yd)[0][2][0], (*dt19yd)[1][2][0], (*dt19yd)[2][2][0], (*dt19yd)[3][2][0], (*dt19yd)[4][2][0], (*dt19yd)[5][2][0], (*dt19yd)[6][2][0], (*dt19yd)[0][3][0], (*dt19yd)[1][3][0], (*dt19yd)[2][3][0], (*dt19yd)[3][3][0], (*dt19yd)[4][3][0], (*dt19yd)[5][3][0], (*dt19yd)[6][3][0], (*dt19yd)[0][0][1], (*dt19yd)[1][0][1], (*dt19yd)[2][0][1], (*dt19yd)[3][0][1], (*dt19yd)[4][0][1], (*dt19yd)[5][0][1], (*dt19yd)[6][0][1], (*dt19yd)[0][1][1], (*dt19yd)[1][1][1], (*dt19yd)[2][1][1], (*dt19yd)[3][1][1], (*dt19yd)[4][1][1], (*dt19yd)[5][1][1], (*dt19yd)[6][1][1], (*dt19yd)[0][2][1], (*dt19yd)[1][2][1], (*dt19yd)[2][2][1], (*dt19yd)[3][2][1], (*dt19yd)[4][2][1], (*dt19yd)[5][2][1], (*dt19yd)[6][2][1], (*dt19yd)[0][3][1], (*dt19yd)[1][3][1], (*dt19yd)[2][3][1], (*dt19yd)[3][3][1], (*dt19yd)[4][3][1], (*dt19yd)[5][3][1], (*dt19yd)[6][3][1], (*dt19yd)[0][0][2], (*dt19yd)[1][0][2], (*dt19yd)[2][0][2], (*dt19yd)[3][0][2], (*dt19yd)[4][0][2], (*dt19yd)[5][0][2], (*dt19yd)[6][0][2], (*dt19yd)[0][1][2], (*dt19yd)[1][1][2], (*dt19yd)[2][1][2], (*dt19yd)[3][1][2], (*dt19yd)[4][1][2], (*dt19yd)[5][1][2], (*dt19yd)[6][1][2], (*dt19yd)[0][2][2], (*dt19yd)[1][2][2], (*dt19yd)[2][2][2], (*dt19yd)[3][2][2], (*dt19yd)[4][2][2], (*dt19yd)[5][2][2], (*dt19yd)[6][2][2], (*dt19yd)[0][3][2], (*dt19yd)[1][3][2], (*dt19yd)[2][3][2], (*dt19yd)[3][3][2], (*dt19yd)[4][3][2], (*dt19yd)[5][3][2], (*dt19yd)[6][3][2], (*dt19yd)[0][0][3], (*dt19yd)[1][0][3], (*dt19yd)[2][0][3], (*dt19yd)[3][0][3], (*dt19yd)[4][0][3], (*dt19yd)[5][0][3], (*dt19yd)[6][0][3], (*dt19yd)[0][1][3], (*dt19yd)[1][1][3], (*dt19yd)[2][1][3], (*dt19yd)[3][1][3], (*dt19yd)[4][1][3], (*dt19yd)[5][1][3], (*dt19yd)[6][1][3], (*dt19yd)[0][2][3], (*dt19yd)[1][2][3], (*dt19yd)[2][2][3], (*dt19yd)[3][2][3], (*dt19yd)[4][2][3], (*dt19yd)[5][2][3], (*dt19yd)[6][2][3], (*dt19yd)[0][3][3], (*dt19yd)[1][3][3], (*dt19yd)[2][3][3], (*dt19yd)[3][3][3], (*dt19yd)[4][3][3], (*dt19yd)[5][3][3], (*dt19yd)[6][3][3] = .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .5, 0., 0., 0., 0., 0., 0., .7, 0., 0., 0., 0., 0., 0., 1.7, 0., 0., 0., 0., 0., 0., -2.6, 0., 0., 0., 0., 0., 0., .5, -.9, .3, 0., 0., 0., 0., .7, -.9, 1.2, 0., 0., 0., 0., 1.7, -.9, .5, 0., 0., 0., 0., -2.6, -.9, -1.3, 0., 0., 0., 0., .5, -.9, .3, .7, -.6, .2, .8, .7, -.9, 1.2, .7, -1.5, .2, 1.6, 1.7, -.9, .5, .7, -1.6, .2, 2.4, -2.6, -.9, -1.3, .7, 2.9, .2, -4.0
	//*
	//*     .. Executable Statements ..
	//*
	for (*ki) = 1; (*ki) <= 4; (*ki)++ {
		(*incx) = (*incxs)[(*ki)-1]
		(*incy) = (*incys)[(*ki)-1]
		(*mx) = (ABS((*incx)))
		(*my) = (ABS((*incy)))
		//*
		for (*kn) = 1; (*kn) <= 4; (*kn)++ {
			(*n) = (*ns)[(*kn)-1]
			(*ksize) = (MIN(int(2), (*kn)))
			(*lenx) = (*lens)[(*kn)-1][(*mx)-1]
			(*leny) = (*lens)[(*kn)-1][(*my)-1]
			//*           .. Initialize all argument arrays ..
			for (*i) = 1; (*i) <= 7; (*i)++ {
				(*sx)[(*i)-1] = (*dx1)[(*i)-1]
				(*sy)[(*i)-1] = (*dy1)[(*i)-1]
			//Label20:
			}
			//*
			if (*icase) == 1 {
				//*              .. Ddot ..
				Stest1(Ddot(n, sx, incx, sy, incy), &((*dt7)[(*kn)-1][(*ki)-1]), &((*ssize1)[(*kn)-1]), sfac)
			} else if (*icase) == 2 {
				//*              .. Daxpy ..
				Daxpy(n, sa, sx, incx, sy, incy)
				for (*j) = 1; (*j) <= (*leny); (*j)++ {
					(*sty)[(*j)-1] = (*dt8)[(*j)-1][(*kn)-1][(*ki)-1]
				//Label40:
				}
				Stest(leny, sy, sty, &((*ssize2)[0][(*ksize)-1]), sfac)
			} else if (*icase) == 5 {
				//*              .. Dcopy ..
				for (*i) = 1; (*i) <= 7; (*i)++ {
					(*sty)[(*i)-1] = (*dt10y)[(*i)-1][(*kn)-1][(*ki)-1]
				//Label60:
				}
				Dcopy(n, sx, incx, sy, incy)
				Stest(leny, sy, sty, &((*ssize2)[0][0]), func() *float64{y := 1.0; return &y}())
			} else if (*icase) == 6 {
				//*              .. Dswap ..
				Dswap(n, sx, incx, sy, incy)
				for (*i) = 1; (*i) <= 7; (*i)++ {
					(*stx)[(*i)-1] = (*dt10x)[(*i)-1][(*kn)-1][(*ki)-1]
					(*sty)[(*i)-1] = (*dt10y)[(*i)-1][(*kn)-1][(*ki)-1]
				//Label80:
				}
				Stest(lenx, sx, stx, &((*ssize2)[0][0]), func() *float64{y := 1.0; return &y}())
				Stest(leny, sy, sty, &((*ssize2)[0][0]), func() *float64{y := 1.0; return &y}())
			} else if (*icase) == 12 {
				//*              .. Drotm ..
				(*kni) = (*kn) + 4*((*ki)-1)
				for (*kpar) = 1; (*kpar) <= 4; (*kpar)++ {
					for (*i) = 1; (*i) <= 7; (*i)++ {
						(*sx)[(*i)-1] = (*dx1)[(*i)-1]
						(*sy)[(*i)-1] = (*dy1)[(*i)-1]
						(*stx)[(*i)-1] = (*dt19x)[(*i)-1][(*kpar)-1][(*kni)-1]
						(*sty)[(*i)-1] = (*dt19y)[(*i)-1][(*kpar)-1][(*kni)-1]
					}
					//*
					for (*i) = 1; (*i) <= 5; (*i)++ {
						(*dtemp)[(*i)-1] = (*dpar)[(*i)-1][(*kpar)-1]
					}
					//*
					for (*i) = 1; (*i) <= (*lenx); (*i)++ {
						(*ssize)[(*i)-1] = (*stx)[(*i)-1]
					}
					//*                   SEE REMARK ABOVE ABOUT DT11X(1,2,7)
					//*                       AND DT11X(5,3,8).
					if ((*kpar) == 2) && ((*kni) == 7) {
						(*ssize)[0] = 2.4
					}
					if ((*kpar) == 3) && ((*kni) == 8) {
						(*ssize)[4] = 1.8
					}
					//*
					Drotm(n, sx, incx, sy, incy, dtemp)
					Stest(lenx, sx, stx, ssize, sfac)
					Stest(leny, sy, sty, sty, sfac)
				}
			} else if (*icase) == 13 {
				//*              .. Dsdot ..
				TESTDsdot(real(Dsdot(n, real((*sx)), incx, real((*sy)), incy)), real(((*dt7)[(*kn)-1][(*ki)-1])), real(((*ssize1)[(*kn)-1])), func() *float64{y := .3125e-1; return &y}())
			} else {
				WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}(), *func() *[]byte{y := []byte(" shouldn't be here in check2"); return &y}())
				panic("")
			}
		//Label100:
		}
	//Label120:
	}
	return
}

func Check3(sfac *float64) {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sc := new(float64)
	ss := new(float64)
	i := new(int)
	k := new(int)
	ki := new(int)
	kn := new(int)
	ksize := new(int)
	lenx := new(int)
	leny := new(int)
	mx := new(int)
	my := new(int)
	copyx := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	copyy := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	dt9x := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dt9y := func() *[][][]float64 {
		arr := make([][][]float64, 7)
		for u := 0; u < 7; u++ {
			arr[u] = make([][]float64, 4)
			for w := 0; w < 4; w++ {
				arr[u][w] = make([]float64, 4)
			}
		}
		return &arr
	}()
	dx1 := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	dy1 := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	mwpc := func() *[]float64 {
		arr := make([]float64, 11)
		return &arr
	}()
	mwps := func() *[]float64 {
		arr := make([]float64, 11)
		return &arr
	}()
	mwpstx := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	mwpsty := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	mwptx := func() *[][]float64 {
		arr := make([][]float64, 11)
		for u := 0; u < 11; u++ {
			arr[u] = make([]float64, 5)
		}
		return &arr
	}()
	mwpty := func() *[][]float64 {
		arr := make([][]float64, 11)
		for u := 0; u < 11; u++ {
			arr[u] = make([]float64, 5)
		}
		return &arr
	}()
	mwpx := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	mwpy := func() *[]float64 {
		arr := make([]float64, 5)
		return &arr
	}()
	ssize2 := func() *[][]float64 {
		arr := make([][]float64, 14)
		for u := 0; u < 14; u++ {
			arr[u] = make([]float64, 2)
		}
		return &arr
	}()
	stx := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sty := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sx := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	sy := func() *[]float64 {
		arr := make([]float64, 7)
		return &arr
	}()
	incxs := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	incys := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	lens := func() *[][]int {
		arr := make([][]int, 4)
		for u := 0; u < 4; u++ {
			arr[u] = make([]int, 2)
		}
		return &arr
	}()
	mwpinx := func() *[]int {
		arr := make([]int, 11)
		return &arr
	}()
	mwpiny := func() *[]int {
		arr := make([]int, 11)
		return &arr
	}()
	mwpn := func() *[]int {
		arr := make([]int, 11)
		return &arr
	}()
	ns := func() *[]int {
		arr := make([]int, 4)
		return &arr
	}()
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Data statements ..
	(*incxs)[0], (*incxs)[1], (*incxs)[2], (*incxs)[3] = 1, 2, -2, -1
	(*incys)[0], (*incys)[1], (*incys)[2], (*incys)[3] = 1, -2, 1, -2
	(*lens)[0][0], (*lens)[1][0], (*lens)[2][0], (*lens)[3][0], (*lens)[0][1], (*lens)[1][1], (*lens)[2][1], (*lens)[3][1] = 1, 1, 2, 4, 1, 1, 3, 7
	(*ns)[0], (*ns)[1], (*ns)[2], (*ns)[3] = 0, 1, 2, 4
	(*dx1)[0], (*dx1)[1], (*dx1)[2], (*dx1)[3], (*dx1)[4], (*dx1)[5], (*dx1)[6] = 0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4
	(*dy1)[0], (*dy1)[1], (*dy1)[2], (*dy1)[3], (*dy1)[4], (*dy1)[5], (*dy1)[6] = 0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8
	(*sc), (*ss) = 0.8, 0.6
	(*dt9x)[0][0][0], (*dt9x)[1][0][0], (*dt9x)[2][0][0], (*dt9x)[3][0][0], (*dt9x)[4][0][0], (*dt9x)[5][0][0], (*dt9x)[6][0][0], (*dt9x)[0][1][0], (*dt9x)[1][1][0], (*dt9x)[2][1][0], (*dt9x)[3][1][0], (*dt9x)[4][1][0], (*dt9x)[5][1][0], (*dt9x)[6][1][0], (*dt9x)[0][2][0], (*dt9x)[1][2][0], (*dt9x)[2][2][0], (*dt9x)[3][2][0], (*dt9x)[4][2][0], (*dt9x)[5][2][0], (*dt9x)[6][2][0], (*dt9x)[0][3][0], (*dt9x)[1][3][0], (*dt9x)[2][3][0], (*dt9x)[3][3][0], (*dt9x)[4][3][0], (*dt9x)[5][3][0], (*dt9x)[6][3][0], (*dt9x)[0][0][1], (*dt9x)[1][0][1], (*dt9x)[2][0][1], (*dt9x)[3][0][1], (*dt9x)[4][0][1], (*dt9x)[5][0][1], (*dt9x)[6][0][1], (*dt9x)[0][1][1], (*dt9x)[1][1][1], (*dt9x)[2][1][1], (*dt9x)[3][1][1], (*dt9x)[4][1][1], (*dt9x)[5][1][1], (*dt9x)[6][1][1], (*dt9x)[0][2][1], (*dt9x)[1][2][1], (*dt9x)[2][2][1], (*dt9x)[3][2][1], (*dt9x)[4][2][1], (*dt9x)[5][2][1], (*dt9x)[6][2][1], (*dt9x)[0][3][1], (*dt9x)[1][3][1], (*dt9x)[2][3][1], (*dt9x)[3][3][1], (*dt9x)[4][3][1], (*dt9x)[5][3][1], (*dt9x)[6][3][1], (*dt9x)[0][0][2], (*dt9x)[1][0][2], (*dt9x)[2][0][2], (*dt9x)[3][0][2], (*dt9x)[4][0][2], (*dt9x)[5][0][2], (*dt9x)[6][0][2], (*dt9x)[0][1][2], (*dt9x)[1][1][2], (*dt9x)[2][1][2], (*dt9x)[3][1][2], (*dt9x)[4][1][2], (*dt9x)[5][1][2], (*dt9x)[6][1][2], (*dt9x)[0][2][2], (*dt9x)[1][2][2], (*dt9x)[2][2][2], (*dt9x)[3][2][2], (*dt9x)[4][2][2], (*dt9x)[5][2][2], (*dt9x)[6][2][2], (*dt9x)[0][3][2], (*dt9x)[1][3][2], (*dt9x)[2][3][2], (*dt9x)[3][3][2], (*dt9x)[4][3][2], (*dt9x)[5][3][2], (*dt9x)[6][3][2], (*dt9x)[0][0][3], (*dt9x)[1][0][3], (*dt9x)[2][0][3], (*dt9x)[3][0][3], (*dt9x)[4][0][3], (*dt9x)[5][0][3], (*dt9x)[6][0][3], (*dt9x)[0][1][3], (*dt9x)[1][1][3], (*dt9x)[2][1][3], (*dt9x)[3][1][3], (*dt9x)[4][1][3], (*dt9x)[5][1][3], (*dt9x)[6][1][3], (*dt9x)[0][2][3], (*dt9x)[1][2][3], (*dt9x)[2][2][3], (*dt9x)[3][2][3], (*dt9x)[4][2][3], (*dt9x)[5][2][3], (*dt9x)[6][2][3], (*dt9x)[0][3][3], (*dt9x)[1][3][3], (*dt9x)[2][3][3], (*dt9x)[3][3][3], (*dt9x)[4][3][3], (*dt9x)[5][3][3], (*dt9x)[6][3][3] = 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, -0.46, -0.22, 1.06, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.66, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.96, 0.1, -0.76, 0.8, 0.90, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0, 0.90, 0.1, -0.22, 0.8, 0.18, -0.3, -0.02, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.78, 0.26, -0.76, 1.12, 0.0, 0.0, 0.0
	(*dt9y)[0][0][0], (*dt9y)[1][0][0], (*dt9y)[2][0][0], (*dt9y)[3][0][0], (*dt9y)[4][0][0], (*dt9y)[5][0][0], (*dt9y)[6][0][0], (*dt9y)[0][1][0], (*dt9y)[1][1][0], (*dt9y)[2][1][0], (*dt9y)[3][1][0], (*dt9y)[4][1][0], (*dt9y)[5][1][0], (*dt9y)[6][1][0], (*dt9y)[0][2][0], (*dt9y)[1][2][0], (*dt9y)[2][2][0], (*dt9y)[3][2][0], (*dt9y)[4][2][0], (*dt9y)[5][2][0], (*dt9y)[6][2][0], (*dt9y)[0][3][0], (*dt9y)[1][3][0], (*dt9y)[2][3][0], (*dt9y)[3][3][0], (*dt9y)[4][3][0], (*dt9y)[5][3][0], (*dt9y)[6][3][0], (*dt9y)[0][0][1], (*dt9y)[1][0][1], (*dt9y)[2][0][1], (*dt9y)[3][0][1], (*dt9y)[4][0][1], (*dt9y)[5][0][1], (*dt9y)[6][0][1], (*dt9y)[0][1][1], (*dt9y)[1][1][1], (*dt9y)[2][1][1], (*dt9y)[3][1][1], (*dt9y)[4][1][1], (*dt9y)[5][1][1], (*dt9y)[6][1][1], (*dt9y)[0][2][1], (*dt9y)[1][2][1], (*dt9y)[2][2][1], (*dt9y)[3][2][1], (*dt9y)[4][2][1], (*dt9y)[5][2][1], (*dt9y)[6][2][1], (*dt9y)[0][3][1], (*dt9y)[1][3][1], (*dt9y)[2][3][1], (*dt9y)[3][3][1], (*dt9y)[4][3][1], (*dt9y)[5][3][1], (*dt9y)[6][3][1], (*dt9y)[0][0][2], (*dt9y)[1][0][2], (*dt9y)[2][0][2], (*dt9y)[3][0][2], (*dt9y)[4][0][2], (*dt9y)[5][0][2], (*dt9y)[6][0][2], (*dt9y)[0][1][2], (*dt9y)[1][1][2], (*dt9y)[2][1][2], (*dt9y)[3][1][2], (*dt9y)[4][1][2], (*dt9y)[5][1][2], (*dt9y)[6][1][2], (*dt9y)[0][2][2], (*dt9y)[1][2][2], (*dt9y)[2][2][2], (*dt9y)[3][2][2], (*dt9y)[4][2][2], (*dt9y)[5][2][2], (*dt9y)[6][2][2], (*dt9y)[0][3][2], (*dt9y)[1][3][2], (*dt9y)[2][3][2], (*dt9y)[3][3][2], (*dt9y)[4][3][2], (*dt9y)[5][3][2], (*dt9y)[6][3][2], (*dt9y)[0][0][3], (*dt9y)[1][0][3], (*dt9y)[2][0][3], (*dt9y)[3][0][3], (*dt9y)[4][0][3], (*dt9y)[5][0][3], (*dt9y)[6][0][3], (*dt9y)[0][1][3], (*dt9y)[1][1][3], (*dt9y)[2][1][3], (*dt9y)[3][1][3], (*dt9y)[4][1][3], (*dt9y)[5][1][3], (*dt9y)[6][1][3], (*dt9y)[0][2][3], (*dt9y)[1][2][3], (*dt9y)[2][2][3], (*dt9y)[3][2][3], (*dt9y)[4][2][3], (*dt9y)[5][2][3], (*dt9y)[6][2][3], (*dt9y)[0][3][3], (*dt9y)[1][3][3], (*dt9y)[2][3][3], (*dt9y)[3][3][3], (*dt9y)[4][3][3], (*dt9y)[5][3][3], (*dt9y)[6][3][3] = 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.78, 0.54, 0.08, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -0.9, -0.12, 0.0, 0.0, 0.0, 0.0, 0.64, -0.9, -0.30, 0.7, -0.18, 0.2, 0.28, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, -1.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.64, -1.26, 0.54, 0.20, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.0, 0.0, 0.0, 0.0, 0.04, -0.9, 0.18, 0.7, -0.18, 0.2, 0.16
	(*ssize2)[0][0], (*ssize2)[1][0], (*ssize2)[2][0], (*ssize2)[3][0], (*ssize2)[4][0], (*ssize2)[5][0], (*ssize2)[6][0], (*ssize2)[7][0], (*ssize2)[8][0], (*ssize2)[9][0], (*ssize2)[10][0], (*ssize2)[11][0], (*ssize2)[12][0], (*ssize2)[13][0], (*ssize2)[0][1], (*ssize2)[1][1], (*ssize2)[2][1], (*ssize2)[3][1], (*ssize2)[4][1], (*ssize2)[5][1], (*ssize2)[6][1], (*ssize2)[7][1], (*ssize2)[8][1], (*ssize2)[9][1], (*ssize2)[10][1], (*ssize2)[11][1], (*ssize2)[12][1], (*ssize2)[13][1] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17
	//*     .. Executable Statements ..
	//*
	for (*ki) = 1; (*ki) <= 4; (*ki)++ {
		(*incx) = (*incxs)[(*ki)-1]
		(*incy) = (*incys)[(*ki)-1]
		(*mx) = (ABS((*incx)))
		(*my) = (ABS((*incy)))
		//*
		for (*kn) = 1; (*kn) <= 4; (*kn)++ {
			(*n) = (*ns)[(*kn)-1]
			(*ksize) = (MIN(int(2), (*kn)))
			(*lenx) = (*lens)[(*kn)-1][(*mx)-1]
			(*leny) = (*lens)[(*kn)-1][(*my)-1]
			//*
			if (*icase) == 4 {
				//*              .. Drot ..
				for (*i) = 1; (*i) <= 7; (*i)++ {
					(*sx)[(*i)-1] = (*dx1)[(*i)-1]
					(*sy)[(*i)-1] = (*dy1)[(*i)-1]
					(*stx)[(*i)-1] = (*dt9x)[(*i)-1][(*kn)-1][(*ki)-1]
					(*sty)[(*i)-1] = (*dt9y)[(*i)-1][(*kn)-1][(*ki)-1]
				//Label20:
				}
				Drot(n, sx, incx, sy, incy, sc, ss)
				Stest(lenx, sx, stx, &((*ssize2)[0][(*ksize)-1]), sfac)
				Stest(leny, sy, sty, &((*ssize2)[0][(*ksize)-1]), sfac)
			} else {
				WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}(), *func() *[]byte{y := []byte(" shouldn't be here in check3"); return &y}())
				panic("")
			}
		//Label40:
		}
	//Label60:
	}
	//*
	(*mwpc)[0] = 1
	for (*i) = 2; (*i) <= 11; (*i)++ {
		(*mwpc)[(*i)-1] = 0
	//Label80:
	}
	(*mwps)[0] = 0
	for (*i) = 2; (*i) <= 6; (*i)++ {
		(*mwps)[(*i)-1] = 1
	//Label100:
	}
	for (*i) = 7; (*i) <= 11; (*i)++ {
		(*mwps)[(*i)-1] = -1
	//Label120:
	}
	(*mwpinx)[0] = 1
	(*mwpinx)[1] = 1
	(*mwpinx)[2] = 1
	(*mwpinx)[3] = -1
	(*mwpinx)[4] = 1
	(*mwpinx)[5] = -1
	(*mwpinx)[6] = 1
	(*mwpinx)[7] = 1
	(*mwpinx)[8] = -1
	(*mwpinx)[9] = 1
	(*mwpinx)[10] = -1
	(*mwpiny)[0] = 1
	(*mwpiny)[1] = 1
	(*mwpiny)[2] = -1
	(*mwpiny)[3] = -1
	(*mwpiny)[4] = 2
	(*mwpiny)[5] = 1
	(*mwpiny)[6] = 1
	(*mwpiny)[7] = -1
	(*mwpiny)[8] = -1
	(*mwpiny)[9] = 2
	(*mwpiny)[10] = 1
	for (*i) = 1; (*i) <= 11; (*i)++ {
		(*mwpn)[(*i)-1] = 5
	//Label140:
	}
	(*mwpn)[4] = 3
	(*mwpn)[9] = 3
	for (*i) = 1; (*i) <= 5; (*i)++ {
		(*mwpx)[(*i)-1] = (*i)
		(*mwpy)[(*i)-1] = (*i)
		(*mwptx)[0][(*i)-1] = (*i)
		(*mwpty)[0][(*i)-1] = (*i)
		(*mwptx)[1][(*i)-1] = (*i)
		(*mwpty)[1][(*i)-1] = -(*i)
		(*mwptx)[2][(*i)-1] = 6 - (*i)
		(*mwpty)[2][(*i)-1] = (*i) - 6
		(*mwptx)[3][(*i)-1] = (*i)
		(*mwpty)[3][(*i)-1] = -(*i)
		(*mwptx)[5][(*i)-1] = 6 - (*i)
		(*mwpty)[5][(*i)-1] = (*i) - 6
		(*mwptx)[6][(*i)-1] = -(*i)
		(*mwpty)[6][(*i)-1] = (*i)
		(*mwptx)[7][(*i)-1] = (*i) - 6
		(*mwpty)[7][(*i)-1] = 6 - (*i)
		(*mwptx)[8][(*i)-1] = -(*i)
		(*mwpty)[8][(*i)-1] = (*i)
		(*mwptx)[10][(*i)-1] = (*i) - 6
		(*mwpty)[10][(*i)-1] = 6 - (*i)
	//Label160:
	}
	(*mwptx)[4][0] = 1
	(*mwptx)[4][1] = 3
	(*mwptx)[4][2] = 5
	(*mwptx)[4][3] = 4
	(*mwptx)[4][4] = 5
	(*mwpty)[4][0] = -1
	(*mwpty)[4][1] = 2
	(*mwpty)[4][2] = -2
	(*mwpty)[4][3] = 4
	(*mwpty)[4][4] = -3
	(*mwptx)[9][0] = -1
	(*mwptx)[9][1] = -3
	(*mwptx)[9][2] = -5
	(*mwptx)[9][3] = 4
	(*mwptx)[9][4] = 5
	(*mwpty)[9][0] = 1
	(*mwpty)[9][1] = 2
	(*mwpty)[9][2] = 2
	(*mwpty)[9][3] = 4
	(*mwpty)[9][4] = 3
	for (*i) = 1; (*i) <= 11; (*i)++ {
		(*incx) = (*mwpinx)[(*i)-1]
		(*incy) = (*mwpiny)[(*i)-1]
		for (*k) = 1; (*k) <= 5; (*k)++ {
			(*copyx)[(*k)-1] = (*mwpx)[(*k)-1]
			(*copyy)[(*k)-1] = (*mwpy)[(*k)-1]
			(*mwpstx)[(*k)-1] = (*mwptx)[(*i)-1][(*k)-1]
			(*mwpsty)[(*k)-1] = (*mwpty)[(*i)-1][(*k)-1]
		//Label180:
		}
		Drot(&((*mwpn)[(*i)-1]), copyx, incx, copyy, incy, &((*mwpc)[(*i)-1]), &((*mwps)[(*i)-1]))
		Stest(func() *int{y := 5; return &y}(), copyx, mwpstx, mwpstx, sfac)
		Stest(func() *int{y := 5; return &y}(), copyy, mwpsty, mwpsty, sfac)
	//Label200:
	}
	return
}

func Stest(_len *int, scomp *[]float64, strue *[]float64, ssize *[]float64, sfac *float64) {
	nout := new(int)
	zero := new(float64)
	len := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sd := new(float64)
	i := new(int)
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     ********************************* Stest **************************
	//*
	//*     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
	//*     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
	//*     NEGLIGIBLE.
	//*
	//*     C. L. LAWSON, JPL, 1974 DEC 10
	//*
	//*     .. Parameters ..
	(*nout) = 6
	(*zero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. External Functions ..
	//*     .. Intrinsic Functions ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Executable Statements ..
	//*
	for (*i) = 1; (*i) <= (*len); (*i)++ {
		(*sd) = (*scomp)[(*i)-1] - (*strue)[(*i)-1]
		if (*ABS((*sfac) * (*sd)) <= ABS(&((*ssize)[(*i)-1]))) * (EPSILON((*zero))) {
			goto Label40
		}
		//*
		//*                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
		//*
		if !(*pass) {
			goto Label20
		}
		//*                             PRINT FAIL MESSAGE AND Header.
		(*pass) = false
		WRITE((*nout), *func() *[]byte{y := []byte("                                       fail\n"); return &y}())
		WRITE((*nout), *func() *[]byte{y := []byte("\n case  n incx incy  i                             comp(i)                             true(i)  difference     size(i)\n \n"); return &y}())
	Label20:
		;
		WRITE((*nout), *func() *[]byte{y := []byte(" %4d%3d%5d%3d%v%v%v%v\n"); return &y}(), (*icase), (*n), (*incx), (*incy), (*i), (*scomp)[(*i)-1], (*strue)[(*i)-1], (*sd), (*ssize)[(*i)-1])
	Label40:
	}
	return
	//*
}

func Testdsdot(scomp *float64, strue *float64, ssize *float64, sfac *float64) {
	nout := new(int)
	zero := new(float64)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	sd := new(float64)
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     ********************************* Stest **************************
	//*
	//*     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
	//*     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
	//*     NEGLIGIBLE.
	//*
	//*     C. L. LAWSON, JPL, 1974 DEC 10
	//*
	//*     .. Parameters ..
	(*nout) = 6
	(*zero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Intrinsic Functions ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Executable Statements ..
	//*
	(*sd) = (*scomp) - (*strue)
	if (*ABS((*sfac) * (*sd)) <= ABS(ssize)) * (EPSILON((*zero))) {
		goto Label40
	}
	//*
	//*                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
	//*
	if !(*pass) {
		goto Label20
	}
	//*                             PRINT FAIL MESSAGE AND Header.
	(*pass) = false
	WRITE((*nout), *func() *[]byte{y := []byte("                                       fail\n"); return &y}())
	WRITE((*nout), *func() *[]byte{y := []byte("\n case  n incx incy                            comp(i)                             true(i)  difference     size(i)\n \n"); return &y}())
Label20:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" %4d%3d%5d%3d%v%v%v%v\n"); return &y}(), (*icase), (*n), (*incx), (*incy), (*scomp), (*strue), (*sd), (*ssize))
Label40:
	;
	return
	//*
}

func Stest1(scomp1 *float64, strue1 *float64, ssize *[]float64, sfac *float64) {
	scomp := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	strue := func() *[]float64 {
		arr := make([]float64, 1)
		return &arr
	}()
	//*     ************************* Stest1 *****************************
	//*
	//*     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//*     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//*     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//*
	//*     C.L. LAWSON, JPL, 1978 DEC 6
	//*
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Arrays ..
	//*     .. External Subroutines ..
	//*     .. Executable Statements ..
	//*
	(*scomp)[0] = (*scomp1)
	(*strue)[0] = (*strue1)
	Stest(func() *int{y := 1; return &y}(), scomp, strue, ssize, sfac)
	//*
	return
}

func Sdiff(sa *float64, sb *float64) (sdiffReturn *float64) {
	sdiffreturn := new(float64)
	//*     ********************************* SDIFF **************************
	//*     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
	//*
	//*     .. Scalar Arguments ..
	//*     .. Executable Statements ..
	(*sdiff) = (*sa) - (*sb)
	return
}

func Itest1(icomp *int, itrue *int) {
	nout := new(int)
	icase := new(int)
	incx := new(int)
	incy := new(int)
	n := new(int)
	pass := new(bool)
	id := new(int)
	common.combla.pass = new(bool)
	common.combla.incy = new(int)
	common.combla.incx = new(int)
	common.combla.n = new(int)
	common.combla.icase = new(int)
	//*     ********************************* Itest1 *************************
	//*
	//*     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
	//*     EQUALITY.
	//*     C. L. LAWSON, JPL, 1974 DEC 10
	//*
	//*     .. Parameters ..
	(*nout) = 6
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Common blocks ..
	icase = common.combla.icase
	n = common.combla.n
	incx = common.combla.incx
	incy = common.combla.incy
	pass = common.combla.pass
	//*     .. Executable Statements ..
	//*
	if (*icomp) == (*itrue) {
		goto Label40
	}
	//*
	//*                            HERE ICOMP IS NOT EQUAL TO ITRUE.
	//*
	if !(*pass) {
		goto Label20
	}
	//*                             PRINT FAIL MESSAGE AND Header.
	(*pass) = false
	WRITE((*nout), *func() *[]byte{y := []byte("                                       fail\n"); return &y}())
	WRITE((*nout), *func() *[]byte{y := []byte("\n case  n incx incy                                comp                                true     difference\n \n"); return &y}())
Label20:
	;
	(*id) = (*icomp) - (*itrue)
	WRITE((*nout), *func() *[]byte{y := []byte(" %4d%3d%5d%36d%12d\n"); return &y}(), (*icase), (*n), (*incx), (*incy), (*icomp), (*itrue), (*id))
Label40:
	;
	return
	//*
}
