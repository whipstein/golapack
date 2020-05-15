package main

// Sblat1 ...
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       PROGRAM SBLAT1
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//    Test program for the REAL Level 1 BLAS.
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
// \ingroup single_blas_testing
//
//  =====================================================================
func Sblat1() {
	var sfac float32

	sfac = 9.765625e-4

	write(6, []byte(" Real BLAS Test Program Results\n \n"))
	for ic := 1; ic <= 13; ic++ {
		common.combla.icase = ic
		header()
		common.combla.pass = true
		common.combla.incx = 9999
		common.combla.incy = 9999
		switch {
		case common.combla.icase == 3 || common.combla.icase == 11:
			check0(sfac)
		case common.combla.icase == 7 || common.combla.icase == 8 || common.combla.icase == 9 || common.combla.icase == 10:
			check1(sfac)
		case common.combla.icase == 1 || common.combla.icase == 2 || common.combla.icase == 5 || common.combla.icase == 6 || common.combla.icase == 12 || common.combla.icase == 13:
			check2(sfac)
		case common.combla.icase == 4:
			check3(sfac)
		default:
			panic("Test case out of range!")
		}
		if common.combla.pass {
			write(6, []byte("                                    ----- PASS -----\n"))
		}
	}
}

func header() {
	l := []string{
		" SDOT ",
		"SAXPY ",
		"SROTG ",
		" SROT ",
		"SCOPY ",
		"SSWAP ",
		"SNRM2 ",
		"SASUM ",
		"SSCAL ",
		"ISAMAX",
		"SROTMG",
		"SROTM ",
		"SDSDOT",
	}

	write(6, []byte("\n Test of subprogram number%3d            %6s\n"), common.combla.icase, l[common.combla.icase-1])
	return
}

func check0(sfac float32) {
	var sa, sb, sc, ss float32
	da1 := []float32{0.3, 0.4, -0.3, -0.4, -0.3, 0, 0, 1}
	datrue := []float32{0.5, 0.5, 0.5, -0.5, -0.5, 0, 1, 1}
	db1 := []float32{0.4, 0.3, 0.4, 0.3, -0.4, 0, 1, 0}
	dbtrue := []float32{0, 0.6, 0, -0.6, 0, 0, 1, 0}
	dc1 := []float32{0.6, 0.8, -0.6, 0.8, 0.6, 1, 0, 1}
	ds1 := []float32{0.8, 0.6, 0.8, -0.6, 0.8, 0, 1, 0}
	//     INPUT FOR MODIFIED GIVENS
	dab := [][]float32{
		{0.1, 0.3, 1.2, 0.2},
		{0.7, 0.2, 0.6, 4.2},
		{0.0, 0.0, 0.0, 0.0},
		{4.0, -1.0, 2.0, 4.0},
		{6e-10, 2.e-2, 1e5, 10.0},
		{4e10, 2e-2, 1e-5, 10.0},
		{2e-10, 4e-2, 1e5, 10.0},
		{2e10, 4e-2, 1e-5, 10.0},
		{4.0, -2.0, 8.0, 4.0},
	}
	//    TRUE RESULTS FOR MODIFIED GIVENS
	dtrue := [][]float32{
		{0.0, 0.0, 1.3, 0.2, 0.0, 0.0, 0.0, 0.5, 0.0},
		{0.0, 0.0, 4.5, 4.2, 1.0, 0.5, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 4.0, -1.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 15e-3, 0.0, 10.0, -1.0, 0.0, -1e-4, 0.0, 1.0},
		{0.0, 0.0, 6144e-5, 10.0, -1.0, 4096, -1e6, 0.0, 1.0},
		{0.0, 0.0, 15.0, 10.0, -1.0, 5e-5, 0.0, 1.0, 0.0},
		{0.0, 0.0, 15.0, 10.0, -1.0, 5e5, -4096, 1.0, 4096e-6},
		{0.0, 0.0, 7.0, 4.0, 0.0, 0.0, -0.5, -0.25, 0.0},
	}
	//                   4096 = 2 ** 12
	d12 := float32(4096)
	dtemp := []float32{0, 0, 0, 0, 0, 0, 0, 0, 0}

	dtrue[0][0] = 12 / 130
	dtrue[0][1] = 36 / 130
	dtrue[0][6] = -1 / 6
	dtrue[1][0] = 14 / 75
	dtrue[1][1] = 49 / 75
	dtrue[1][8] = 1 / 7
	dtrue[4][0] = 45e-11 * (d12 * d12)
	dtrue[4][2] = 4e5 / (3. * d12)
	dtrue[4][5] = 1 / d12
	dtrue[4][7] = 1e4 / (3. * d12)
	dtrue[5][0] = 4e10 / (1.5 * d12 * d12)
	dtrue[5][1] = 2e-2 / 1.5
	dtrue[5][7] = 5e-7 * d12
	dtrue[6][0] = 4 / 150
	dtrue[6][1] = (2e-10 / 1.5) * (d12 * d12)
	dtrue[6][6] = -dtrue[4][5]
	dtrue[6][8] = 1e4 / d12
	dtrue[7][0] = dtrue[6][0]
	dtrue[7][1] = 2e10 / (1.5 * d12 * d12)
	dtrue[8][0] = 32 / 7
	dtrue[8][1] = -16 / 7

	//     Compute true values which cannot be prestored
	//     in decimal notation
	dbtrue[0] = 1.0 / 0.6
	dbtrue[2] = -1.0 / 0.6
	dbtrue[4] = 1.0 / 0.6
	//
	for k := 0; k < 8; k++ {
		//        .. Set N=K for identification in output if any ..
		common.combla.n = k
		if common.combla.icase == 3 {
			//           .. SROTG ..
			sa = da1[k]
			sb = db1[k]
			Srotg(&sa, &sb, &sc, &ss)
			stest1(sa, datrue[k], datrue[k], sfac, []byte("sa mismatch"))
			stest1(sb, dbtrue[k], dbtrue[k], sfac, []byte("sb mismatch"))
			stest1(sc, dc1[k], dc1[k], sfac, []byte("sc mismatch"))
			stest1(ss, ds1[k], ds1[k], sfac, []byte("sd mismatch"))
		} else if common.combla.icase == 11 {
			//           .. SROTMG ..
			dtempx := make([]float32, 5)
			dtrue1 := make([]float32, 9)
			for i := 0; i < 4; i++ {
				dtemp[i] = dab[k][i]
				dtempx[i] = 0.0
			}
			dtempx[5] = 0.0
			for i := 0; i < 5; i++ {
				dtempx[i] = dtemp[i+4]
			}
			Srotmg(dtemp[0], dtemp[1], dtemp[2], dtemp[3], &dtempx)
			for i := 0; i < 9; i++ {
				if i <= 5 {
					dtemp[i+4] = dtempx[i]
				}
				dtrue1[i] = dtrue[k][i]
			}
			stest(9, &dtemp, &dtrue1, &dtrue1, sfac, []byte("dtemp mismatch"))
		} else {
			write(6, []byte(" %v\n"), []byte(" Shouldn't be here in CHECK0"))
			panic("")
		}
	}
	return
}

func check1(sfac float32) {
	var len int
	dtrue1 := []float32{0.0, 0.3, 0.5, 0.7, 0.6}
	dtrue3 := []float32{0.0, 0.3, 0.7, 1.1, 1.0}
	dtrue5 := [][][]float32{
		{{0.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, {-0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0}, {0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}, {0.20, -0.60, 0.30, 5.0, 5.0, 5.0, 5.0, 5.0}, {0.03, -0.09, 0.15, -0.03, 6.0, 6.0, 6.0, 6.0}},
		{{0.10, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0}, {0.09, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0}, {0.09, 2.0, -0.12, 2.0, 2.0, 2.0, 2.0, 2.0}, {0.06, 3.0, -0.18, 5.0, 0.09, 2.0, 2.0, 2.0}, {0.03, 4.0, -0.09, 6.0, -0.15, 7.0, -0.03, 3.0}},
	}
	dv := [][][]float32{
		{{0.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, {0.3, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0}, {0.3, -0.4, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}, {0.2, -0.6, 0.3, 5.0, 5.0, 5.0, 5.0, 5.0}, {0.1, -0.3, 0.5, -0.1, 6.0, 6.0, 6.0, 6.0}},
		{{0.1, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0}, {0.3, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0}, {0.3, 2.0, -0.4, 2.0, 2.0, 2.0, 2.0, 2.0}, {0.2, 3.0, -0.6, 5.0, 0.3, 2.0, 2.0, 2.0}, {0.1, 4.0, -0.3, 6.0, -0.5, 7.0, -0.1, 3.0}},
	}
	sa := []float32{0.3, -1.0, 0.0, 1.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3}
	strue := make([]float32, 8)
	sx := make([]float32, 8)
	itrue2 := []int{0, 1, 2, 2, 3}

	for common.combla.incx = 1; common.combla.incx <= 2; common.combla.incx++ {
		for np1 := 0; np1 < 5; np1++ {
			common.combla.n = np1
			len = 2 * max(common.combla.n, 0)
			//           .. Set vector arguments ..
			for i := 0; i < len; i++ {
				sx[i] = dv[common.combla.incx-1][np1][i]
			}
			//
			if common.combla.icase == 7 {
				//              .. SNRM2 ..
				stest1(Snrm2(common.combla.n, &sx, common.combla.incx), dtrue1[np1], dtrue1[np1], sfac, []byte("Snrm2 mismatch"))
			} else if common.combla.icase == 8 {
				//              .. SASUM ..
				stest1(Sasum(common.combla.n, &sx, common.combla.incx), dtrue3[np1], dtrue3[np1], sfac, []byte("Sasum mismatch"))
			} else if common.combla.icase == 9 {
				//              .. SSCAL ..
				Sscal(common.combla.n, sa[(common.combla.incx-1)*5+np1], &sx, common.combla.incx)
				for i := 0; i < len; i++ {
					strue[i] = dtrue5[common.combla.incx-1][np1][i]
				}
				stest(len, &sx, &strue, &strue, sfac, []byte("sx mismatch"))
			} else if common.combla.icase == 10 {
				//              .. ISAMAX ..
				itest1(Isamax(common.combla.n, &sx, common.combla.incx), itrue2[np1])
			} else {
				write(6, []byte(" %v\n"), []byte(" Shouldn't be here in CHECK1"))
				panic("")
			}
		}
	}
	return
}

func check2(sfac float32) {
	var kni, ksize, lenx, leny, mx, my int
	dt10x := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.3, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0}, {0.8, 0.1, -0.6, 0.8, 0.3, -0.3, 0.5}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0}, {0.7, 0.1, 0.3, 0.8, -0.9, -0.3, 0.5}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.3, -0.6, 0.8, 0.0, 0.0, 0.0}},
	}
	dt10y := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.5, -0.9, 0.6, 0.0, 0.0, 0.0, 0.0}, {-0.4, -0.9, 0.9, 0.7, -0.5, 0.2, 0.6}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.5, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.4, 0.9, -0.5, 0.6, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, -0.9, 0.1, 0.0, 0.0, 0.0, 0.0}, {0.6, -0.9, 0.1, 0.7, -0.5, 0.2, 0.8}},
	}
	dt7 := [][]float32{
		{0.0, 0.30, 0.21, 0.62},
		{0.0, 0.30, -0.07, 0.85},
		{0.0, 0.30, -0.79, -0.74},
		{0.0, 0.30, 0.33, 1.27},
	}
	dt8 := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, -0.87, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, -0.87, 0.15, 0.94, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.35, -0.9, 0.48, 0.0, 0.0, 0.0, 0.0}, {0.38, -0.9, 0.57, 0.7, -0.75, 0.2, 0.98}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.35, -0.72, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.38, -0.63, 0.15, 0.88, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.68, -0.9, 0.33, 0.0, 0.0, 0.0, 0.0}, {0.68, -0.9, 0.33, 0.7, -0.75, 0.2, 1.04}},
	}
	dx1 := []float32{0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4}
	dy1 := []float32{0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8}
	ssize1 := []float32{0, 0.3, 1.6, 3.2}
	ssize2 := [][]float32{
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17},
	}
	ssize3 := []float32{.1, .4, 1.7, 3.3}
	ssize := make([]float32, 7)
	stx := make([]float32, 7)
	sty := make([]float32, 7)
	sx := make([]float32, 7)
	sy := make([]float32, 7)
	//                         FOR DROTM
	dpar := [][]float32{
		{-2.0, 0.0, 0.0, 0.0, 0.0},
		{-1.0, 2.0, -3.0, -4.0, 5.0},
		{0.0, 0.0, 2.0, -3.0, 0.0},
		{1.0, 5.0, 2.0, 0.0, -4.0},
	}
	//                        TRUE X RESULTS F0R ROTATIONS DROTM
	dt19xa := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, 3.8, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 2.8, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0}, {-0.8, 3.8, -2.2, -1.2, 0.0, 0.0, 0.0}, {-0.9, 2.8, -1.4, -1.3, 0.0, 0.0, 0.0}, {3.5, -0.4, -2.2, 4.7, 0.0, 0.0, 0.0}},
	}
	dt19xb := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.1, -3.0, 0.0, 0.0, 0.0, 0.0}, {-0.3, 0.1, -2.0, 0.0, 0.0, 0.0, 0.0}, {3.3, 0.1, -2.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4}, {-2.0, 0.1, 1.4, 0.8, 0.6, -0.3, -2.8}, {-1.8, 0.1, 1.3, 0.8, 0.0, -0.3, -1.9}, {3.8, 0.1, -3.1, 0.8, 4.8, -0.3, -1.5}},
	}
	dt19xc := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.0, 0.0, 0.0, 0.0}, {4.8, 0.1, -3.0, 0.0, 0.0, 0.0, 0.0}, {3.3, 0.1, -2.0, 0.0, 0.0, 0.0, 0.0}, {2.1, 0.1, -2.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4}, {-1.6, 0.1, -2.2, 0.8, 5.4, -0.3, -2.8}, {-1.5, 0.1, -1.4, 0.8, 3.6, -0.3, -1.9}, {3.7, 0.1, -2.2, 0.8, 3.6, -0.3, -1.5}},
	}
	dt19xd := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.8, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.9, -0.8, 0.0, 0.0, 0.0, 0.0, 0.0}, {3.5, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.6, 0.1, -0.5, 0.8, 0.0, 0.0, 0.0}, {-0.8, -1.0, 1.4, -1.6, 0.0, 0.0, 0.0}, {-0.9, -0.8, 1.3, -1.6, 0.0, 0.0, 0.0}, {3.5, 0.8, -3.1, 4.8, 0.0, 0.0, 0.0}},
	}
	//                        TRUE Y RESULTS FOR ROTATIONS DROTM
	dt19ya := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {1.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, -4.8, 0.0, 0.0, 0.0, 0.0, 0.0}, {1.7, -0.7, 0.0, 0.0, 0.0, 0.0, 0.0}, {-2.6, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0}, {0.7, -4.8, 3.0, 1.1, 0.0, 0.0, 0.0}, {1.7, -0.7, -0.7, 2.3, 0.0, 0.0, 0.0}, {-2.6, 3.5, -0.7, -3.6, 0.0, 0.0, 0.0}},
	}
	dt19yb := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {1.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.0, 0.0, 0.0, 0.0}, {4.0, -0.9, -0.3, 0.0, 0.0, 0.0, 0.0}, {-0.5, -0.9, 1.5, 0.0, 0.0, 0.0, 0.0}, {-1.5, -0.9, -1.8, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8}, {3.7, -0.9, -1.2, 0.7, -1.5, 0.2, 2.2}, {-0.3, -0.9, 2.1, 0.7, -1.6, 0.2, 2.0}, {-1.6, -0.9, -2.1, 0.7, 2.9, 0.2, -3.8}},
	}
	dt19yc := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {1.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0}, {4.0, -6.3, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0}, {-1.5, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.7, 0.0, 0.0, 0.0}, {3.7, -7.2, 3.0, 1.7, 0.0, 0.0, 0.0}, {-0.3, 0.9, -0.7, 1.9, 0.0, 0.0, 0.0}, {-1.6, 2.7, -0.7, -3.4, 0.0, 0.0, 0.0}},
	}
	dt19yd := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {1.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-2.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.0, 0.0, 0.0, 0.0}, {0.7, -0.9, 1.2, 0.0, 0.0, 0.0, 0.0}, {1.7, -0.9, 0.5, 0.0, 0.0, 0.0, 0.0}, {-2.6, -0.9, -1.3, 0.0, 0.0, 0.0, 0.0}},
		{{0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8}, {0.7, -0.9, 1.2, 0.7, -1.5, 0.2, 1.6}, {1.7, -0.9, 0.5, 0.7, -1.6, 0.2, 2.4}, {-2.6, -0.9, -1.3, 0.7, 2.9, 0.2, -4.0}},
	}
	var dt19x [16][4][7]float32
	var dt19y [16][4][7]float32
	lendt19x0 := len(dt19xa)
	lendt19x1 := len(dt19xa[0])
	lendt19x2 := len(dt19xa[0][0])
	for k := 0; k < lendt19x0; k++ {
		for j := 0; j < lendt19x1; j++ {
			for i := 0; i < lendt19x2; i++ {
				dt19x[k][j][i] = dt19xa[k][j][i]
				dt19x[k+lendt19x0][j][i] = dt19xb[k][j][i]
				dt19x[k+2*lendt19x0][j][i] = dt19xc[k][j][i]
				dt19x[k+3*lendt19x0][j][i] = dt19xd[k][j][i]
				dt19y[k][j][i] = dt19ya[k][j][i]
				dt19y[k+lendt19x0][j][i] = dt19yb[k][j][i]
				dt19y[k+2*lendt19x0][j][i] = dt19yc[k][j][i]
				dt19y[k+3*lendt19x0][j][i] = dt19yd[k][j][i]
			}
		}
	}

	dtemp := make([]float32, 5)
	st7b := [][]float32{
		{.1, 0.4, 0.31, 0.72},
		{0.1, 0.4, 0.03, 0.95},
		{0.1, 0.4, -0.69, -0.64},
		{0.1, 0.4, 0.43, 1.37},
	}
	incxs := []int{1, 2, -2, -1}
	incys := []int{1, -2, 1, -2}
	lens := [][]int{
		{1, 1},
		{1, 1},
		{2, 3},
		{4, 7},
	}
	ns := []int{0, 1, 2, 4}
	sa := float32(0.3)

	for ki := 0; ki < 4; ki++ {
		common.combla.incx = incxs[ki]
		common.combla.incy = incys[ki]
		mx = absint(common.combla.incx) - 1
		my = absint(common.combla.incy) - 1

		for kn := 0; kn < 4; kn++ {
			common.combla.n = ns[kn]
			ksize = min(1, kn)
			lenx = lens[kn][mx]
			leny = lens[kn][my]

			//           .. Initialize all argument arrays ..
			for i := 0; i < 7; i++ {
				sx[i] = dx1[i]
				sy[i] = dy1[i]
			}

			switch common.combla.icase {
			case 1:
				//              .. SDOT ..
				stest1(Sdot(common.combla.n, &sx, common.combla.incx, &sy, common.combla.incy), dt7[ki][kn], ssize1[kn], sfac, []byte("Sdot mismatch"))
			case 2:
				//              .. SAXPY ..
				Saxpy(common.combla.n, sa, &sx, common.combla.incx, &sy, common.combla.incy)
				for j := 0; j < leny; j++ {
					sty[j] = dt8[ki][kn][j]
					//Label40:
				}
				ssize2x := ssize2[ksize:][0]
				stest(leny, &sy, &sty, &ssize2x, sfac, []byte("sy mismatch"))
			case 5:
				//              .. SCOPY ..
				for i := 0; i < 7; i++ {
					sty[i] = dt10y[ki][kn][i]
					//Label60:
				}
				Scopy(common.combla.n, &sx, common.combla.incx, &sy, common.combla.incy)
				stest(leny, &sy, &sty, &ssize2[:][0], 1, []byte("sy mismatch"))
			case 6:
				//              .. SSWAP ..
				Sswap(common.combla.n, &sx, common.combla.incx, &sy, common.combla.incy)
				for i := 0; i < 7; i++ {
					stx[i] = dt10x[ki][kn][i]
					sty[i] = dt10y[ki][kn][i]
					//Label80:
				}
				stest(lenx, &sx, &stx, &ssize2[:][0], 1, []byte("sx mismatch"))
				stest(leny, &sy, &sty, &ssize2[:][0], 1, []byte("sy mismatch"))
			case 12:
				//              .. SROTM ..
				kni = kn + 4*ki
				for kpar := 0; kpar < 4; kpar++ {
					for i := 0; i < 7; i++ {
						sx[i] = dx1[i]
						sy[i] = dy1[i]
						stx[i] = dt19x[kni][kpar][i]
						sty[i] = dt19y[kni][kpar][i]
					}
					//
					for i := 0; i < 5; i++ {
						dtemp[i] = dpar[kpar][i]
					}
					//
					for i := 0; i < lenx; i++ {
						ssize[i] = stx[i]
					}
					//                   SEE REMARK ABOVE ABOUT DT11X(1,2,7)
					//                       AND DT11X(5,3,8).
					if (kpar == 2) && (kni == 7) {
						ssize[0] = 2.4
					}
					if (kpar == 3) && (kni == 8) {
						ssize[4] = 1.8
					}
					//
					Srotm(common.combla.n, &sx, common.combla.incx, &sy, common.combla.incy, &dtemp)
					stest(lenx, &sx, &stx, &ssize, sfac, []byte("sx mismatch"))
					stest(leny, &sy, &sty, &sty, sfac, []byte("sy mismatch"))
				}
			case 13:
				//              .. SDSROT ..
				stest1(Sdsdot(common.combla.n, 0.1, &sx, common.combla.incx, &sy, common.combla.incy), st7b[ki][kn], ssize3[kn], sfac, []byte("Sdsdot mismatch"))
			default:
				write(6, []byte(" %v\n"), []byte(" Shouldn't be here in CHECK2"))
				panic("")
			}
		}
	}
	return
}

func check3(sfac float32) {
	var ksize, lenx, leny, mx, my int
	copyx := make([]float32, 5)
	copyy := make([]float32, 5)
	dt9x := [][][]float32{
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, -0.46, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, -0.46, -0.22, 1.06, 0.0, 0.0, 0.0}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.66, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0}, {0.96, 0.1, -0.76, 0.8, 0.90, -0.3, -0.02}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {-0.06, 0.1, -0.1, 0.0, 0.0, 0.0, 0.0}, {0.90, 0.1, -0.22, 0.8, 0.18, -0.3, -0.02}},
		{{0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.26, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.78, 0.26, -0.76, 1.12, 0.0, 0.0, 0.0}},
	}
	dt9y := [][][]float32{
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, -0.78, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, -0.78, 0.54, 0.08, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, -0.9, -0.12, 0.0, 0.0, 0.0, 0.0}, {0.64, -0.9, -0.30, 0.7, -0.18, 0.2, 0.28}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.7, -1.08, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.64, -1.26, 0.54, 0.20, 0.0, 0.0, 0.0}},
		{{0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.04, -0.9, 0.18, 0.0, 0.0, 0.0, 0.0}, {0.04, -0.9, 0.18, 0.7, -0.18, 0.2, 0.16}},
	}
	dx1 := []float32{0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4}
	dy1 := []float32{0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8}
	mwpc := make([]float32, 11)
	mwps := make([]float32, 11)
	mwpstx := make([]float32, 5)
	mwpsty := make([]float32, 5)
	var mwptx [5][11]float32
	var mwpty [5][11]float32
	mwpx := make([]float32, 5)
	mwpy := make([]float32, 5)
	ssize2 := [][]float32{
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17},
	}
	sx := make([]float32, 7)
	sy := make([]float32, 7)
	stx := make([]float32, 7)
	sty := make([]float32, 7)
	incxs := []int{1, 2, -2, -1}
	incys := []int{1, -2, 1, -2}
	lens := [][]int{
		{1, 1, 2, 4},
		{1, 1, 3, 7},
	}
	mwpinx := make([]int, 11)
	mwpiny := make([]int, 11)
	mwpn := make([]int, 11)
	ns := []int{0, 1, 2, 4}
	sc := float32(0.8)
	ss := float32(0.6)

	for ki := 0; ki < 4; ki++ {
		common.combla.incx = incxs[ki]
		common.combla.incy = incys[ki]
		mx = absint(common.combla.incx) - 1
		my = absint(common.combla.incy) - 1
		//
		for kn := 0; kn < 4; kn++ {
			common.combla.n = ns[kn]
			ksize = min(1, kn)
			lenx = lens[mx][kn]
			leny = lens[my][kn]
			//
			if common.combla.icase == 4 {
				//              .. SROT ..
				for i := 0; i < 7; i++ {
					sx[i] = dx1[i]
					sy[i] = dy1[i]
					stx[i] = dt9x[ki][kn][i]
					sty[i] = dt9y[ki][kn][i]
				}
				Srot(common.combla.n, &sx, common.combla.incx, &sy, common.combla.incy, sc, ss)
				ssize2x := ssize2[ksize:][0]
				stest(lenx, &sx, &stx, &ssize2x, sfac, []byte("sx mismatch"))
				stest(leny, &sy, &sty, &ssize2x, sfac, []byte("sy mismatch"))
			} else {
				write(6, []byte(" %v\n"), []byte(" Shouldn't be here in CHECK3"))
				panic("")
			}
		}
	}
	//
	mwpc[0] = 1
	for i := 1; i < 11; i++ {
		mwpc[i] = 0
		//Label80:
	}
	mwps[0] = 0
	for i := 1; i < 6; i++ {
		mwps[i] = 1
		//Label100:
	}
	for i := 6; i < 11; i++ {
		mwps[i] = -1
		//Label120:
	}
	mwpinx[0] = 1
	mwpinx[1] = 1
	mwpinx[2] = 1
	mwpinx[3] = -1
	mwpinx[4] = 1
	mwpinx[5] = -1
	mwpinx[6] = 1
	mwpinx[7] = 1
	mwpinx[8] = -1
	mwpinx[9] = 1
	mwpinx[10] = -1
	mwpiny[0] = 1
	mwpiny[1] = 1
	mwpiny[2] = -1
	mwpiny[3] = -1
	mwpiny[4] = 2
	mwpiny[5] = 1
	mwpiny[6] = 1
	mwpiny[7] = -1
	mwpiny[8] = -1
	mwpiny[9] = 2
	mwpiny[10] = 1
	for i := 0; i < 11; i++ {
		mwpn[i] = 5
		//Label140:
	}
	mwpn[4] = 3
	mwpn[9] = 3
	for i := 0; i < 5; i++ {
		ix := float32(i + 1)
		mwpx[i] = ix
		mwpy[i] = ix
		mwptx[i][0] = ix
		mwpty[i][0] = ix
		mwptx[i][1] = ix
		mwpty[i][1] = -ix
		mwptx[i][2] = 6 - ix
		mwpty[i][2] = ix - 6
		mwptx[i][3] = ix
		mwpty[i][3] = -ix
		mwptx[i][5] = 6 - ix
		mwpty[i][5] = ix - 6
		mwptx[i][6] = -ix
		mwpty[i][6] = ix
		mwptx[i][7] = ix - 6
		mwpty[i][7] = 6 - ix
		mwptx[i][8] = -ix
		mwpty[i][8] = ix
		mwptx[i][10] = ix - 6
		mwpty[i][10] = 6 - ix
		//Label160:
	}
	mwptx[0][4] = 1
	mwptx[1][4] = 3
	mwptx[2][4] = 5
	mwptx[3][4] = 4
	mwptx[4][4] = 5
	mwpty[0][4] = -1
	mwpty[1][4] = 2
	mwpty[2][4] = -2
	mwpty[3][4] = 4
	mwpty[4][4] = -3
	mwptx[0][9] = -1
	mwptx[1][9] = -3
	mwptx[2][9] = -5
	mwptx[3][9] = 4
	mwptx[4][9] = 5
	mwpty[0][9] = 1
	mwpty[1][9] = 2
	mwpty[2][9] = 2
	mwpty[3][9] = 4
	mwpty[4][9] = 3
	for i := 0; i < 11; i++ {
		common.combla.incx = mwpinx[i]
		common.combla.incy = mwpiny[i]
		for k := 0; k < 5; k++ {
			copyx[k] = mwpx[k]
			copyy[k] = mwpy[k]
			mwpstx[k] = mwptx[k][i]
			mwpsty[k] = mwpty[k][i]
			//Label180:
		}
		mwpcx := mwpc[i]
		mwpsx := mwps[i]
		Srot(mwpn[i], &copyx, common.combla.incx, &copyy, common.combla.incy, mwpcx, mwpsx)
		stest(5, &copyx, &mwpstx, &mwpstx, sfac, []byte("copyx mismatch"))
		stest(5, &copyy, &mwpsty, &mwpsty, sfac, []byte("copyy mismatch"))
		//Label200:
	}
	return
}

func stest(len int, scomp *[]float32, strue *[]float32, ssize *[]float32, sfac float32, info []byte) {
	var sd float32

	for i := 0; i < len; i++ {
		sd = (*scomp)[i] - (*strue)[i]
		if absf32(sfac*sd) <= absf32((*ssize)[i])*epsilonf32() {
			// if absf32(sd) <= absf32((*ssize)[i]) {
			continue
		}
		if !(common.combla.pass) {
			write(6, []byte("   %2d  %2d  %4d %4d  %2d   % .8f   % .8f   % .8f   % .8f    %v\n"), common.combla.icase, common.combla.n, common.combla.incx, common.combla.incy, i, (*scomp)[i], (*strue)[i], sd, (*ssize)[i], info)
			continue
		}
		common.combla.pass = false
		write(6, []byte("                                    ***** FAIL *****\n"))
		write(6, []byte("\n CASE   N  INCX INCY   I       COMP(I)       TRUE(I)    DIFFERENCE       SIZE(I)    INFO\n\n"))
	}
	return
}

func stest1(scomp1 float32, strue1 float32, ssize1 float32, sfac float32, info1 []byte) {
	//     ************************* STEST1 *****************************
	//
	//     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
	//     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
	//     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
	//
	//     C.L. LAWSON, JPL, 1978 DEC 6
	//
	//     .. Scalar Arguments ..
	//     .. Array Arguments ..
	//     .. Local Arrays ..
	//     .. External Subroutines ..
	//     .. Executable Statements ..
	//
	scomp := []float32{scomp1}
	ssize := []float32{ssize1}
	strue := []float32{strue1}
	stest(1, &scomp, &strue, &ssize, sfac, info1)
	//
	return
}

func sdiff(sa float32, sb float32) float32 {
	//     ********************************* SDIFF **************************
	//     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
	//
	//     .. Scalar Arguments ..
	//     .. Executable Statements ..
	return (sa - sb)
}

func itest1(icomp int, itrue int) {
	var id int
	//     ********************************* ITEST1 *************************
	//
	//     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
	//     EQUALITY.
	//     C. L. LAWSON, JPL, 1974 DEC 10
	//
	//     .. Parameters ..
	//     .. Scalar Arguments ..
	//     .. Scalars in Common ..
	//     .. Local Scalars ..
	//     .. Common blocks ..
	//     .. Executable Statements ..
	//
	if icomp == itrue {
		goto Label40
	}
	//
	//                            HERE ICOMP IS NOT EQUAL TO ITRUE.
	//
	if !common.combla.pass {
		goto Label20
	}
	//                             PRINT FAIL MESSAGE AND HEADER.
	common.combla.pass = false
	write(6, []byte("                                       FAIL\n"))
	write(6, []byte("\n CASE  N INCX INCY                                COMP                                TRUE     DIFFERENCE\n \n"))
Label20:
	;
	id = icomp - itrue
	write(6, []byte(" %4d%3d%5d%36d%12d\n"), common.combla.icase, common.combla.n, common.combla.incx, common.combla.incy, icomp, itrue, id)
Label40:
	;
	return
	//
}
