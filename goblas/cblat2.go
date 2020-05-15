package goblas

import "testing"

// \brief \b Cblat2
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       PROGRAM Cblat2
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Test program for the COMPLEX          Level 2 Blas.
//
// The program must be driven by a short data file. The first 18 records
// of the file are read using list-directed input, the last 17 records
// are read using the format ( A6, L2). An annotated example of a data
// file can be obtained by deleting the first 3 characters from the
// following 35 lines:
// 'cblat2.out'      NAME OF SUMMARY OUTPUT FILE
// 6                 UNIT NUMBER OF SUMMARY FILE
// 'CBLA2T.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
// -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
// F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
// F        LOGICAL FLAG, T TO STOP ON FAILURES.
// T        LOGICAL FLAG, T TO TEST ERROR EXITS.
// 16.0     THRESHOLD VALUE OF TEST RATIO
// 6                 NUMBER OF VALUES OF N
// 0 1 2 3 5 9       VALUES OF N
// 4                 NUMBER OF VALUES OF K
// 0 1 2 4           VALUES OF K
// 4                 NUMBER OF VALUES OF incx AND incy
// 1 2 -1 -2         VALUES OF incx AND incy
// 3                 NUMBER OF VALUES OF ALPHA
// (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
// 3                 NUMBER OF VALUES OF BETA
// (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
// Cgemv  T PUT F FOR NO TEST. SAME COLUMNS.
// Cgbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Chemv  T PUT F FOR NO TEST. SAME COLUMNS.
// Chbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Chpmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctrmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctbmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctpmv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctrsv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctbsv  T PUT F FOR NO TEST. SAME COLUMNS.
// Ctpsv  T PUT F FOR NO TEST. SAME COLUMNS.
// Cgerc  T PUT F FOR NO TEST. SAME COLUMNS.
// Cgeru  T PUT F FOR NO TEST. SAME COLUMNS.
// Cher   T PUT F FOR NO TEST. SAME COLUMNS.
// Chpr   T PUT F FOR NO TEST. SAME COLUMNS.
// Cher2  T PUT F FOR NO TEST. SAME COLUMNS.
// Chpr2  T PUT F FOR NO TEST. SAME COLUMNS.
//
// Further Details
// ===============
//
//    See:
//
//       Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
//       An  extended  set of Fortran  Basic Linear Algebra Subprograms.
//
//       Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
//       and  Computer Science  Division,  Argonne  National Laboratory,
//       9700 South Cass Avenue, Argonne, Illinois 60439, US.
//
//       Or
//
//       NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
//       Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
//       OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
//       Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.
//
//
// -- Written on 10-August-1987.
//    Richard Hanson, Sandia National Labs.
//    Jeremy Du Croz, NAG Central Office.
//
//    10-9-00:  Change STATUS='NEW' to 'UNKNOWN' so that the testers
//              can be run multiple times without deleting generated
//              output files (susan)
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
// \ingroup complex_blas_testing
//
//  =====================================================================
func TestCblat2(t *testing.T) {
	nin := new(int)
	nsubs := new(int)
	zero := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	nmax := new(int)
	incmax := new(int)
	ninmax := new(int)
	nidmax := new(int)
	nkbmax := new(int)
	nalmax := new(int)
	nbemax := new(int)
	eps := new(float64)
	err := new(float64)
	thresh := new(float64)
	i := new(int)
	isnum := new(int)
	j := new(int)
	n := new(int)
	nalf := new(int)
	nbet := new(int)
	nidim := new(int)
	ninc := new(int)
	nkb := new(int)
	nout := new(int)
	ntra := new(int)
	fatal := new(bool)
	ltestt := new(bool)
	rewi := new(bool)
	same := new(bool)
	sfatal := new(bool)
	trace := new(bool)
	tsterr := new(bool)
	trans := new(byte)
	snamet := func() *[]byte {
		arr := make([]byte, 6)
		return &arr
	}()
	snaps := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	summry := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	a := func() *[][]complex128 {
		arr := make([][]complex128, -1)
		for u := 0; u < -1; u++ {
			arr[u] = make([]complex128, -1)
		}
		return &arr
	}()
	aa := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	alf := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	as := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	bet := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	x := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	xs := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	xx := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	y := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	ys := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	yt := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	yy := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	z := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	g := func() *[]float64 {
		arr := make([]float64, -1)
		return &arr
	}()
	idim := func() *[]int {
		arr := make([]int, -1)
		return &arr
	}()
	inc := func() *[]int {
		arr := make([]int, -1)
		return &arr
	}()
	kb := func() *[]int {
		arr := make([]int, -1)
		return &arr
	}()
	ltest := func() *[]bool {
		arr := make([]bool, -1)
		return &arr
	}()
	snames := func() *[][]byte {
		arr := make([][]byte, 6)
		for u := 0; u < 6; u++ {
			arr[u] = make([]byte, -1)
		}
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 6)
		return &arr
	}()
	common.infoc.lerr = new(float64)
	common.infoc.ok = new(float64)
	common.infoc.noutc = new(float64)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new(int)
	//*
	//*  -- Reference BLAS test routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     April 2012
	//*
	//*  =====================================================================
	//*
	//*     .. Parameters ..
	(*nin) = 5
	(*nsubs) = 17
	(*zero) = (0.0 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	(*nmax) = 65
	(*incmax) = 2
	(*ninmax) = 7
	(*nidmax) = 9
	(*nkbmax) = 7
	(*nalmax) = 7
	(*nbemax) = 7
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//*     .. Data statements ..
	//F4GO: NOT IMPLEMENTED :data snames / "cgemv ", "cgbmv ", "chemv ", "chbmv ", "chpmv ", "ctrmv ", "ctbmv ", "ctpmv ", "ctrsv ", "ctbsv ", "ctpsv ", "cgerc ", "cgeru ", "cher  ", "chpr  ", "cher2 ", "chpr2 " /
	//*     .. Executable Statements ..
	//*
	//*     Read name and unit number for summary output file and open file.
	//*
	READ(nin, []byte(" %v\n"), summry)
	READ(nin, []byte(" %v\n"), nout)
	OPEN(nout, SUMMRY)
	(*noutc) = (*nout)
	//*
	//*     Read name and unit number for snapshot output file and open file.
	//*
	READ(nin, []byte(" %v\n"), snaps)
	READ(nin, []byte(" %v\n"), ntra)
	(*trace) = (*ntra) >= 0
	if *trace {
		OPEN(ntra, snaps)
	}
	//*     Read the flag that directs rewinding of the snapshot file.
	READ(nin, []byte(" %v\n"), rewi)
	(*rewi) = (*rewi) && (*trace)
	//*     Read the flag that directs stopping on any failure.
	READ(nin, []byte(" %v\n"), sfatal)
	//*     Read the flag that indicates whether error exits are to be tested.
	READ(nin, []byte(" %v\n"), tsterr)
	//*     Read the threshold value of the test ratio
	READ(nin, []byte(" %v\n"), thresh)
	//*
	//*     Read and check the parameter values for the tests.
	//*
	//*     Values of N
	READ(nin, []byte(" %v\n"), nidim)
	if (*nidim) < 1 || (*nidim) > (*nidmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), 'n', (*nidmax))
		goto Label230
	}
	WRITE((*nin), " %v\n", ((*idim)[(*i)-1], (*i) = 1, (*nidim)))
	for (*i) = 1; (*i) <= (*nidim); (*i)++ {
		if (*idim)[(*i)-1] < 0 || (*idim)[(*i)-1] > (*nmax) {
			WRITE((*nout), *func() *[]byte{y := []byte(" value of n is less than 0 or greater than %2d\n"); return &y}(), (*nmax))
			goto Label230
		}
	//Label10:
	}
	//*     Values of K
	READ(nin, []byte(" %v\n"), nkb)
	if (*nkb) < 1 || (*nkb) > (*nkbmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), 'k', (*nkbmax))
		goto Label230
	}
	WRITE((*nin), " %v\n", ((*kb)[(*i)-1], (*i) = 1, (*nkb)))
	for (*i) = 1; (*i) <= (*nkb); (*i)++ {
		if (*kb)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte{y := []byte(" value of k is less than 0\n"); return &y}())
			goto Label230
		}
	//Label20:
	}
	//*     Values of incx and incy
	READ(nin, []byte(" %v\n"), ninc)
	if (*ninc) < 1 || (*ninc) > (*ninmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), *func() *[]byte{y := []byte("incx and incy"); return &y}(), (*ninmax))
		goto Label230
	}
	WRITE((*nin), " %v\n", ((*inc)[(*i)-1], (*i) = 1, (*ninc)))
	for (*i) = 1; (*i) <= (*ninc); (*i)++ {
		if (*inc)[(*i)-1] == 0 || abs ((*inc)[(*i)-1]) > (*incmax) {
			WRITE((*nout), *func() *[]byte{y := []byte(" absolute value of incx or incy is 0 or greater than %2d\n"); return &y}(), (*incmax))
			goto Label230
		}
	//Label30:
	}
	//*     Values of ALPHA
	READ(nin, []byte(" %v\n"), nalf)
	if (*nalf) < 1 || (*nalf) > (*nalmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), *func() *[]byte{y := []byte("alpha"); return &y}(), (*nalmax))
		goto Label230
	}
	WRITE((*nin), " %v\n", ((*alf)[(*i)-1], (*i) = 1, (*nalf)))
	//*     Values of BETA
	READ(nin, []byte(" %v\n"), nbet)
	if (*nbet) < 1 || (*nbet) > (*nbemax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), *func() *[]byte{y := []byte("beta"); return &y}(), (*nbemax))
		goto Label230
	}
	WRITE((*nin), " %v\n", ((*bet)[(*i)-1], (*i) = 1, (*nbet)))
	//*
	//*     Report values of parameters.
	//*
	WRITE((*nout), *func() *[]byte{y := []byte(" tests of the complex          level 2 blas%v the following parameter values will be used:\n"); return &y}())
	WRITE((*nout), "   for n              %6d\n", ((*idim)[(*i)-1], (*i) = 1, (*nidim)))
	WRITE((*nout), "   for k              %6d\n", ((*kb)[(*i)-1], (*i) = 1, (*nkb)))
	WRITE((*nout), "   for incx and incy  %6d\n", ((*inc)[(*i)-1], (*i) = 1, (*ninc)))
	WRITE((*nout), "   for alpha          (%4.1f,%4.1f)  %v\n", ((*alf)[(*i)-1], (*i) = 1, (*nalf)))
	WRITE((*nout), "   for beta           (%4.1f,%4.1f)  %v\n", ((*bet)[(*i)-1], (*i) = 1, (*nbet)))
	if !(*tsterr) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}())
		WRITE((*nout), *func() *[]byte{y := []byte(" error-exits will not be tested\n"); return &y}())
	}
	WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}())
	WRITE((*nout), *func() *[]byte{y := []byte(" routines pass computational tests if test ratio is less than%8.2f\n"); return &y}(), (*thresh))
	WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}())
	//*
	//*     Read names of subroutines and flags which indicate
	//*     whether they are to be tested.
	//*
	for (*i) = 1; (*i) <= (*nsubs); (*i)++ {
		(*ltest)[(*i)-1] = false
	//Label40:
	}
Label50:
	READ(nin, []byte("%6s  %t\n"), snamet, ltestt)
	for (*i) = 1; (*i) <= (*nsubs); (*i)++ {
		if (*snamet) == (*snames)[(*i)-1] {
			goto Label70
		}
	//Label60:
	}
	WRITE((*nout), *func() *[]byte{y := []byte(" subprogram name %6s not recognized\n ******* tests abandoned *******\n"); return &y}(), (*snamet))
	t.Error("")
Label70:
	(*ltest)[(*i)-1] = (*ltestt)
	goto Label50
	//*
//Label80:
	;
	CLOSE(NIN)
	//*
	//*     Compute EPS (the machine precision).
	//*
	(*eps) = (EPSILON((*rzero)))
	WRITE((*nout), *func() *[]byte{y := []byte(" relative machine precision is taken to be%f%9.1e\n"); return &y}(), (*eps))
	//*
	//*     Check the reliability of CMVCH using exact data.
	//*
	(*n) = (MIN(int(32), (*nmax)))
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*a)[(*i)-1][(*j)-1] = (MAX((*i)-(*j)+1, 0))
		//Label110:
		}
		(*x)[(*j)-1] = (*j)
		(*y)[(*j)-1] = (*zero)
	//Label120:
	}
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		(*yy)[(*j)-1] = (*j)*(((*j)+1)*(*j))/2 - (((*j)+1)*(*j)*((*j)-1))/3
	//Label130:
	}
	//*     YY holds the exact result. On exit from CMVCH YT holds
	//*     the result computed by CMVCH.
	(*trans) = 'n'
	Cmvch((*trans), (*n), (*n), (*one), (*a), (*nmax), (*x), 1, (*zero), (*y), 1, (*yt), (*g), (*yy), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*Lce(yy, yt, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in cmvch -  in-line dot products are being evaluated wrongly.\n cmvch was called with trans = %c and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*trans), (*same), (*err))
		t.Error("")
	}
	(*trans) = 't'
	Cmvch((*trans), (*n), (*n), (*one), (*a), (*nmax), (*x), -1, (*zero), (*y), -1, (*yt), (*g), (*yy), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*Lce(yy, yt, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in cmvch -  in-line dot products are being evaluated wrongly.\n cmvch was called with trans = %c and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*trans), (*same), (*err))
		t.Error("")
	}
	//*
	//*     Test each subroutine in turn.
	//*
	for (*isnum) = 1; (*isnum) <= (*nsubs); (*isnum)++ {
		WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}())
		if !(*ltest)[(*isnum) - 1] {
			//*           Subprogram is not to be tested.
			WRITE((*nout), *func() *[]byte{y := []byte(" %6s was not tested\n"); return &y}(), (*snames)[(*isnum)-1])
		} else {
			(*srnamt) = (*snames)[(*isnum)-1]
			//*           Test error exits.
			if *tsterr {
				Cchke(isnum, &((*snames)[(*isnum)-1]), nout)
				WRITE((*nout), *func() *[]byte{y := []byte(" %v\n"); return &y}())
			}
			//*           Test computations.
			(*infot) = 0
			(*ok) = true
			(*fatal) = false
			switch *isnum {
			case 1:
				goto Label140
			case 2:
				goto Label140
			case 3:
				goto Label150
			case 4:
				goto Label150
			case 5:
				goto Label150
			case 6:
				goto Label160
			case 7:
				goto Label160
			case 8:
				goto Label160
			case 9:
				goto Label160
			case 10:
				goto Label160
			case 11:
				goto Label160
			case 12:
				goto Label170
			case 13:
				goto Label170
			case 14:
				goto Label180
			case 15:
				goto Label180
			case 16:
				goto Label190
			case 17:
				goto Label190
			}
			//*           Test Cgemv, 01, and Cgbmv, 02.
		Label140:
			;
			Cchk1(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nkb, kb, nalf, alf, nbet, bet, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g)
			goto Label200
			//*           Test Chemv, 03, Chbmv, 04, and Chpmv, 05.
		Label150:
			;
			Cchk2(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nkb, kb, nalf, alf, nbet, bet, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g)
			goto Label200
			//*           Test Ctrmv, 06, Ctbmv, 07, Ctpmv, 08,
			//*           Ctrsv, 09, Ctbsv, 10, and Ctpsv, 11.
		Label160:
			;
			Cchk3(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nkb, kb, ninc, inc, nmax, incmax, a, aa, as, y, yy, ys, yt, g, z)
			goto Label200
			//*           Test Cgerc, 12, Cgeru, 13.
		Label170:
			;
			Cchk4(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			goto Label200
			//*           Test Cher, 14, and Chpr, 15.
		Label180:
			;
			Cchk5(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			goto Label200
			//*           Test Cher2, 16, and Chpr2, 17.
		Label190:
			;
			Cchk6(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			//*
		Label200:
			;
			if (*fatal) && (*sfatal) {
				goto Label220
			}
		}
	//Label210:
	}
	WRITE((*nout), *func() *[]byte{y := []byte("\n end of tests\n"); return &y}())
	goto Label240
	//*
Label220:
	;
	WRITE((*nout), *func() *[]byte{y := []byte("\n ******* fatal error - tests abandoned *******\n"); return &y}())
	goto Label240
	//*
Label230:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" amend data file or increase array sizes in program\n ******* tests abandoned *******\n"); return &y}())
	//*
Label240:
	;
	if *trace {
		CLOSE(ntra)
	}
	CLOSE(nout)
	t.Error("")
	//*
	//*
	//*     End of Cblat2.
	//*
}

func Cchk1(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nkb *int, kb *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, y *[]complex128, yy *[]complex128, ys *[]complex128, yt *[]complex128, g *[]float64) {
	zero := new(complex128)
	half := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bls := new(complex128)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ic := new(int)
	iku := new(int)
	im := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	incy := new(int)
	incys := new(int)
	ix := new(int)
	iy := new(int)
	kl := new(int)
	kls := new(int)
	ku := new(int)
	kus := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lx := new(int)
	ly := new(int)
	m := new(int)
	ml := new(int)
	ms := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	nd := new(int)
	nk := new(int)
	nl := new(int)
	ns := new(int)
	banded := new(bool)
	full := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	tran := new(bool)
	trans := new(byte)
	transs := new(byte)
	ich := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	//*
	//*  Tests Cgemv and Cgbmv.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Data statements ..
	(*ich)[0], (*ich)[1], (*ich)[2] = 'n', 't', 'c'
	//*     .. Executable Statements ..
	(*full) = (*sname)[3     - 1] == "e"
	(*banded) = (*sname)[3     - 1] == "b"
	//*     Define the number of arguments.
	if *full {
		(*nargs) = 11
	} else if *banded {
		(*nargs) = 13
	}
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		(*nd) = (*n)/2 + 1
		//*
		for (*im) = 1; (*im) <= 2; (*im)++ {
			if (*im) == 1 {
				(*m) = (MAX((*n)-(*nd), 0))
			}
			if (*im) == 2 {
				(*m) = (MIN((*n)+(*nd), (*nmax)))
			}
			//*
			if *banded {
				(*nk) = (*nkb)
			} else {
				(*nk) = 1
			}
			for (*iku) = 1; (*iku) <= (*nk); (*iku)++ {
				if *banded {
					(*ku) = (*kb)[(*iku)-1]
					(*kl) = (MAX((*ku)-1, 0))
				} else {
					(*ku) = (*n) - 1
					(*kl) = (*m) - 1
				}
				//*              Set LDA to 1 more than minimum value if room.
				if *banded {
					(*lda) = (*kl) + (*ku) + 1
				} else {
					(*lda) = (*m)
				}
				if (*lda) < (*nmax) {
					(*lda) = (*lda) + 1
				}
				//*              Skip tests if not enough room.
				if (*lda) > (*nmax) {
					goto Label100
				}
				(*laa) = (*lda) * (*n)
				(*null) = (*n) <= 0 || (*m) <= 0
				//*
				//*              Generate the matrix A.
				//*
				(*transl) = (*zero)
				Cmake(&((*sname)[1]), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, a, nmax, aa, lda, kl, ku, reset, transl)
				//*
				for (*ic) = 1; (*ic) <= 3; (*ic)++ {
					(*trans) = (*ich)[(*ic)-1]
					(*tran) = (*trans) == "t" || (*trans) == "c"
					//*
					if *tran {
						(*ml) = (*n)
						(*nl) = (*m)
					} else {
						(*ml) = (*m)
						(*nl) = (*n)
					}
					//*
					for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
						(*incx) = (*inc)[(*ix)-1]
						(*lx) = ABS((*incx)) * (*nl)
						//*
						//*                    Generate the vector X.
						//*
						(*transl) = (*half)
						Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), nl, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*nl)-1, reset, transl)
						if (*nl) > 1 {
							(*x)[(*nl)/1] = (*zero)
							(*xx)[1+ABS((*incx))*((*nl)/2-1)-1] = (*zero)
						}
						//*
						for (*iy) = 1; (*iy) <= (*ninc); (*iy)++ {
							(*incy) = (*inc)[(*iy)-1]
							(*ly) = ABS((*incy)) * (*ml)
							//*
							for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
								(*alpha) = (*alf)[(*ia)-1]
								//*
								for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
									(*beta) = (*bet)[(*ib)-1]
									//*
									//*                             Generate the vector Y.
									//*
									(*transl) = (*zero)
									Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), ml, y, func() *int{y := 1; return &y}(), yy, ABS((*incy)), func() *int{y := 0; return &y}(), (*ml)-1, reset, transl)
									//*
									(*nc) = (*nc) + 1
									//*
									//*                             Save every datum before calling the
									//*                             subroutine.
									//*
									(*transs) = (*trans)
									(*ms) = (*m)
									(*ns) = (*n)
									(*kls) = (*kl)
									(*kus) = (*ku)
									(*als) = (*alpha)
									for (*i) = 1; (*i) <= (*laa); (*i)++ {
										(*as)[(*i)-1] = (*aa)[(*i)-1]
									//Label10:
									}
									(*ldas) = (*lda)
									for (*i) = 1; (*i) <= (*lx); (*i)++ {
										(*xs)[(*i)-1] = (*xx)[(*i)-1]
									//Label20:
									}
									(*incxs) = (*incx)
									(*bls) = (*beta)
									for (*i) = 1; (*i) <= (*ly); (*i)++ {
										(*ys)[(*i)-1] = (*yy)[(*i)-1]
									//Label30:
									}
									(*incys) = (*incy)
									//*
									//*                             Call the subroutine.
									//*
									if *full {
										if *trace {
											WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)         .\n"); return &y}(), (*nc), (*sname), (*trans), (*m), (*n), (*alpha), (*lda), (*incx), (*beta), (*incy))
										}
										if *rewi {
											REWIND(ntra)
										}
										Cgemv(trans, m, n, alpha, aa, lda, xx, incx, beta, yy, incy)
									} else if *banded {
										if *trace {
											WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d) .\n"); return &y}(), (*nc), (*sname), (*trans), (*m), (*n), (*kl), (*ku), (*alpha), (*lda), (*incx), (*beta), (*incy))
										}
										if *rewi {
											REWIND(ntra)
										}
										Cgbmv(trans, m, n, kl, ku, alpha, aa, lda, xx, incx, beta, yy, incy)
									}
									//*
									//*                             Check if error-exit was taken incorrectly.
									//*
									if !(*ok) {
										WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
										(*fatal) = true
										goto Label130
									}
									//*
									//*                             See what data changed inside subroutines.
									//*
									(*isame)[0] = (*trans) == (*transs)
									(*isame)[1] = (*ms) == (*m)
									(*isame)[2] = (*ns) == (*n)
									if *full {
										(*isame)[3] = (*als) == (*alpha)
										(*isame)[4] = (*Lce(as, aa, laa))
										(*isame)[5] = (*ldas) == (*lda)
										(*isame)[6] = (*Lce(xs, xx, lx))
										(*isame)[7] = (*incxs) == (*incx)
										(*isame)[8] = (*bls) == (*beta)
										if *null {
											(*isame)[9] = (*Lce(ys, yy, ly))
										} else {
											(*isame)[9] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), ml, ys, yy, ABS((*incy))))
										}
										(*isame)[10] = (*incys) == (*incy)
									} else if *banded {
										(*isame)[3] = (*kls) == (*kl)
										(*isame)[4] = (*kus) == (*ku)
										(*isame)[5] = (*als) == (*alpha)
										(*isame)[6] = (*Lce(as, aa, laa))
										(*isame)[7] = (*ldas) == (*lda)
										(*isame)[8] = (*Lce(xs, xx, lx))
										(*isame)[9] = (*incxs) == (*incx)
										(*isame)[10] = (*bls) == (*beta)
										if *null {
											(*isame)[11] = (*Lce(ys, yy, ly))
										} else {
											(*isame)[11] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), ml, ys, yy, ABS((*incy))))
										}
										(*isame)[12] = (*incys) == (*incy)
									}
									//*
									//*                             If data was incorrectly changed, report
									//*                             and return.
									//*
									(*same) = true
									for (*i) = 1; (*i) <= (*nargs); (*i)++ {
										(*same) = (*same) && (*isame)[(*i)-1]
										if !(*isame)[(*i)-1] {
											WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
										}
									//Label40:
									}
									if !(*same) {
										(*fatal) = true
										goto Label130
									}
									//*
									if !(*null) {
										//*
										//*                                Check the result.
										//*
										Cmvch((*trans), (*m), (*n), (*alpha), (*a), (*nmax), (*x), (*incx), (*beta), (*y), (*incy), (*yt), (*g), (*yy), (*eps), (*err), (*fatal), (*nout), true)
										(*errmax) = (MAX((*errmax), (*err)))
										//*                                If got really bad answer, report and
										//*                                return.
										if *fatal {
											goto Label130
										}
									} else {
										//*                                Avoid repeating tests with M.le.0 or
										//*                                N.le.0.
										goto Label110
									}
									//*
								//Label50:
								}
								//*
							//Label60:
							}
							//*
						//Label70:
						}
						//*
					//Label80:
					}
					//*
				//Label90:
				}
				//*
			Label100:
			}
			//*
		Label110:
		}
		//*
	//Label120:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label140
	//*
Label130:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *full {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)         .\n"); return &y}(), (*nc), (*sname), (*trans), (*m), (*n), (*alpha), (*lda), (*incx), (*beta), (*incy))
	} else if *banded {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d) .\n"); return &y}(), (*nc), (*sname), (*trans), (*m), (*n), (*kl), (*ku), (*alpha), (*lda), (*incx), (*beta), (*incy))
	}
	//*
Label140:
	;
	return
	//*
	//*
	//*     End of Cchk1.
	//*
}

func Cchk2(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nkb *int, kb *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, y *[]complex128, yy *[]complex128, ys *[]complex128, yt *[]complex128, g *[]float64) {
	zero := new(complex128)
	half := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bls := new(complex128)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ic := new(int)
	ik := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	incy := new(int)
	incys := new(int)
	ix := new(int)
	iy := new(int)
	k := new(int)
	ks := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lx := new(int)
	ly := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	nk := new(int)
	ns := new(int)
	banded := new(bool)
	full := new(bool)
	null := new(bool)
	packed := new(bool)
	reset := new(bool)
	same := new(bool)
	uplo := new(byte)
	uplos := new(byte)
	ich := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Chemv, Chbmv and Chpmv.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Data statements ..
	(*ich)[0], (*ich)[1] = 'u', 'l'
	//*     .. Executable Statements ..
	(*full) = (*sname)[3     - 1] == "e"
	(*banded) = (*sname)[3     - 1] == "b"
	(*packed) = (*sname)[3     - 1] == "p"
	//*     Define the number of arguments.
	if *full {
		(*nargs) = 10
	} else if *banded {
		(*nargs) = 11
	} else if *packed {
		(*nargs) = 9
	}
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*
		if *banded {
			(*nk) = (*nkb)
		} else {
			(*nk) = 1
		}
		for (*ik) = 1; (*ik) <= (*nk); (*ik)++ {
			if *banded {
				(*k) = (*kb)[(*ik)-1]
			} else {
				(*k) = (*n) - 1
			}
			//*           Set LDA to 1 more than minimum value if room.
			if *banded {
				(*lda) = (*k) + 1
			} else {
				(*lda) = (*n)
			}
			if (*lda) < (*nmax) {
				(*lda) = (*lda) + 1
			}
			//*           Skip tests if not enough room.
			if (*lda) > (*nmax) {
				goto Label100
			}
			if *packed {
				(*laa) = ((*n) * ((*n) + 1)) / 2
			} else {
				(*laa) = (*lda) * (*n)
			}
			(*null) = (*n) <= 0
			//*
			for (*ic) = 1; (*ic) <= 2; (*ic)++ {
				(*uplo) = (*ich)[(*ic)-1]
				//*
				//*              Generate the matrix A.
				//*
				(*transl) = (*zero)
				Cmake(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), n, n, a, nmax, aa, lda, k, k, reset, transl)
				//*
				for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
					(*incx) = (*inc)[(*ix)-1]
					(*lx) = ABS((*incx)) * (*n)
					//*
					//*                 Generate the vector X.
					//*
					(*transl) = (*half)
					Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
					if (*n) > 1 {
						(*x)[(*n)/1] = (*zero)
						(*xx)[1+ABS((*incx))*((*n)/2-1)-1] = (*zero)
					}
					//*
					for (*iy) = 1; (*iy) <= (*ninc); (*iy)++ {
						(*incy) = (*inc)[(*iy)-1]
						(*ly) = ABS((*incy)) * (*n)
						//*
						for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
							(*alpha) = (*alf)[(*ia)-1]
							//*
							for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
								(*beta) = (*bet)[(*ib)-1]
								//*
								//*                          Generate the vector Y.
								//*
								(*transl) = (*zero)
								Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, y, func() *int{y := 1; return &y}(), yy, ABS((*incy)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
								//*
								(*nc) = (*nc) + 1
								//*
								//*                          Save every datum before calling the
								//*                          subroutine.
								//*
								(*uplos) = (*uplo)
								(*ns) = (*n)
								(*ks) = (*k)
								(*als) = (*alpha)
								for (*i) = 1; (*i) <= (*laa); (*i)++ {
									(*as)[(*i)-1] = (*aa)[(*i)-1]
								//Label10:
								}
								(*ldas) = (*lda)
								for (*i) = 1; (*i) <= (*lx); (*i)++ {
									(*xs)[(*i)-1] = (*xx)[(*i)-1]
								//Label20:
								}
								(*incxs) = (*incx)
								(*bls) = (*beta)
								for (*i) = 1; (*i) <= (*ly); (*i)++ {
									(*ys)[(*i)-1] = (*yy)[(*i)-1]
								//Label30:
								}
								(*incys) = (*incy)
								//*
								//*                          Call the subroutine.
								//*
								if *full {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)             .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*lda), (*incx), (*beta), (*incy))
									}
									if *rewi {
										REWIND(ntra)
									}
									Chemv(uplo, n, alpha, aa, lda, xx, incx, beta, yy, incy)
								} else if *banded {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*k), (*alpha), (*lda), (*incx), (*beta), (*incy))
									}
									if *rewi {
										REWIND(ntra)
									}
									Chbmv(uplo, n, k, alpha, aa, lda, xx, incx, beta, yy, incy)
								} else if *packed {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), ap, x,%2d,(%4.1f,%4.1f), y,%2d)                .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*beta), (*incy))
									}
									if *rewi {
										REWIND(ntra)
									}
									Chpmv(uplo, n, alpha, aa, xx, incx, beta, yy, incy)
								}
								//*
								//*                          Check if error-exit was taken incorrectly.
								//*
								if !(*ok) {
									WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
									(*fatal) = true
									goto Label120
								}
								//*
								//*                          See what data changed inside subroutines.
								//*
								(*isame)[0] = (*uplo) == (*uplos)
								(*isame)[1] = (*ns) == (*n)
								if *full {
									(*isame)[2] = (*als) == (*alpha)
									(*isame)[3] = (*Lce(as, aa, laa))
									(*isame)[4] = (*ldas) == (*lda)
									(*isame)[5] = (*Lce(xs, xx, lx))
									(*isame)[6] = (*incxs) == (*incx)
									(*isame)[7] = (*bls) == (*beta)
									if *null {
										(*isame)[8] = (*Lce(ys, yy, ly))
									} else {
										(*isame)[8] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, ys, yy, ABS((*incy))))
									}
									(*isame)[9] = (*incys) == (*incy)
								} else if *banded {
									(*isame)[2] = (*ks) == (*k)
									(*isame)[3] = (*als) == (*alpha)
									(*isame)[4] = (*Lce(as, aa, laa))
									(*isame)[5] = (*ldas) == (*lda)
									(*isame)[6] = (*Lce(xs, xx, lx))
									(*isame)[7] = (*incxs) == (*incx)
									(*isame)[8] = (*bls) == (*beta)
									if *null {
										(*isame)[9] = (*Lce(ys, yy, ly))
									} else {
										(*isame)[9] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, ys, yy, ABS((*incy))))
									}
									(*isame)[10] = (*incys) == (*incy)
								} else if *packed {
									(*isame)[2] = (*als) == (*alpha)
									(*isame)[3] = (*Lce(as, aa, laa))
									(*isame)[4] = (*Lce(xs, xx, lx))
									(*isame)[5] = (*incxs) == (*incx)
									(*isame)[6] = (*bls) == (*beta)
									if *null {
										(*isame)[7] = (*Lce(ys, yy, ly))
									} else {
										(*isame)[7] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, ys, yy, ABS((*incy))))
									}
									(*isame)[8] = (*incys) == (*incy)
								}
								//*
								//*                          If data was incorrectly changed, report and
								//*                          return.
								//*
								(*same) = true
								for (*i) = 1; (*i) <= (*nargs); (*i)++ {
									(*same) = (*same) && (*isame)[(*i)-1]
									if !(*isame)[(*i)-1] {
										WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
									}
								//Label40:
								}
								if !(*same) {
									(*fatal) = true
									goto Label120
								}
								//*
								if !(*null) {
									//*
									//*                             Check the result.
									//*
									Cmvch( "n", (*n), (*n), (*alpha), (*a), (*nmax), (*x), (*incx), (*beta), (*y), (*incy), (*yt), (*g), (*yy), (*eps), (*err), (*fatal), (*nout), true)
									(*errmax) = (MAX((*errmax), (*err)))
									//*                             If got really bad answer, report and
									//*                             return.
									if *fatal {
										goto Label120
									}
								} else {
									//*                             Avoid repeating tests with N.le.0
									goto Label110
								}
								//*
							//Label50:
							}
							//*
						//Label60:
						}
						//*
					//Label70:
					}
					//*
				//Label80:
				}
				//*
			//Label90:
			}
			//*
		Label100:
		}
		//*
	Label110:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label130
	//*
Label120:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *full {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)             .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*lda), (*incx), (*beta), (*incy))
	} else if *banded {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, x,%2d,(%4.1f,%4.1f), y,%2d)         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*k), (*alpha), (*lda), (*incx), (*beta), (*incy))
	} else if *packed {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), ap, x,%2d,(%4.1f,%4.1f), y,%2d)                .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*beta), (*incy))
	}
	//*
Label130:
	;
	return
	//*
	//*
	//*     End of Cchk2.
	//*
}

func Cchk3(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nkb *int, kb *[]int, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, xt *[]complex128, g *[]float64, z *[]complex128) {
	zero := new(complex128)
	half := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	icd := new(int)
	ict := new(int)
	icu := new(int)
	ik := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	ix := new(int)
	k := new(int)
	ks := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lx := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	nk := new(int)
	ns := new(int)
	banded := new(bool)
	full := new(bool)
	null := new(bool)
	packed := new(bool)
	reset := new(bool)
	same := new(bool)
	diag := new(byte)
	diags := new(byte)
	trans := new(byte)
	transs := new(byte)
	uplo := new(byte)
	uplos := new(byte)
	ichd := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	ichu := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	icht := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Ctrmv, Ctbmv, Ctpmv, Ctrsv, Ctbsv and Ctpsv.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Data statements ..
	(*ichu)[0], (*ichu)[1], (*icht)[0], (*icht)[1], (*icht)[2], (*ichd)[0], (*ichd)[1] = 'u', 'l', 'n', 't', 'c', 'u', 'n'
	//*     .. Executable Statements ..
	(*full) = (*sname)[3     - 1] == "r"
	(*banded) = (*sname)[3     - 1] == "b"
	(*packed) = (*sname)[3     - 1] == "p"
	//*     Define the number of arguments.
	if *full {
		(*nargs) = 8
	} else if *banded {
		(*nargs) = 9
	} else if *packed {
		(*nargs) = 7
	}
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*     Set up zero vector for CMVCH.
	for (*i) = 1; (*i) <= (*nmax); (*i)++ {
		(*z)[(*i)-1] = (*zero)
	//Label10:
	}
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*
		if *banded {
			(*nk) = (*nkb)
		} else {
			(*nk) = 1
		}
		for (*ik) = 1; (*ik) <= (*nk); (*ik)++ {
			if *banded {
				(*k) = (*kb)[(*ik)-1]
			} else {
				(*k) = (*n) - 1
			}
			//*           Set LDA to 1 more than minimum value if room.
			if *banded {
				(*lda) = (*k) + 1
			} else {
				(*lda) = (*n)
			}
			if (*lda) < (*nmax) {
				(*lda) = (*lda) + 1
			}
			//*           Skip tests if not enough room.
			if (*lda) > (*nmax) {
				goto Label100
			}
			if *packed {
				(*laa) = ((*n) * ((*n) + 1)) / 2
			} else {
				(*laa) = (*lda) * (*n)
			}
			(*null) = (*n) <= 0
			//*
			for (*icu) = 1; (*icu) <= 2; (*icu)++ {
				(*uplo) = (*ichu)[(*icu)-1]
				//*
				for (*ict) = 1; (*ict) <= 3; (*ict)++ {
					(*trans) = (*icht)[(*ict)-1]
					//*
					for (*icd) = 1; (*icd) <= 2; (*icd)++ {
						(*diag) = (*ichd)[(*icd)-1]
						//*
						//*                    Generate the matrix A.
						//*
						(*transl) = (*zero)
						Cmake(&((*sname)[1]), uplo, diag, n, n, a, nmax, aa, lda, k, k, reset, transl)
						//*
						for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
							(*incx) = (*inc)[(*ix)-1]
							(*lx) = ABS((*incx)) * (*n)
							//*
							//*                       Generate the vector X.
							//*
							(*transl) = (*half)
							Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
							if (*n) > 1 {
								(*x)[(*n)/1] = (*zero)
								(*xx)[1+ABS((*incx))*((*n)/2-1)-1] = (*zero)
							}
							//*
							(*nc) = (*nc) + 1
							//*
							//*                       Save every datum before calling the subroutine.
							//*
							(*uplos) = (*uplo)
							(*transs) = (*trans)
							(*diags) = (*diag)
							(*ns) = (*n)
							(*ks) = (*k)
							for (*i) = 1; (*i) <= (*laa); (*i)++ {
								(*as)[(*i)-1] = (*aa)[(*i)-1]
							//Label20:
							}
							(*ldas) = (*lda)
							for (*i) = 1; (*i) <= (*lx); (*i)++ {
								(*xs)[(*i)-1] = (*xx)[(*i)-1]
							//Label30:
							}
							(*incxs) = (*incx)
							//*
							//*                       Call the subroutine.
							//*
							if (*sname)[4     - 1] == "mv" {
								if *full {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                                   .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*lda), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctrmv(uplo, trans, diag, n, aa, lda, xx, incx)
								} else if *banded {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                               .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*k), (*lda), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctbmv(uplo, trans, diag, n, k, aa, lda, xx, incx)
								} else if *packed {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, ap, x,%2d)                                      .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctpmv(uplo, trans, diag, n, aa, xx, incx)
								}
							} else if (*sname)[4     - 1] == "sv" {
								if *full {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                                   .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*lda), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctrsv(uplo, trans, diag, n, aa, lda, xx, incx)
								} else if *banded {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                               .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*k), (*lda), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctbsv(uplo, trans, diag, n, k, aa, lda, xx, incx)
								} else if *packed {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, ap, x,%2d)                                      .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*incx))
									}
									if *rewi {
										REWIND(ntra)
									}
									Ctpsv(uplo, trans, diag, n, aa, xx, incx)
								}
							}
							//*
							//*                       Check if error-exit was taken incorrectly.
							//*
							if !(*ok) {
								WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
								(*fatal) = true
								goto Label120
							}
							//*
							//*                       See what data changed inside subroutines.
							//*
							(*isame)[0] = (*uplo) == (*uplos)
							(*isame)[1] = (*trans) == (*transs)
							(*isame)[2] = (*diag) == (*diags)
							(*isame)[3] = (*ns) == (*n)
							if *full {
								(*isame)[4] = (*Lce(as, aa, laa))
								(*isame)[5] = (*ldas) == (*lda)
								if *null {
									(*isame)[6] = (*Lce(xs, xx, lx))
								} else {
									(*isame)[6] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, xs, xx, ABS((*incx))))
								}
								(*isame)[7] = (*incxs) == (*incx)
							} else if *banded {
								(*isame)[4] = (*ks) == (*k)
								(*isame)[5] = (*Lce(as, aa, laa))
								(*isame)[6] = (*ldas) == (*lda)
								if *null {
									(*isame)[7] = (*Lce(xs, xx, lx))
								} else {
									(*isame)[7] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, xs, xx, ABS((*incx))))
								}
								(*isame)[8] = (*incxs) == (*incx)
							} else if *packed {
								(*isame)[4] = (*Lce(as, aa, laa))
								if *null {
									(*isame)[5] = (*Lce(xs, xx, lx))
								} else {
									(*isame)[5] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, xs, xx, ABS((*incx))))
								}
								(*isame)[6] = (*incxs) == (*incx)
							}
							//*
							//*                       If data was incorrectly changed, report and
							//*                       return.
							//*
							(*same) = true
							for (*i) = 1; (*i) <= (*nargs); (*i)++ {
								(*same) = (*same) && (*isame)[(*i)-1]
								if !(*isame)[(*i)-1] {
									WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
								}
							//Label40:
							}
							if !(*same) {
								(*fatal) = true
								goto Label120
							}
							//*
							if !(*null) {
								if (*sname)[4     - 1] == "mv" {
									//*
									//*                             Check the result.
									//*
									Cmvch((*trans), (*n), (*n), (*one), (*a), (*nmax), (*x), (*incx), (*zero), (*z), (*incx), (*xt), (*g), (*xx), (*eps), (*err), (*fatal), (*nout), true)
								} else if (*sname)[4     - 1] == "sv" {
									//*
									//*                             Compute approximation to original vector.
									//*
									for (*i) = 1; (*i) <= (*n); (*i)++ {
										(*z)[(*i)-1] = (*xx)[1+((*i)-1)*ABS((*incx))-1]
										(*xx)[1+((*i)-1)*ABS((*incx))-1] = (*x)[(*i)-1]
									//Label50:
									}
									Cmvch((*trans), (*n), (*n), (*one), (*a), (*nmax), (*z), (*incx), (*zero), (*x), (*incx), (*xt), (*g), (*xx), (*eps), (*err), (*fatal), (*nout), false)
								}
								(*errmax) = (MAX((*errmax), (*err)))
								//*                          If got really bad answer, report and return.
								if *fatal {
									goto Label120
								}
							} else {
								//*                          Avoid repeating tests with N.le.0.
								goto Label110
							}
							//*
						//Label60:
						}
						//*
					//Label70:
					}
					//*
				//Label80:
				}
				//*
			//Label90:
			}
			//*
		Label100:
		}
		//*
	Label110:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label130
	//*
Label120:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *full {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                                   .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*lda), (*incx))
	} else if *banded {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, a,%3d, x,%2d)                               .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*k), (*lda), (*incx))
	} else if *packed {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d, ap, x,%2d)                                      .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*diag), (*n), (*incx))
	}
	//*
Label130:
	;
	return
	//*
	//*
	//*     End of Cchk3.
	//*
}

func Cchk4(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, y *[]complex128, yy *[]complex128, ys *[]complex128, yt *[]complex128, g *[]float64, z *[]complex128) {
	zero := new(complex128)
	half := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	im := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	incy := new(int)
	incys := new(int)
	ix := new(int)
	iy := new(int)
	j := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lx := new(int)
	ly := new(int)
	m := new(int)
	ms := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	nd := new(int)
	ns := new(int)
	conj := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	w := func() *[]complex128 {
		arr := make([]complex128, 1)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Cgerc and Cgeru.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Executable Statements ..
	(*conj) = (*sname)[5     - 1] == "c"
	//*     Define the number of arguments.
	(*nargs) = 9
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		(*nd) = (*n)/2 + 1
		//*
		for (*im) = 1; (*im) <= 2; (*im)++ {
			if (*im) == 1 {
				(*m) = (MAX((*n)-(*nd), 0))
			}
			if (*im) == 2 {
				(*m) = (MIN((*n)+(*nd), (*nmax)))
			}
			//*
			//*           Set LDA to 1 more than minimum value if room.
			(*lda) = (*m)
			if (*lda) < (*nmax) {
				(*lda) = (*lda) + 1
			}
			//*           Skip tests if not enough room.
			if (*lda) > (*nmax) {
				goto Label110
			}
			(*laa) = (*lda) * (*n)
			(*null) = (*n) <= 0 || (*m) <= 0
			//*
			for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
				(*incx) = (*inc)[(*ix)-1]
				(*lx) = ABS((*incx)) * (*m)
				//*
				//*              Generate the vector X.
				//*
				(*transl) = (*half)
				Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), m, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*m)-1, reset, transl)
				if (*m) > 1 {
					(*x)[(*m)/1] = (*zero)
					(*xx)[1+ABS((*incx))*((*m)/2-1)-1] = (*zero)
				}
				//*
				for (*iy) = 1; (*iy) <= (*ninc); (*iy)++ {
					(*incy) = (*inc)[(*iy)-1]
					(*ly) = ABS((*incy)) * (*n)
					//*
					//*                 Generate the vector Y.
					//*
					(*transl) = (*zero)
					Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, y, func() *int{y := 1; return &y}(), yy, ABS((*incy)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
					if (*n) > 1 {
						(*y)[(*n)/1] = (*zero)
						(*yy)[1+ABS((*incy))*((*n)/2-1)-1] = (*zero)
					}
					//*
					for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
						(*alpha) = (*alf)[(*ia)-1]
						//*
						//*                    Generate the matrix A.
						//*
						(*transl) = (*zero)
						Cmake(&((*sname)[1]), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, a, nmax, aa, lda, (*m)-1, (*n)-1, reset, transl)
						//*
						(*nc) = (*nc) + 1
						//*
						//*                    Save every datum before calling the subroutine.
						//*
						(*ms) = (*m)
						(*ns) = (*n)
						(*als) = (*alpha)
						for (*i) = 1; (*i) <= (*laa); (*i)++ {
							(*as)[(*i)-1] = (*aa)[(*i)-1]
						//Label10:
						}
						(*ldas) = (*lda)
						for (*i) = 1; (*i) <= (*lx); (*i)++ {
							(*xs)[(*i)-1] = (*xx)[(*i)-1]
						//Label20:
						}
						(*incxs) = (*incx)
						for (*i) = 1; (*i) <= (*ly); (*i)++ {
							(*ys)[(*i)-1] = (*yy)[(*i)-1]
						//Label30:
						}
						(*incys) = (*incy)
						//*
						//*                    Call the subroutine.
						//*
						if *trace {
							WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s(%3d,(%4.1f,%4.1f), x,%2d, y,%2d, a,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*m), (*n), (*alpha), (*incx), (*incy), (*lda))
						}
						if *conj {
							if *rewi {
								REWIND(ntra)
							}
							Cgerc(m, n, alpha, xx, incx, yy, incy, aa, lda)
						} else {
							if *rewi {
								REWIND(ntra)
							}
							Cgeru(m, n, alpha, xx, incx, yy, incy, aa, lda)
						}
						//*
						//*                    Check if error-exit was taken incorrectly.
						//*
						if !(*ok) {
							WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
							(*fatal) = true
							goto Label140
						}
						//*
						//*                    See what data changed inside subroutine.
						//*
						(*isame)[0] = (*ms) == (*m)
						(*isame)[1] = (*ns) == (*n)
						(*isame)[2] = (*als) == (*alpha)
						(*isame)[3] = (*Lce(xs, xx, lx))
						(*isame)[4] = (*incxs) == (*incx)
						(*isame)[5] = (*Lce(ys, yy, ly))
						(*isame)[6] = (*incys) == (*incy)
						if *null {
							(*isame)[7] = (*Lce(as, aa, laa))
						} else {
							(*isame)[7] = (*Lceres(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, as, aa, lda))
						}
						(*isame)[8] = (*ldas) == (*lda)
						//*
						//*                    If data was incorrectly changed, report and return.
						//*
						(*same) = true
						for (*i) = 1; (*i) <= (*nargs); (*i)++ {
							(*same) = (*same) && (*isame)[(*i)-1]
							if !(*isame)[(*i)-1] {
								WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
							}
						//Label40:
						}
						if !(*same) {
							(*fatal) = true
							goto Label140
						}
						//*
						if !(*null) {
							//*
							//*                       Check the result column by column.
							//*
							if (*incx) > 0 {
								for (*i) = 1; (*i) <= (*m); (*i)++ {
									(*z)[(*i)-1] = (*x)[(*i)-1]
								//Label50:
								}
							} else {
								for (*i) = 1; (*i) <= (*m); (*i)++ {
									(*z)[(*i)-1] = (*x)[(*m)-(*i)+0]
								//Label60:
								}
							}
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								if (*incy) > 0 {
									(*w)[0] = (*y)[(*j)-1]
								} else {
									(*w)[0] = (*y)[(*n)-(*j)+0]
								}
								if *conj {
									(*w)[0] = (CONJG(((*w)[0])))
								}
								Cmvch( "n", (*m), 1, (*alpha), (*z), (*nmax), (*w), 1, (*one), (*a)[0][(*j) - 1], 1, (*yt), (*g), (*aa)[1 + ((*j) - 1) * (*lda) - 1], (*eps), (*err), (*fatal), (*nout), true)
								(*errmax) = (MAX((*errmax), (*err)))
								//*                          If got really bad answer, report and return.
								if *fatal {
									goto Label130
								}
							//Label70:
							}
						} else {
							//*                       Avoid repeating tests with M.le.0 or N.le.0.
							goto Label110
						}
						//*
					//Label80:
					}
					//*
				//Label90:
				}
				//*
			//Label100:
			}
			//*
		Label110:
		}
		//*
	//Label120:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label150
	//*
Label130:
	;
	WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	//*
Label140:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s(%3d,(%4.1f,%4.1f), x,%2d, y,%2d, a,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*m), (*n), (*alpha), (*incx), (*incy), (*lda))
	//*
Label150:
	;
	return
	//*
	//*
	//*     End of Cchk4.
	//*
}

func Cchk5(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, y *[]complex128, yy *[]complex128, ys *[]complex128, yt *[]complex128, g *[]float64, z *[]complex128) {
	zero := new(complex128)
	half := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	ralpha := new(float64)
	rals := new(float64)
	i := new(int)
	ia := new(int)
	ic := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	ix := new(int)
	j := new(int)
	ja := new(int)
	jj := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lj := new(int)
	lx := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	full := new(bool)
	null := new(bool)
	packed := new(bool)
	reset := new(bool)
	same := new(bool)
	upper := new(bool)
	uplo := new(byte)
	uplos := new(byte)
	ich := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	w := func() *[]complex128 {
		arr := make([]complex128, 1)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Cher and Chpr.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Data statements ..
	(*ich)[0], (*ich)[1] = 'u', 'l'
	//*     .. Executable Statements ..
	(*full) = (*sname)[3     - 1] == "e"
	(*packed) = (*sname)[3     - 1] == "p"
	//*     Define the number of arguments.
	if *full {
		(*nargs) = 7
	} else if *packed {
		(*nargs) = 6
	}
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*        Set LDA to 1 more than minimum value if room.
		(*lda) = (*n)
		if (*lda) < (*nmax) {
			(*lda) = (*lda) + 1
		}
		//*        Skip tests if not enough room.
		if (*lda) > (*nmax) {
			goto Label100
		}
		if *packed {
			(*laa) = ((*n) * ((*n) + 1)) / 2
		} else {
			(*laa) = (*lda) * (*n)
		}
		//*
		for (*ic) = 1; (*ic) <= 2; (*ic)++ {
			(*uplo) = (*ich)[(*ic)-1]
			(*upper) = (*uplo) == "u"
			//*
			for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
				(*incx) = (*inc)[(*ix)-1]
				(*lx) = ABS((*incx)) * (*n)
				//*
				//*              Generate the vector X.
				//*
				(*transl) = (*half)
				Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
				if (*n) > 1 {
					(*x)[(*n)/1] = (*zero)
					(*xx)[1+ABS((*incx))*((*n)/2-1)-1] = (*zero)
				}
				//*
				for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
					(*ralpha) = (*real(((*alf)[(*ia)-1])))
					(*alpha) = (*cmplx(ralpha, rzero))
					(*null) = (*n) <= 0 || (*ralpha) == (*rzero)
					//*
					//*                 Generate the matrix A.
					//*
					(*transl) = (*zero)
					Cmake(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), n, n, a, nmax, aa, lda, (*n)-1, (*n)-1, reset, transl)
					//*
					(*nc) = (*nc) + 1
					//*
					//*                 Save every datum before calling the subroutine.
					//*
					(*uplos) = (*uplo)
					(*ns) = (*n)
					(*rals) = (*ralpha)
					for (*i) = 1; (*i) <= (*laa); (*i)++ {
						(*as)[(*i)-1] = (*aa)[(*i)-1]
					//Label10:
					}
					(*ldas) = (*lda)
					for (*i) = 1; (*i) <= (*lx); (*i)++ {
						(*xs)[(*i)-1] = (*xx)[(*i)-1]
					//Label20:
					}
					(*incxs) = (*incx)
					//*
					//*                 Call the subroutine.
					//*
					if *full {
						if *trace {
							WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, x,%2d, a,%3d)                                      .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*ralpha), (*incx), (*lda))
						}
						if *rewi {
							REWIND(ntra)
						}
						Cher(uplo, n, ralpha, xx, incx, aa, lda)
					} else if *packed {
						if *trace {
							WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, x,%2d, ap)                                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*ralpha), (*incx))
						}
						if *rewi {
							REWIND(ntra)
						}
						Chpr(uplo, n, ralpha, xx, incx, aa)
					}
					//*
					//*                 Check if error-exit was taken incorrectly.
					//*
					if !(*ok) {
						WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
						(*fatal) = true
						goto Label120
					}
					//*
					//*                 See what data changed inside subroutines.
					//*
					(*isame)[0] = (*uplo) == (*uplos)
					(*isame)[1] = (*ns) == (*n)
					(*isame)[2] = (*rals) == (*ralpha)
					(*isame)[3] = (*Lce(xs, xx, lx))
					(*isame)[4] = (*incxs) == (*incx)
					if *null {
						(*isame)[5] = (*Lce(as, aa, laa))
					} else {
						(*isame)[5] = (*Lceres(&((*sname)[1]), uplo, n, n, as, aa, lda))
					}
					if !(*packed) {
						(*isame)[6] = (*ldas) == (*lda)
					}
					//*
					//*                 If data was incorrectly changed, report and return.
					//*
					(*same) = true
					for (*i) = 1; (*i) <= (*nargs); (*i)++ {
						(*same) = (*same) && (*isame)[(*i)-1]
						if !(*isame)[(*i)-1] {
							WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
						}
					//Label30:
					}
					if !(*same) {
						(*fatal) = true
						goto Label120
					}
					//*
					if !(*null) {
						//*
						//*                    Check the result column by column.
						//*
						if (*incx) > 0 {
							for (*i) = 1; (*i) <= (*n); (*i)++ {
								(*z)[(*i)-1] = (*x)[(*i)-1]
							//Label40:
							}
						} else {
							for (*i) = 1; (*i) <= (*n); (*i)++ {
								(*z)[(*i)-1] = (*x)[(*n)-(*i)+0]
							//Label50:
							}
						}
						(*ja) = 1
						for (*j) = 1; (*j) <= (*n); (*j)++ {
							(*w)[0] = (CONJG(((*z)[(*j)-1])))
							if *upper {
								(*jj) = 1
								(*lj) = (*j)
							} else {
								(*jj) = (*j)
								(*lj) = (*n) - (*j) + 1
							}
							Cmvch( "n", (*lj), 1, (*alpha), (*z)[(*jj) - 1], (*lj), (*w), 1, (*one), (*a)[(*jj) - 1][(*j) - 1], 1, (*yt), (*g), (*aa)[(*ja) - 1], (*eps), (*err), (*fatal), (*nout), true)
							if *full {
								if *upper {
									(*ja) = (*ja) + (*lda)
								} else {
									(*ja) = (*ja) + (*lda) + 1
								}
							} else {
								(*ja) = (*ja) + (*lj)
							}
							(*errmax) = (MAX((*errmax), (*err)))
							//*                       If got really bad answer, report and return.
							if *fatal {
								goto Label110
							}
						//Label60:
						}
					} else {
						//*                    Avoid repeating tests if N.le.0.
						if (*n) <= 0 {
							goto Label100
						}
					}
					//*
				//Label70:
				}
				//*
			//Label80:
			}
			//*
		//Label90:
		}
		//*
	Label100:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label130
	//*
Label110:
	;
	WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	//*
Label120:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *full {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, x,%2d, a,%3d)                                      .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*ralpha), (*incx), (*lda))
	} else if *packed {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, x,%2d, ap)                                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*ralpha), (*incx))
	}
	//*
Label130:
	;
	return
	//*
	//*
	//*     End of Cchk5.
	//*
}

func Cchk6(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, ninc *int, inc *[]int, nmax *int, incmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, x *[]complex128, xx *[]complex128, xs *[]complex128, y *[]complex128, yy *[]complex128, ys *[]complex128, yt *[]complex128, g *[]float64, z *[][]complex128) {
	zero := new(complex128)
	half := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	transl := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	ic := new(int)
	in := new(int)
	incx := new(int)
	incxs := new(int)
	incy := new(int)
	incys := new(int)
	ix := new(int)
	iy := new(int)
	j := new(int)
	ja := new(int)
	jj := new(int)
	laa := new(int)
	lda := new(int)
	ldas := new(int)
	lj := new(int)
	lx := new(int)
	ly := new(int)
	n := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	full := new(bool)
	null := new(bool)
	packed := new(bool)
	reset := new(bool)
	same := new(bool)
	upper := new(bool)
	uplo := new(byte)
	uplos := new(byte)
	ich := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	w := func() *[]complex128 {
		arr := make([]complex128, 2)
		return &arr
	}()
	isame := func() *[]bool {
		arr := make([]bool, 13)
		return &arr
	}()
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Cher2 and Chpr2.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*half) = (0.5 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Functions ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Data statements ..
	(*ich)[0], (*ich)[1] = 'u', 'l'
	//*     .. Executable Statements ..
	(*full) = (*sname)[3     - 1] == "e"
	(*packed) = (*sname)[3     - 1] == "p"
	//*     Define the number of arguments.
	if *full {
		(*nargs) = 9
	} else if *packed {
		(*nargs) = 8
	}
	//*
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*        Set LDA to 1 more than minimum value if room.
		(*lda) = (*n)
		if (*lda) < (*nmax) {
			(*lda) = (*lda) + 1
		}
		//*        Skip tests if not enough room.
		if (*lda) > (*nmax) {
			goto Label140
		}
		if *packed {
			(*laa) = ((*n) * ((*n) + 1)) / 2
		} else {
			(*laa) = (*lda) * (*n)
		}
		//*
		for (*ic) = 1; (*ic) <= 2; (*ic)++ {
			(*uplo) = (*ich)[(*ic)-1]
			(*upper) = (*uplo) == "u"
			//*
			for (*ix) = 1; (*ix) <= (*ninc); (*ix)++ {
				(*incx) = (*inc)[(*ix)-1]
				(*lx) = ABS((*incx)) * (*n)
				//*
				//*              Generate the vector X.
				//*
				(*transl) = (*half)
				Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, x, func() *int{y := 1; return &y}(), xx, ABS((*incx)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
				if (*n) > 1 {
					(*x)[(*n)/1] = (*zero)
					(*xx)[1+ABS((*incx))*((*n)/2-1)-1] = (*zero)
				}
				//*
				for (*iy) = 1; (*iy) <= (*ninc); (*iy)++ {
					(*incy) = (*inc)[(*iy)-1]
					(*ly) = ABS((*incy)) * (*n)
					//*
					//*                 Generate the vector Y.
					//*
					(*transl) = (*zero)
					Cmake(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *int{y := 1; return &y}(), n, y, func() *int{y := 1; return &y}(), yy, ABS((*incy)), func() *int{y := 0; return &y}(), (*n)-1, reset, transl)
					if (*n) > 1 {
						(*y)[(*n)/1] = (*zero)
						(*yy)[1+ABS((*incy))*((*n)/2-1)-1] = (*zero)
					}
					//*
					for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
						(*alpha) = (*alf)[(*ia)-1]
						(*null) = (*n) <= 0 || (*alpha) == (*zero)
						//*
						//*                    Generate the matrix A.
						//*
						(*transl) = (*zero)
						Cmake(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), n, n, a, nmax, aa, lda, (*n)-1, (*n)-1, reset, transl)
						//*
						(*nc) = (*nc) + 1
						//*
						//*                    Save every datum before calling the subroutine.
						//*
						(*uplos) = (*uplo)
						(*ns) = (*n)
						(*als) = (*alpha)
						for (*i) = 1; (*i) <= (*laa); (*i)++ {
							(*as)[(*i)-1] = (*aa)[(*i)-1]
						//Label10:
						}
						(*ldas) = (*lda)
						for (*i) = 1; (*i) <= (*lx); (*i)++ {
							(*xs)[(*i)-1] = (*xx)[(*i)-1]
						//Label20:
						}
						(*incxs) = (*incx)
						for (*i) = 1; (*i) <= (*ly); (*i)++ {
							(*ys)[(*i)-1] = (*yy)[(*i)-1]
						//Label30:
						}
						(*incys) = (*incy)
						//*
						//*                    Call the subroutine.
						//*
						if *full {
							if *trace {
								WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), x,%2d, y,%2d, a,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*incy), (*lda))
							}
							if *rewi {
								REWIND(ntra)
							}
							Cher2(uplo, n, alpha, xx, incx, yy, incy, aa, lda)
						} else if *packed {
							if *trace {
								WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), x,%2d, y,%2d, ap)                            .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*incy))
							}
							if *rewi {
								REWIND(ntra)
							}
							Chpr2(uplo, n, alpha, xx, incx, yy, incy, aa)
						}
						//*
						//*                    Check if error-exit was taken incorrectly.
						//*
						if !(*ok) {
							WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
							(*fatal) = true
							goto Label160
						}
						//*
						//*                    See what data changed inside subroutines.
						//*
						(*isame)[0] = (*uplo) == (*uplos)
						(*isame)[1] = (*ns) == (*n)
						(*isame)[2] = (*als) == (*alpha)
						(*isame)[3] = (*Lce(xs, xx, lx))
						(*isame)[4] = (*incxs) == (*incx)
						(*isame)[5] = (*Lce(ys, yy, ly))
						(*isame)[6] = (*incys) == (*incy)
						if *null {
							(*isame)[7] = (*Lce(as, aa, laa))
						} else {
							(*isame)[7] = (*Lceres(&((*sname)[1]), uplo, n, n, as, aa, lda))
						}
						if !(*packed) {
							(*isame)[8] = (*ldas) == (*lda)
						}
						//*
						//*                    If data was incorrectly changed, report and return.
						//*
						(*same) = true
						for (*i) = 1; (*i) <= (*nargs); (*i)++ {
							(*same) = (*same) && (*isame)[(*i)-1]
							if !(*isame)[(*i)-1] {
								WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - parameter number %2d was changed incorrectly *******\n"); return &y}(), (*i))
							}
						//Label40:
						}
						if !(*same) {
							(*fatal) = true
							goto Label160
						}
						//*
						if !(*null) {
							//*
							//*                       Check the result column by column.
							//*
							if (*incx) > 0 {
								for (*i) = 1; (*i) <= (*n); (*i)++ {
									(*z)[(*i)-1][0] = (*x)[(*i)-1]
								//Label50:
								}
							} else {
								for (*i) = 1; (*i) <= (*n); (*i)++ {
									(*z)[(*i)-1][0] = (*x)[(*n)-(*i)+0]
								//Label60:
								}
							}
							if (*incy) > 0 {
								for (*i) = 1; (*i) <= (*n); (*i)++ {
									(*z)[(*i)-1][1] = (*y)[(*i)-1]
								//Label70:
								}
							} else {
								for (*i) = 1; (*i) <= (*n); (*i)++ {
									(*z)[(*i)-1][1] = (*y)[(*n)-(*i)+0]
								//Label80:
								}
							}
							(*ja) = 1
							for (*j) = 1; (*j) <= (*n); (*j)++ {
								(*w)[0] = (*alpha) * CONJG(((*z)[(*j)-1][1]))
								(*w)[1] = CONJG((*alpha)) * CONJG(((*z)[(*j)-1][0]))
								if *upper {
									(*jj) = 1
									(*lj) = (*j)
								} else {
									(*jj) = (*j)
									(*lj) = (*n) - (*j) + 1
								}
								Cmvch( "n", (*lj), 2, (*one), (*z)[(*jj) - 1][0], (*nmax), (*w), 1, (*one), (*a)[(*jj) - 1][(*j) - 1], 1, (*yt), (*g), (*aa)[(*ja) - 1], (*eps), (*err), (*fatal), (*nout), true)
								if *full {
									if *upper {
										(*ja) = (*ja) + (*lda)
									} else {
										(*ja) = (*ja) + (*lda) + 1
									}
								} else {
									(*ja) = (*ja) + (*lj)
								}
								(*errmax) = (MAX((*errmax), (*err)))
								//*                          If got really bad answer, report and return.
								if *fatal {
									goto Label150
								}
							//Label90:
							}
						} else {
							//*                       Avoid repeating tests with N.le.0.
							if (*n) <= 0 {
								goto Label140
							}
						}
						//*
					//Label100:
					}
					//*
				//Label110:
				}
				//*
			//Label120:
			}
			//*
		//Label130:
		}
		//*
	Label140:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label170
	//*
Label150:
	;
	WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	//*
Label160:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *full {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), x,%2d, y,%2d, a,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*incy), (*lda))
	} else if *packed {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), x,%2d, y,%2d, ap)                            .\n"); return &y}(), (*nc), (*sname), (*uplo), (*n), (*alpha), (*incx), (*incy))
	}
	//*
Label170:
	;
	return
	//*
	//*
	//*     End of Cchk6.
	//*
}

func Cchke(isnum *int, srnamt *[]byte, nout *int) {
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	alpha := new(complex128)
	beta := new(complex128)
	ralpha := new(float64)
	a := func() *[][]complex128 {
		arr := make([][]complex128, 1)
		for u := 0; u < 1; u++ {
			arr[u] = make([]complex128, 1)
		}
		return &arr
	}()
	x := func() *[]complex128 {
		arr := make([]complex128, 1)
		return &arr
	}()
	y := func() *[]complex128 {
		arr := make([]complex128, 1)
		return &arr
	}()
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests the error exits from the Level 2 Blas.
	//*  Requires a special version of the error-handling routine Xerbla.
	//*  ALPHA, RALPHA, BETA, A, X and Y should not need to be defined.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Subroutines ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	noutc = common.infoc.noutc
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	//*     .. Executable Statements ..
	//*     ok is set to .FALSE. by the special version of Xerbla or by CHKXER
	//*     if anything is wrong.
	(*ok) = true
	//*     lerr is set to .TRUE. by the special version of Xerbla each time
	//*     it is called, and is then tested and re-set by CHKXER.
	(*lerr) = false
	switch *isnum {
	case 1:
		goto Label10
	case 2:
		goto Label20
	case 3:
		goto Label30
	case 4:
		goto Label40
	case 5:
		goto Label50
	case 6:
		goto Label60
	case 7:
		goto Label70
	case 8:
		goto Label80
	case 9:
		goto Label90
	case 10:
		goto Label100
	case 11:
		goto Label110
	case 12:
		goto Label120
	case 13:
		goto Label130
	case 14:
		goto Label140
	case 15:
		goto Label150
	case 16:
		goto Label160
	case 17:
		goto Label170
	}
Label10:
	;
	(*infot) = 1
	Cgemv(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cgemv(func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Cgemv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Cgemv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Cgemv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Cgemv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label20:
	;
	(*infot) = 1
	Cgbmv(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cgbmv(func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 1; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Cgbmv(func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label30:
	;
	(*infot) = 1
	Chemv(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Chemv(func() *byte{y := byte('u'); return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Chemv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Chemv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Chemv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label40:
	;
	(*infot) = 1
	Chbmv(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Chbmv(func() *byte{y := byte('u'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Chbmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Chbmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 1; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Chbmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Chbmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label50:
	;
	(*infot) = 1
	Chpmv(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Chpmv(func() *byte{y := byte('u'); return &y}(), -1, alpha, a, x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Chpmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, x, func() *int{y := 0; return &y}(), beta, y, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Chpmv(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, a, x, func() *int{y := 1; return &y}(), beta, y, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label60:
	;
	(*infot) = 1
	Ctrmv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctrmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctrmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctrmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ctrmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Ctrmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label70:
	;
	(*infot) = 1
	Ctbmv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ctbmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label80:
	;
	(*infot) = 1
	Ctpmv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctpmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctpmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctpmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Ctpmv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label90:
	;
	(*infot) = 1
	Ctrsv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctrsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctrsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctrsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ctrsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Ctrsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label100:
	;
	(*infot) = 1
	Ctbsv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ctbsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}(), x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label110:
	;
	(*infot) = 1
	Ctpsv(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ctpsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ctpsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ctpsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, a, x, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Ctpsv(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), a, x, func() *int{y := 0; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label120:
	;
	(*infot) = 1
	Cgerc(-1, func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cgerc(func() *int{y := 0; return &y}(), -1, alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Cgerc(func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 0; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Cgerc(func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Cgerc(func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label130:
	;
	(*infot) = 1
	Cgeru(-1, func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cgeru(func() *int{y := 0; return &y}(), -1, alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Cgeru(func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 0; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Cgeru(func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Cgeru(func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label140:
	;
	(*infot) = 1
	Cher(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), ralpha, x, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cher(func() *byte{y := byte('u'); return &y}(), -1, ralpha, x, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Cher(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), ralpha, x, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Cher(func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), ralpha, x, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label150:
	;
	(*infot) = 1
	Chpr(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), ralpha, x, func() *int{y := 1; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Chpr(func() *byte{y := byte('u'); return &y}(), -1, ralpha, x, func() *int{y := 1; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Chpr(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), ralpha, x, func() *int{y := 0; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label160:
	;
	(*infot) = 1
	Cher2(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Cher2(func() *byte{y := byte('u'); return &y}(), -1, alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Cher2(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 0; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Cher2(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 0; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Cher2(func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label180
Label170:
	;
	(*infot) = 1
	Chpr2(func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Chpr2(func() *byte{y := byte('u'); return &y}(), -1, alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 1; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Chpr2(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 0; return &y}(), y, func() *int{y := 1; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Chpr2(func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), alpha, x, func() *int{y := 1; return &y}(), y, func() *int{y := 0; return &y}(), a)
	CHKXER(srnamt, infot, nout, lerr, ok)
	//*
Label180:
	;
	if *ok {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the tests of error-exits\n"); return &y}(), (*srnamt))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed the tests of error-exits *******\n"); return &y}(), (*srnamt))
	}
	return
	//*
	//*
	//*     End of Cchke.
	//*
}

func Cmake(_type *int, uplo *byte, diag *byte, m *int, n *int, a *[][]complex128, nmax *int, aa *[]complex128, lda *int, kl *int, ku *int, reset *bool, transl *complex128) {
	zero := new(complex128)
	one := new(complex128)
	rogue := new(complex128)
	rzero := new(float64)
	rrogue := new(float64)
	type := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	i := new(int)
	i1 := new(int)
	i2 := new(int)
	i3 := new(int)
	ibeg := new(int)
	iend := new(int)
	ioff := new(int)
	j := new(int)
	jj := new(int)
	kk := new(int)
	gen := new(bool)
	lower := new(bool)
	sym := new(bool)
	tri := new(bool)
	unit := new(bool)
	upper := new(bool)
	//*
	//*  Generates values for an M by N matrix A within the bandwidth
	//*  defined by KL and KU.
	//*  Stores the values in the array AA in the data structure required
	//*  by the routine, with unwanted elements set to rogue value.
	//*
	//*  TYPE is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rogue) = (-1.0e10 + (1.0e10)*1i)
	(*rzero) = 0.0
	(*rrogue) = -1.0e10
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. External Functions ..
	//*     .. Intrinsic Functions ..
	//*     .. Executable Statements ..
	(*gen) = (*type)[1     - 1] == "g"
	(*sym) = (*type)[1     - 1] == "h"
	(*tri) = (*type)[1     - 1] == "t"
	(*upper) = ((*sym) || (*tri)) && (*uplo) == "u"
	(*lower) = ((*sym) || (*tri)) && (*uplo) == "l"
	(*unit) = (*tri) && (*diag) == "u"
	//*
	//*     Generate data in array A.
	//*
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		for (*i) = 1; (*i) <= (*m); (*i)++ {
			if (*gen) || ((*upper) && (*i) <= (*j)) || ((*lower) && (*i) >= (*j)) {
				if ((*i) <= (*j) && (*j) - (*i) <= (*ku)) || ((*i) >= (*j) && (*i) - (*j) <= (*kl)) {
					(*a)[(*i)-1][(*j)-1] = CBEG(reset) + (*transl)
				} else {
					(*a)[(*i)-1][(*j)-1] = (*zero)
				}
				if (*i) != (*j) {
					if *sym {
						(*a)[(*j)-1][(*i)-1] = (CONJG(((*a)[(*i)-1][(*j)-1])))
					} else if *tri {
						(*a)[(*j)-1][(*i)-1] = (*zero)
					}
				}
			}
		//Label10:
		}
		if *sym {
			(*a)[(*j)-1][(*j)-1] = (*cmplx(real(((*a)[(*j)-1][(*j)-1])), rzero))
		}
		if *tri {
			(*a)[(*j)-1][(*j)-1] = (*a)[(*j)-1][(*j)-1] + (*one)
		}
		if *unit {
			(*a)[(*j)-1][(*j)-1] = (*one)
		}
	//Label20:
	}
	//*
	//*     Store elements in array AS in data structure required by routine.
	//*
	if (*type) == "ge" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			for (*i) = 1; (*i) <= (*m); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*a)[(*i)-1][(*j)-1]
			//Label30:
			}
			for (*i) = (*m) + 1; (*i) <= (*lda); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label40:
			}
		//Label50:
		}
	} else if (*type) == "gb" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			for (*i1) = 1; (*i1) <= (*ku)+1-(*j); (*i1)++ {
				(*aa)[(*i1)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label60:
			}
			for (*i2) = (*i1); (*i2) <= (MIN((*kl)+(*ku)+1, (*ku)+1+(*m)-(*j))); (*i2)++ {
				(*aa)[(*i2)+((*j)-1)*(*lda)-1] = (*a)[(*i2)+(*j)-(*ku)-0][(*j)-1]
			//Label70:
			}
			for (*i3) = (*i2); (*i3) <= (*lda); (*i3)++ {
				(*aa)[(*i3)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label80:
			}
		//Label90:
		}
	} else if (*type) == "he" || (*type) == "tr" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			if *upper {
				(*ibeg) = 1
				if *unit {
					(*iend) = (*j) - 1
				} else {
					(*iend) = (*j)
				}
			} else {
				if *unit {
					(*ibeg) = (*j) + 1
				} else {
					(*ibeg) = (*j)
				}
				(*iend) = (*n)
			}
			for (*i) = 1; (*i) <= (*ibeg)-1; (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label100:
			}
			for (*i) = (*ibeg); (*i) <= (*iend); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*a)[(*i)-1][(*j)-1]
			//Label110:
			}
			for (*i) = (*iend) + 1; (*i) <= (*lda); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label120:
			}
			if *sym {
				(*jj) = (*j) + ((*j)-1)*(*lda)
				(*aa)[(*jj)-1] = (*cmplx(real(((*aa)[(*jj)-1])), rrogue))
			}
		//Label130:
		}
	} else if (*type) == "hb" || (*type) == "tb" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			if *upper {
				(*kk) = (*kl) + 1
				(*ibeg) = (MAX(1, (*kl)+2-(*j)))
				if *unit {
					(*iend) = (*kl)
				} else {
					(*iend) = (*kl) + 1
				}
			} else {
				(*kk) = 1
				if *unit {
					(*ibeg) = 2
				} else {
					(*ibeg) = 1
				}
				(*iend) = (MIN((*kl)+1, 1+(*m)-(*j)))
			}
			for (*i) = 1; (*i) <= (*ibeg)-1; (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label140:
			}
			for (*i) = (*ibeg); (*i) <= (*iend); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*a)[(*i)+(*j)-(*kk)-1][(*j)-1]
			//Label150:
			}
			for (*i) = (*iend) + 1; (*i) <= (*lda); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label160:
			}
			if *sym {
				(*jj) = (*kk) + ((*j)-1)*(*lda)
				(*aa)[(*jj)-1] = (*cmplx(real(((*aa)[(*jj)-1])), rrogue))
			}
		//Label170:
		}
	} else if (*type) == "hp" || (*type) == "tp" {
		(*ioff) = 0
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			if *upper {
				(*ibeg) = 1
				(*iend) = (*j)
			} else {
				(*ibeg) = (*j)
				(*iend) = (*n)
			}
			for (*i) = (*ibeg); (*i) <= (*iend); (*i)++ {
				(*ioff) = (*ioff) + 1
				(*aa)[(*ioff)-1] = (*a)[(*i)-1][(*j)-1]
				if (*i) == (*j) {
					if *unit {
						(*aa)[(*ioff)-1] = (*rogue)
					}
					if *sym {
						(*aa)[(*ioff)-1] = (*cmplx(real(((*aa)[(*ioff)-1])), rrogue))
					}
				}
			//Label180:
			}
		//Label190:
		}
	}
	return
	//*
	//*     End of Cmake.
	//*
}

func Cmvch(trans *byte, m *int, n *int, alpha *complex128, a *[][]complex128, nmax *int, x *[]complex128, incx *int, beta *complex128, y *[]complex128, incy *int, yt *[]complex128, g *[]float64, yy *[]complex128, eps *float64, err *float64, fatal *bool, nout *int, mv *bool) {
	zero := new(complex128)
	rzero := new(float64)
	rone := new(float64)
	c := new(complex128)
	erri := new(float64)
	i := new(int)
	incxl := new(int)
	incyl := new(int)
	iy := new(int)
	j := new(int)
	jx := new(int)
	kx := new(int)
	ky := new(int)
	ml := new(int)
	nl := new(int)
	ctran := new(bool)
	tran := new(bool)
	abs1 := new(float64)
	//*
	//*  Checks the results of the computational tests.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*rzero) = 0.0
	(*rone) = 1.0
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Intrinsic Functions ..
	//*     .. Statement Functions ..
	//*     .. Statement Function definitions ..
	ABS1(c) = ABS(real((*c))) + ABS(imag((*c)))
	//*     .. Executable Statements ..
	(*tran) = (*trans) == "t"
	(*ctran) = (*trans) == "c"
	if (*tran) || (*ctran) {
		(*ml) = (*n)
		(*nl) = (*m)
	} else {
		(*ml) = (*m)
		(*nl) = (*n)
	}
	if (*incx) < 0 {
		(*kx) = (*nl)
		(*incxl) = -1
	} else {
		(*kx) = 1
		(*incxl) = 1
	}
	if (*incy) < 0 {
		(*ky) = (*ml)
		(*incyl) = -1
	} else {
		(*ky) = 1
		(*incyl) = 1
	}
	//*
	//*     Compute expected result in YT using data in A, X and Y.
	//*     Compute gauges in G.
	//*
	(*iy) = (*ky)
	for (*i) = 1; (*i) <= (*ml); (*i)++ {
		(*yt)[(*iy)-1] = (*zero)
		(*g)[(*iy)-1] = (*rzero)
		(*jx) = (*kx)
		if *tran {
			for (*j) = 1; (*j) <= (*nl); (*j)++ {
				(*yt)[(*iy)-1] = (*yt)[(*iy)-1] + (*a)[(*j)-1][(*i)-1]*(*x)[(*jx)-1]
				(*g)[(*iy)-1] = (*g)[(*iy)-1] + ABS1(&((*a)[(*j)-1][(*i)-1]))*ABS1(&((*x)[(*jx)-1]))
				(*jx) = (*jx) + (*incxl)
			//Label10:
			}
		} else if *ctran {
			for (*j) = 1; (*j) <= (*nl); (*j)++ {
				(*yt)[(*iy)-1] = (*yt)[(*iy)-1] + CONJG(((*a)[(*j)-1][(*i)-1]))*(*x)[(*jx)-1]
				(*g)[(*iy)-1] = (*g)[(*iy)-1] + ABS1(&((*a)[(*j)-1][(*i)-1]))*ABS1(&((*x)[(*jx)-1]))
				(*jx) = (*jx) + (*incxl)
			//Label20:
			}
		} else {
			for (*j) = 1; (*j) <= (*nl); (*j)++ {
				(*yt)[(*iy)-1] = (*yt)[(*iy)-1] + (*a)[(*i)-1][(*j)-1]*(*x)[(*jx)-1]
				(*g)[(*iy)-1] = (*g)[(*iy)-1] + ABS1(&((*a)[(*i)-1][(*j)-1]))*ABS1(&((*x)[(*jx)-1]))
				(*jx) = (*jx) + (*incxl)
			//Label30:
			}
		}
		(*yt)[(*iy)-1] = (*alpha)*(*yt)[(*iy)-1] + (*beta)*(*y)[(*iy)-1]
		(*g)[(*iy)-1] = ABS1(alpha)*(*g)[(*iy)-1] + ABS1(beta)*ABS1(&((*y)[(*iy)-1]))
		(*iy) = (*iy) + (*incyl)
	//Label40:
	}
	//*
	//*     Compute the error ratio for this result.
	//*
	(*err) = (*zero)
	for (*i) = 1; (*i) <= (*ml); (*i)++ {
		(*erri) = ABS((*yt)[(*i)-1]-(*yy)[1+((*i)-1)*ABS((*incy))-1]) / (*eps)
		if (*g)[(*i)-1] != (*rzero) {
			(*erri) = (*erri) / (*g)[(*i)-1]
		}
		(*err) = (MAX((*err), (*erri)))
		if (*err) * SQRT((*eps)) >= (*rone) {
			goto Label60
		}
	//Label50:
	}
	//*     If the loop completes, all results are at least half accurate.
	goto Label80
	//*
	//*     Report fatal error.
	//*
Label60:
	;
	(*fatal) = true
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - computed result is less than half accurate *******\n                       expected result                    computed result\n"); return &y}())
	for (*i) = 1; (*i) <= (*ml); (*i)++ {
		if *mv {
			WRITE((*nout), *func() *[]byte{y := []byte(" %7d  (%15.6f,%15.6f)\n"); return &y}(), (*i), (*yt)[(*i)-1], (*yy)[1+((*i)-1)*ABS((*incy))-1])
		} else {
			WRITE((*nout), *func() *[]byte{y := []byte(" %7d  (%15.6f,%15.6f)\n"); return &y}(), (*i), (*yy)[1+((*i)-1)*ABS((*incy))-1], (*yt)[(*i)-1])
		}
	//Label70:
	}
	//*
Label80:
	;
	return
	//*
	//*
	//*     End of CMVCH.
	//*
}

func Lce(ri *[]complex128, rj *[]complex128, lr *int) (lceReturn *bool) {
	lcereturn := new(bool)
	i := new(int)
	//*
	//*  Tests if two arrays are identical.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Executable Statements ..
	for (*i) = 1; (*i) <= (*lr); (*i)++ {
		if (*ri)[(*i)-1] != (*rj)[(*i)-1] {
			goto Label20
		}
	//Label10:
	}
	(*lce) = true
	goto Label30
Label20:
	;
	(*lce) = false
Label30:
	;
	return
	//*
	//*     End of Lce.
	//*
}

func Lceres(_type *int, uplo *byte, m *int, n *int, aa *[][]complex128, as *[][]complex128, lda *int) (lceresReturn *bool) {
	lceresreturn := new(bool)
	type := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	i := new(int)
	ibeg := new(int)
	iend := new(int)
	j := new(int)
	upper := new(bool)
	//*
	//*  Tests if selected elements in two arrays are equal.
	//*
	//*  TYPE is 'GE', 'HE' or 'HP'.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Array Arguments ..
	//*     .. Local Scalars ..
	//*     .. Executable Statements ..
	(*upper) = (*uplo) == "u"
	if (*type) == "ge" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			for (*i) = (*m) + 1; (*i) <= (*lda); (*i)++ {
				if (*aa)[(*i)-1][(*j)-1] != (*as)[(*i)-1][(*j)-1] {
					goto Label70
				}
			//Label10:
			}
		//Label20:
		}
	} else if (*type) == "he" {
		for (*j) = 1; (*j) <= (*n); (*j)++ {
			if *upper {
				(*ibeg) = 1
				(*iend) = (*j)
			} else {
				(*ibeg) = (*j)
				(*iend) = (*n)
			}
			for (*i) = 1; (*i) <= (*ibeg)-1; (*i)++ {
				if (*aa)[(*i)-1][(*j)-1] != (*as)[(*i)-1][(*j)-1] {
					goto Label70
				}
			//Label30:
			}
			for (*i) = (*iend) + 1; (*i) <= (*lda); (*i)++ {
				if (*aa)[(*i)-1][(*j)-1] != (*as)[(*i)-1][(*j)-1] {
					goto Label70
				}
			//Label40:
			}
		//Label50:
		}
	}
	//*
	(*lceres) = true
	goto Label80
Label70:
	;
	(*lceres) = false
Label80:
	;
	return
	//*
	//*     End of Lceres.
	//*
}

func Cbeg(reset *bool) (cbegReturn *complex128) {
	cbegreturn := new(complex128)
	i := new(int)
	ic := new(int)
	j := new(int)
	mi := new(int)
	mj := new(int)
	//*
	//*  Generates complex numbers as pairs of random numbers uniformly
	//*  distributed between -0.5 and 0.5.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Local Scalars ..
	//*     .. Save statement ..
	//*     .. Intrinsic Functions ..
	//*     .. Executable Statements ..
	if *reset {
		//*        Initialize local variables.
		(*mi) = 891
		(*mj) = 457
		(*i) = 7
		(*j) = 7
		(*ic) = 0
		(*reset) = false
	}
	//*
	//*     The sequence of values of I or J is bounded between 1 and 999.
	//*     If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
	//*     If initial I or J = 4 or 8, the period will be 25.
	//*     If initial I or J = 5, the period will be 10.
	//*     IC is used to break up the period by skipping 1 value of I or J
	//*     in 6.
	//*
	(*ic) = (*ic) + 1
Label10:
	;
	(*i) = (*i) * (*mi)
	(*j) = (*j) * (*mj)
	(*i) = (*i) - 1000*((*i)/1000)
	(*j) = (*j) - 1000*((*j)/1000)
	if (*ic) >= 5 {
		(*ic) = 0
		goto Label10
	}
	(*cbeg) = (*cmplx(((*i)-500)/1001.0, ((*j)-500)/1001.0))
	return
	//*
	//*     End of CBEG.
	//*
}

func Sdiff(x *float64, y *float64) (sdiffReturn *float64) {
	sdiffreturn := new(float64)
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Executable Statements ..
	(*sdiff) = (*x) - (*y)
	return
	//*
	//*     End of SDIFF.
	//*
}

func Chkxer(srnamt *[]byte, infot *int, nout *int, lerr *bool, ok *bool) {
	//*
	//*  Tests whether Xerbla has detected an error when it should.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Executable Statements ..
	if !(*lerr) {
		WRITE((*nout), *func() *[]byte{y := []byte(" ***** illegal value of parameter number %2d not detected by %6s *****\n"); return &y}(), (*infot), (*srnamt))
		(*ok) = false
	}
	(*lerr) = false
	return
	//*
	//*
	//*     End of CHKXER.
	//*
}

func Xerbla(srname *[]byte, info *int) {
	infot := new(int)
	nout := new(int)
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 6)
		return &arr
	}()
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.nout := new(float64)
	common.infoc.infot := new(int)
	common.srnamc.srnamt := new(int)
	//*
	//*  This is a special version of Xerbla to be used only as part of
	//*  the test program for testing error exits from the Level 2 BLAS
	//*  routines.
	//*
	//*  Xerbla  is an error handler for the Level 2 BLAS routines.
	//*
	//*  It is called by the Level 2 BLAS routines if an input parameter is
	//*  invalid.
	//*
	//*  Auxiliary routine for test program for Level 2 Blas.
	//*
	//*  -- Written on 10-August-1987.
	//*     Richard Hanson, Sandia National Labs.
	//*     Jeremy Du Croz, NAG Central Office.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Common blocks ..
	infot = common.infoc.infot
	nout = common.infoc.nout
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	//*     .. Executable Statements ..
	(*lerr) = true
	if (*info) != (*infot) {
		if (*infot) != 0 {
			WRITE((*nout), *func() *[]byte{y := []byte(" ******* xerbla was called with info = %6d instead of %2d *******\n"); return &y}(), (*info), (*infot))
		} else {
			WRITE((*nout), *func() *[]byte{y := []byte(" ******* xerbla was called with info = %6d *******\n"); return &y}(), (*info))
		}
		(*ok) = false
	}
	if (*srname) != (*srnamt) {
		WRITE((*nout), *func() *[]byte{y := []byte(" ******* xerbla was called with srname = %6s instead of %6s *******\n"); return &y}(), (*srname), (*srnamt))
		(*ok) = false
	}
	return
	//*
	//*
	//*     End of Xerbla
	//*
}
