package goblas

import 

// \brief \b Zblat3
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       PROGRAM Zblat3
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Test program for the COMPLEX//16       Level 3 Blas.
//
// The program must be driven by a short data file. The first 14 records
// of the file are read using list-directed input, the last 9 records
// are read using the format ( A6, L2). An annotated example of a data
// file can be obtained by deleting the first 3 characters from the
// following 23 lines:
// 'zblat3.out'      NAME OF SUMMARY OUTPUT FILE
// 6                 UNIT NUMBER OF SUMMARY FILE
// 'Zblat3.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
// -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
// F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
// F        LOGICAL FLAG, T TO STOP ON FAILURES.
// T        LOGICAL FLAG, T TO TEST ERROR EXITS.
// 16.0     THRESHOLD VALUE OF TEST RATIO
// 6                 NUMBER OF VALUES OF N
// 0 1 2 3 5 9       VALUES OF N
// 3                 NUMBER OF VALUES OF ALPHA
// (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
// 3                 NUMBER OF VALUES OF BETA
// (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
// Zgemm  T PUT F FOR NO TEST. SAME COLUMNS.
// Zhemm  T PUT F FOR NO TEST. SAME COLUMNS.
// Zsymm  T PUT F FOR NO TEST. SAME COLUMNS.
// Ztrmm  T PUT F FOR NO TEST. SAME COLUMNS.
// Ztrsm  T PUT F FOR NO TEST. SAME COLUMNS.
// Zherk  T PUT F FOR NO TEST. SAME COLUMNS.
// Zsyrk  T PUT F FOR NO TEST. SAME COLUMNS.
// Zher2k T PUT F FOR NO TEST. SAME COLUMNS.
// Zsyr2k T PUT F FOR NO TEST. SAME COLUMNS.
//
//
// Further Details
// ===============
//
// See:
//
//    Dongarra J. J., Du Croz J. J., Duff I. S. and Hammarling S.
//    A Set of Level 3 Basic Linear Algebra Subprograms.
//
//    Technical Memorandum No.88 (Revision 1), Mathematics and
//    Computer Science Division, Argonne National Laboratory, 9700
//    South Cass Avenue, Argonne, Illinois 60439, US.
//
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.
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
// \ingroup complex16_blas_testing
//
//  =====================================================================
func main() {
	nin := new(int)
	nsubs := new(int)
	zero := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	nmax := new(int)
	nidmax := new(int)
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
	nout := new(int)
	ntra := new(int)
	fatal := new(bool)
	ltestt := new(bool)
	rewi := new(bool)
	same := new(bool)
	sfatal := new(bool)
	trace := new(bool)
	tsterr := new(bool)
	transa := new(byte)
	transb := new(byte)
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
	aa := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	ab := func() *[][]complex128 {
		arr := make([][]complex128, -1)
		for u := 0; u < -1; u++ {
			arr[u] = make([]complex128, -1)
		}
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
	bb := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	bet := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	bs := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	c := func() *[][]complex128 {
		arr := make([][]complex128, -1)
		for u := 0; u < -1; u++ {
			arr[u] = make([]complex128, -1)
		}
		return &arr
	}()
	cc := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	cs := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	ct := func() *[]complex128 {
		arr := make([]complex128, -1)
		return &arr
	}()
	w := func() *[]complex128 {
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
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	common.srnamc.srnamt := new(int)
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
	(*nsubs) = 9
	(*zero) = (0.0 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rzero) = 0.0
	(*nmax) = 65
	(*nidmax) = 9
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
	//F4GO: NOT IMPLEMENTED :data snames / "zgemm ", "zhemm ", "zsymm ", "ztrmm ", "ztrsm ", "zherk ", "zsyrk ", "zher2k", "zsyr2k" /
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
		OPEN(NTRA, snaps)
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
	if (*nidim) < 1|| (*nidim) > (*nidmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), 'n', (*nidmax))
		goto Label220
	}
	WRITE((*nin), " %v\n", ((*idim)[(*i)-1], (*i) = 1, (*nidim)))
	for (*i) = 1; (*i) <= (*nidim); (*i)++ {
		if (*idim)[(*i)-1] < 0|| (*idim)[(*i)-1] > (*nmax) {
			WRITE((*nout), *func() *[]byte{y := []byte(" value of n is less than 0 or greater than %2d\n"); return &y}(), (*nmax))
			goto Label220
		}
	//Label10:
	}
	//*     Values of ALPHA
	READ(nin, []byte(" %v\n"), nalf)
	if (*nalf) < 1|| (*nalf) > (*nalmax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), *func() *[]byte{y := []byte("alpha"); return &y}(), (*nalmax))
		goto Label220
	}
	WRITE((*nin), " %v\n", ((*alf)[(*i)-1], (*i) = 1, (*nalf)))
	//*     Values of BETA
	READ(nin, []byte(" %v\n"), nbet)
	if (*nbet) < 1|| (*nbet) > (*nbemax) {
		WRITE((*nout), *func() *[]byte{y := []byte(" number of values of %s is less than 1 or greater than %2d\n"); return &y}(), *func() *[]byte{y := []byte("beta"); return &y}(), (*nbemax))
		goto Label220
	}
	WRITE((*nin), " %v\n", ((*bet)[(*i)-1], (*i) = 1, (*nbet)))
	//*
	//*     Report values of parameters.
	//*
	WRITE((*nout), *func() *[]byte{y := []byte(" tests of the complex*16       level 3 blas%v the following parameter values will be used:\n"); return &y}())
	WRITE((*nout), "   for n              %6d\n", ((*idim)[(*i)-1], (*i) = 1, (*nidim)))
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
	//Label20:
	}
Label30:
	;
	READ(nin, []byte("%6s  %t\n"), snamet, ltestt)
	for (*i) = 1; (*i) <= (*nsubs); (*i)++ {
		if (*snamet) == (*snames)[(*i)-1] {
			goto Label50
		}
	//Label40:
	}
	WRITE((*nout), *func() *[]byte{y := []byte(" subprogram name %6s not recognized\n ******* tests abandoned *******\n"); return &y}(), (*snamet))
	panic("")
Label50:
	;
	(*ltest)[(*i)-1] = (*ltestt)
	goto Label30
	//*
//Label60:
	;
	CLOSE(NIN)
	//*
	//*     Compute EPS (the machine precision).
	//*
	(*eps) = (EPSILON((*rzero)))
	WRITE((*nout), *func() *[]byte{y := []byte(" relative machine precision is taken to be%f%9.1e\n"); return &y}(), (*eps))
	//*
	//*     Check the reliability of ZMMCH using exact data.
	//*
	(*n) = (MIN(int(32), (*nmax)))
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		for (*i) = 1; (*i) <= (*n); (*i)++ {
			(*ab)[(*i)-1][(*j)-1] = (MAX((*i)-(*j)+1, 0))
		//Label90:
		}
		(*ab)[(*j)-1][(*nmax)+0] = (*j)
		(*ab)[0][(*nmax)+(*j)-1] = (*j)
		(*c)[(*j)-1][0] = (*zero)
	//Label100:
	}
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		(*cc)[(*j)-1] = (*j)*(((*j)+1)*(*j))/2 - (((*j)+1)*(*j)*((*j)-1))/3
	//Label110:
	}
	//*     CC holds the exact result. On exit from ZMMCH CT holds
	//*     the result computed by ZMMCH.
	(*transa) = 'n'
	(*transb) = 'n'
	zmmch ((*transa), (*transb), (*n), 1, (*n), (*one), (*ab), (*nmax), (*ab)[0][(*nmax) + 0], (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*cc), (*nmax), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*LZE(cc, ct, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in zmmch -  in-line dot products are being evaluated wrongly.\n zmmch was called with transa = %c and transb = %c\n and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*transa), (*transb), (*same), (*err))
		panic("")
	}
	(*transb) = 'c'
	zmmch ((*transa), (*transb), (*n), 1, (*n), (*one), (*ab), (*nmax), (*ab)[0][(*nmax) + 0], (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*cc), (*nmax), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*LZE(cc, ct, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in zmmch -  in-line dot products are being evaluated wrongly.\n zmmch was called with transa = %c and transb = %c\n and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*transa), (*transb), (*same), (*err))
		panic("")
	}
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		(*ab)[(*j)-1][(*nmax)+0] = (*n) - (*j) + 1
		(*ab)[0][(*nmax)+(*j)-1] = (*n) - (*j) + 1
	//Label120:
	}
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		(*cc)[(*n)-(*j)+0] = (*j)*(((*j)+1)*(*j))/2 - (((*j)+1)*(*j)*((*j)-1))/3
	//Label130:
	}
	(*transa) = 'c'
	(*transb) = 'n'
	zmmch ((*transa), (*transb), (*n), 1, (*n), (*one), (*ab), (*nmax), (*ab)[0][(*nmax) + 0], (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*cc), (*nmax), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*LZE(cc, ct, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in zmmch -  in-line dot products are being evaluated wrongly.\n zmmch was called with transa = %c and transb = %c\n and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*transa), (*transb), (*same), (*err))
		panic("")
	}
	(*transb) = 'c'
	zmmch ((*transa), (*transb), (*n), 1, (*n), (*one), (*ab), (*nmax), (*ab)[0][(*nmax) + 0], (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*cc), (*nmax), (*eps), (*err), (*fatal), (*nout), true)
	(*same) = (*LZE(cc, ct, n))
	if !(*same) || (*err) != (*rzero) {
		WRITE((*nout), *func() *[]byte{y := []byte(" error in zmmch -  in-line dot products are being evaluated wrongly.\n zmmch was called with transa = %c and transb = %c\n and returned same =  %t and err = %12.3f.\n this may be due to faults in the arithmetic or the compiler.\n ******* tests abandoned *******\n"); return &y}(), (*transa), (*transb), (*same), (*err))
		panic("")
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
				ZCHKE(isnum, &((*snames)[(*isnum)-1]), nout)
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
				goto Label150
			case 3:
				goto Label150
			case 4:
				goto Label160
			case 5:
				goto Label160
			case 6:
				goto Label170
			case 7:
				goto Label170
			case 8:
				goto Label180
			case 9:
				goto Label180
			}
			//*           Test Zgemm, 01.
		Label140:
			;
			ZCHK1(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, nbet, bet, nmax, ab, aa, as, &((*ab)[0][(*nmax)+0]), bb, bs, c, cc, cs, ct, g)
			goto Label190
			//*           Test Zhemm, 02, Zsymm, 03.
		Label150:
			;
			ZCHK2(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, nbet, bet, nmax, ab, aa, as, &((*ab)[0][(*nmax)+0]), bb, bs, c, cc, cs, ct, g)
			goto Label190
			//*           Test Ztrmm, 04, Ztrsm, 05.
		Label160:
			;
			ZCHK3(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, nmax, ab, aa, as, &((*ab)[0][(*nmax)+0]), bb, bs, ct, g, c)
			goto Label190
			//*           Test Zherk, 06, Zsyrk, 07.
		Label170:
			;
			ZCHK4(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, nbet, bet, nmax, ab, aa, as, &((*ab)[0][(*nmax)+0]), bb, bs, c, cc, cs, ct, g)
			goto Label190
			//*           Test Zher2k, 08, Zsyr2k, 09.
		Label180:
			;
			ZCHK5(&((*snames)[(*isnum)-1]), eps, thresh, nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, nbet, bet, nmax, ab, aa, as, bb, bs, c, cc, cs, ct, g, w)
			goto Label190
			//*
		Label190:
			;
			if (*fatal) && (*sfatal) {
				goto Label210
			}
		}
	//Label200:
	}
	WRITE((*nout), *func() *[]byte{y := []byte("\n end of tests\n"); return &y}())
	goto Label230
	//*
Label210:
	;
	WRITE((*nout), *func() *[]byte{y := []byte("\n ******* fatal error - tests abandoned *******\n"); return &y}())
	goto Label230
	//*
Label220:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" amend data file or increase array sizes in program\n ******* tests abandoned *******\n"); return &y}())
	//*
Label230:
	;
	if *trace {
		CLOSE(NTRA)
	}
	CLOSE(nout)
	panic("")
	//*
	//*
	//*     End of Zblat3.
	//*
}

func Zchk1(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, nmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, b *[][]complex128, bb *[]complex128, bs *[]complex128, c *[][]complex128, cc *[]complex128, cs *[]complex128, ct *[]complex128, g *[]float64) {
	zero := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bls := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ica := new(int)
	icb := new(int)
	ik := new(int)
	im := new(int)
	in := new(int)
	k := new(int)
	ks := new(int)
	laa := new(int)
	lbb := new(int)
	lcc := new(int)
	lda := new(int)
	ldas := new(int)
	ldb := new(int)
	ldbs := new(int)
	ldc := new(int)
	ldcs := new(int)
	m := new(int)
	ma := new(int)
	mb := new(int)
	ms := new(int)
	n := new(int)
	na := new(int)
	nargs := new(int)
	nb := new(int)
	nc := new(int)
	ns := new(int)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	trana := new(bool)
	tranb := new(bool)
	tranas := new(byte)
	tranbs := new(byte)
	transa := new(byte)
	transb := new(byte)
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
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests Zgemm.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
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
	//*
	(*nargs) = 13
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*im) = 1; (*im) <= (*nidim); (*im)++ {
		(*m) = (*idim)[(*im)-1]
		//*
		for (*in) = 1; (*in) <= (*nidim); (*in)++ {
			(*n) = (*idim)[(*in)-1]
			//*           Set LDC to 1 more than minimum value if room.
			(*ldc) = (*m)
			if (*ldc) < (*nmax) {
				(*ldc) = (*ldc) + 1
			}
			//*           Skip tests if not enough room.
			if (*ldc) > (*nmax) {
				goto Label100
			}
			(*lcc) = (*ldc) * (*n)
			(*null) = (*n) <= 0|| (*m) <= 0
			//*
			for (*ik) = 1; (*ik) <= (*nidim); (*ik)++ {
				(*k) = (*idim)[(*ik)-1]
				//*
				for (*ica) = 1; (*ica) <= 3; (*ica)++ {
					(*transa) = (*ich)[(*ica)-1]
					(*trana) = (*transa) == "t" || (*transa) == "c"
					//*
					if *trana {
						(*ma) = (*k)
						(*na) = (*m)
					} else {
						(*ma) = (*m)
						(*na) = (*k)
					}
					//*                 Set LDA to 1 more than minimum value if room.
					(*lda) = (*ma)
					if (*lda) < (*nmax) {
						(*lda) = (*lda) + 1
					}
					//*                 Skip tests if not enough room.
					if (*lda) > (*nmax) {
						goto Label80
					}
					(*laa) = (*lda) * (*na)
					//*
					//*                 Generate the matrix A.
					//*
					ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, a, nmax, aa, lda, reset, zero)
					//*
					for (*icb) = 1; (*icb) <= 3; (*icb)++ {
						(*transb) = (*ich)[(*icb)-1]
						(*tranb) = (*transb) == "t" || (*transb) == "c"
						//*
						if *tranb {
							(*mb) = (*n)
							(*nb) = (*k)
						} else {
							(*mb) = (*k)
							(*nb) = (*n)
						}
						//*                    Set LDB to 1 more than minimum value if room.
						(*ldb) = (*mb)
						if (*ldb) < (*nmax) {
							(*ldb) = (*ldb) + 1
						}
						//*                    Skip tests if not enough room.
						if (*ldb) > (*nmax) {
							goto Label70
						}
						(*lbb) = (*ldb) * (*nb)
						//*
						//*                    Generate the matrix B.
						//*
						ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), mb, nb, b, nmax, bb, ldb, reset, zero)
						//*
						for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
							(*alpha) = (*alf)[(*ia)-1]
							//*
							for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
								(*beta) = (*bet)[(*ib)-1]
								//*
								//*                          Generate the matrix C.
								//*
								ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, c, nmax, cc, ldc, reset, zero)
								//*
								(*nc) = (*nc) + 1
								//*
								//*                          Save every datum before calling the
								//*                          subroutine.
								//*
								(*tranas) = (*transa)
								(*tranbs) = (*transb)
								(*ms) = (*m)
								(*ns) = (*n)
								(*ks) = (*k)
								(*als) = (*alpha)
								for (*i) = 1; (*i) <= (*laa); (*i)++ {
									(*as)[(*i)-1] = (*aa)[(*i)-1]
								//Label10:
								}
								(*ldas) = (*lda)
								for (*i) = 1; (*i) <= (*lbb); (*i)++ {
									(*bs)[(*i)-1] = (*bb)[(*i)-1]
								//Label20:
								}
								(*ldbs) = (*ldb)
								(*bls) = (*beta)
								for (*i) = 1; (*i) <= (*lcc); (*i)++ {
									(*cs)[(*i)-1] = (*cc)[(*i)-1]
								//Label30:
								}
								(*ldcs) = (*ldc)
								//*
								//*                          Call the subroutine.
								//*
								if *trace {
									WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c','%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d).\n"); return &y}(), (*nc), (*sname), (*transa), (*transb), (*m), (*n), (*k), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
								}
								if *rewi {
									REWIND(NTRA)
								}
								Zgemm(transa, transb, m, n, k, alpha, aa, lda, bb, ldb, beta, cc, ldc)
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
								(*isame)[0] = (*transa) == (*tranas)
								(*isame)[1] = (*transb) == (*tranbs)
								(*isame)[2] = (*ms) == (*m)
								(*isame)[3] = (*ns) == (*n)
								(*isame)[4] = (*ks) == (*k)
								(*isame)[5] = (*als) == (*alpha)
								(*isame)[6] = (*LZE(as, aa, laa))
								(*isame)[7] = (*ldas) == (*lda)
								(*isame)[8] = (*LZE(bs, bb, lbb))
								(*isame)[9] = (*ldbs) == (*ldb)
								(*isame)[10] = (*bls) == (*beta)
								if *null {
									(*isame)[11] = (*LZE(cs, cc, lcc))
								} else {
									(*isame)[11] = (*LZERES(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, cs, cc, ldc))
								}
								(*isame)[12] = (*ldcs) == (*ldc)
								//*
								//*                          If data was incorrectly changed, report
								//*                          and return.
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
									zmmch ((*transa), (*transb), (*m), (*n), (*k), (*alpha), (*a), (*nmax), (*b), (*nmax), (*beta), (*c), (*nmax), (*ct), (*g), (*cc), (*ldc), (*eps), (*err), (*fatal), (*nout), true)
									(*errmax) = (MAX((*errmax), (*err)))
									//*                             If got really bad answer, report and
									//*                             return.
									if *fatal {
										goto Label120
									}
								}
								//*
							//Label50:
							}
							//*
						//Label60:
						}
						//*
					Label70:
					}
					//*
				Label80:
				}
				//*
			//Label90:
			}
			//*
		Label100:
		}
		//*
	//Label110:
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
	WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c','%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d).\n"); return &y}(), (*nc), (*sname), (*transa), (*transb), (*m), (*n), (*k), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
	//*
Label130:
	;
	return
	//*
	//*
	//*     End of ZCHK1.
	//*
}

func Zchk2(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, nmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, b *[][]complex128, bb *[]complex128, bs *[]complex128, c *[][]complex128, cc *[]complex128, cs *[]complex128, ct *[]complex128, g *[]float64) {
	zero := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bls := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ics := new(int)
	icu := new(int)
	im := new(int)
	in := new(int)
	laa := new(int)
	lbb := new(int)
	lcc := new(int)
	lda := new(int)
	ldas := new(int)
	ldb := new(int)
	ldbs := new(int)
	ldc := new(int)
	ldcs := new(int)
	m := new(int)
	ms := new(int)
	n := new(int)
	na := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	conj := new(bool)
	left := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	side := new(byte)
	sides := new(byte)
	uplo := new(byte)
	uplos := new(byte)
	ichs := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	ichu := func() *[]byte {
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
	//*  Tests Zhemm and Zsymm.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
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
	(*ichs)[0], (*ichs)[1], (*ichu)[0], (*ichu)[1] = 'l', 'r', 'u', 'l'
	//*     .. Executable Statements ..
	(*conj) = (*sname)[2     - 1] == "he"
	//*
	(*nargs) = 12
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*im) = 1; (*im) <= (*nidim); (*im)++ {
		(*m) = (*idim)[(*im)-1]
		//*
		for (*in) = 1; (*in) <= (*nidim); (*in)++ {
			(*n) = (*idim)[(*in)-1]
			//*           Set LDC to 1 more than minimum value if room.
			(*ldc) = (*m)
			if (*ldc) < (*nmax) {
				(*ldc) = (*ldc) + 1
			}
			//*           Skip tests if not enough room.
			if (*ldc) > (*nmax) {
				goto Label90
			}
			(*lcc) = (*ldc) * (*n)
			(*null) = (*n) <= 0|| (*m) <= 0
			//*           Set LDB to 1 more than minimum value if room.
			(*ldb) = (*m)
			if (*ldb) < (*nmax) {
				(*ldb) = (*ldb) + 1
			}
			//*           Skip tests if not enough room.
			if (*ldb) > (*nmax) {
				goto Label90
			}
			(*lbb) = (*ldb) * (*n)
			//*
			//*           Generate the matrix B.
			//*
			ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, b, nmax, bb, ldb, reset, zero)
			//*
			for (*ics) = 1; (*ics) <= 2; (*ics)++ {
				(*side) = (*ichs)[(*ics)-1]
				(*left) = (*side) == "l"
				//*
				if *left {
					(*na) = (*m)
				} else {
					(*na) = (*n)
				}
				//*              Set LDA to 1 more than minimum value if room.
				(*lda) = (*na)
				if (*lda) < (*nmax) {
					(*lda) = (*lda) + 1
				}
				//*              Skip tests if not enough room.
				if (*lda) > (*nmax) {
					goto Label80
				}
				(*laa) = (*lda) * (*na)
				//*
				for (*icu) = 1; (*icu) <= 2; (*icu)++ {
					(*uplo) = (*ichu)[(*icu)-1]
					//*
					//*                 Generate the hermitian or symmetric matrix A.
					//*
					ZMAKE(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), na, na, a, nmax, aa, lda, reset, zero)
					//*
					for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
						(*alpha) = (*alf)[(*ia)-1]
						//*
						for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
							(*beta) = (*bet)[(*ib)-1]
							//*
							//*                       Generate the matrix C.
							//*
							ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, c, nmax, cc, ldc, reset, zero)
							//*
							(*nc) = (*nc) + 1
							//*
							//*                       Save every datum before calling the
							//*                       subroutine.
							//*
							(*sides) = (*side)
							(*uplos) = (*uplo)
							(*ms) = (*m)
							(*ns) = (*n)
							(*als) = (*alpha)
							for (*i) = 1; (*i) <= (*laa); (*i)++ {
								(*as)[(*i)-1] = (*aa)[(*i)-1]
							//Label10:
							}
							(*ldas) = (*lda)
							for (*i) = 1; (*i) <= (*lbb); (*i)++ {
								(*bs)[(*i)-1] = (*bb)[(*i)-1]
							//Label20:
							}
							(*ldbs) = (*ldb)
							(*bls) = (*beta)
							for (*i) = 1; (*i) <= (*lcc); (*i)++ {
								(*cs)[(*i)-1] = (*cc)[(*i)-1]
							//Label30:
							}
							(*ldcs) = (*ldc)
							//*
							//*                       Call the subroutine.
							//*
							if *trace {
								WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d)    .\n"); return &y}(), (*nc), (*sname), (*side), (*uplo), (*m), (*n), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
							}
							if *rewi {
								REWIND(NTRA)
							}
							if *conj {
								Zhemm(side, uplo, m, n, alpha, aa, lda, bb, ldb, beta, cc, ldc)
							} else {
								Zsymm(side, uplo, m, n, alpha, aa, lda, bb, ldb, beta, cc, ldc)
							}
							//*
							//*                       Check if error-exit was taken incorrectly.
							//*
							if !(*ok) {
								WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
								(*fatal) = true
								goto Label110
							}
							//*
							//*                       See what data changed inside subroutines.
							//*
							(*isame)[0] = (*sides) == (*side)
							(*isame)[1] = (*uplos) == (*uplo)
							(*isame)[2] = (*ms) == (*m)
							(*isame)[3] = (*ns) == (*n)
							(*isame)[4] = (*als) == (*alpha)
							(*isame)[5] = (*LZE(as, aa, laa))
							(*isame)[6] = (*ldas) == (*lda)
							(*isame)[7] = (*LZE(bs, bb, lbb))
							(*isame)[8] = (*ldbs) == (*ldb)
							(*isame)[9] = (*bls) == (*beta)
							if *null {
								(*isame)[10] = (*LZE(cs, cc, lcc))
							} else {
								(*isame)[10] = (*LZERES(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, cs, cc, ldc))
							}
							(*isame)[11] = (*ldcs) == (*ldc)
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
								goto Label110
							}
							//*
							if !(*null) {
								//*
								//*                          Check the result.
								//*
								if *left {
									zmmch ( "n", "n", (*m), (*n), (*m), (*alpha), (*a), (*nmax), (*b), (*nmax), (*beta), (*c), (*nmax), (*ct), (*g), (*cc), (*ldc), (*eps), (*err), (*fatal), (*nout), true)
								} else {
									zmmch ( "n", "n", (*m), (*n), (*n), (*alpha), (*b), (*nmax), (*a), (*nmax), (*beta), (*c), (*nmax), (*ct), (*g), (*cc), (*ldc), (*eps), (*err), (*fatal), (*nout), true)
								}
								(*errmax) = (MAX((*errmax), (*err)))
								//*                          If got really bad answer, report and
								//*                          return.
								if *fatal {
									goto Label110
								}
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
			Label80:
			}
			//*
		Label90:
		}
		//*
	//Label100:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label120
	//*
Label110:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d)    .\n"); return &y}(), (*nc), (*sname), (*side), (*uplo), (*m), (*n), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
	//*
Label120:
	;
	return
	//*
	//*
	//*     End of ZCHK2.
	//*
}

func Zchk3(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, nmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, b *[][]complex128, bb *[]complex128, bs *[]complex128, ct *[]complex128, g *[]float64, c *[][]complex128) {
	zero := new(complex128)
	one := new(complex128)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	err := new(float64)
	errmax := new(float64)
	i := new(int)
	ia := new(int)
	icd := new(int)
	ics := new(int)
	ict := new(int)
	icu := new(int)
	im := new(int)
	in := new(int)
	j := new(int)
	laa := new(int)
	lbb := new(int)
	lda := new(int)
	ldas := new(int)
	ldb := new(int)
	ldbs := new(int)
	m := new(int)
	ms := new(int)
	n := new(int)
	na := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	left := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	diag := new(byte)
	diags := new(byte)
	side := new(byte)
	sides := new(byte)
	tranas := new(byte)
	transa := new(byte)
	uplo := new(byte)
	uplos := new(byte)
	ichd := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	ichs := func() *[]byte {
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
	//*  Tests Ztrmm and Ztrsm.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
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
	(*ichu)[0], (*ichu)[1], (*icht)[0], (*icht)[1], (*icht)[2], (*ichd)[0], (*ichd)[1], (*ichs)[0], (*ichs)[1] = 'u', 'l', 'n', 't', 'c', 'u', 'n', 'l', 'r'
	//*     .. Executable Statements ..
	//*
	(*nargs) = 11
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*     Set up zero matrix for ZMMCH.
	for (*j) = 1; (*j) <= (*nmax); (*j)++ {
		for (*i) = 1; (*i) <= (*nmax); (*i)++ {
			(*c)[(*i)-1][(*j)-1] = (*zero)
		//Label10:
		}
	//Label20:
	}
	//*
	for (*im) = 1; (*im) <= (*nidim); (*im)++ {
		(*m) = (*idim)[(*im)-1]
		//*
		for (*in) = 1; (*in) <= (*nidim); (*in)++ {
			(*n) = (*idim)[(*in)-1]
			//*           Set LDB to 1 more than minimum value if room.
			(*ldb) = (*m)
			if (*ldb) < (*nmax) {
				(*ldb) = (*ldb) + 1
			}
			//*           Skip tests if not enough room.
			if (*ldb) > (*nmax) {
				goto Label130
			}
			(*lbb) = (*ldb) * (*n)
			(*null) = (*m) <= 0|| (*n) <= 0
			//*
			for (*ics) = 1; (*ics) <= 2; (*ics)++ {
				(*side) = (*ichs)[(*ics)-1]
				(*left) = (*side) == "l"
				if *left {
					(*na) = (*m)
				} else {
					(*na) = (*n)
				}
				//*              Set LDA to 1 more than minimum value if room.
				(*lda) = (*na)
				if (*lda) < (*nmax) {
					(*lda) = (*lda) + 1
				}
				//*              Skip tests if not enough room.
				if (*lda) > (*nmax) {
					goto Label130
				}
				(*laa) = (*lda) * (*na)
				//*
				for (*icu) = 1; (*icu) <= 2; (*icu)++ {
					(*uplo) = (*ichu)[(*icu)-1]
					//*
					for (*ict) = 1; (*ict) <= 3; (*ict)++ {
						(*transa) = (*icht)[(*ict)-1]
						//*
						for (*icd) = 1; (*icd) <= 2; (*icd)++ {
							(*diag) = (*ichd)[(*icd)-1]
							//*
							for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
								(*alpha) = (*alf)[(*ia)-1]
								//*
								//*                          Generate the matrix A.
								//*
								ZMAKE(func() *[]byte{y := []byte("tr"); return &y}(), uplo, diag, na, na, a, nmax, aa, lda, reset, zero)
								//*
								//*                          Generate the matrix B.
								//*
								ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, b, nmax, bb, ldb, reset, zero)
								//*
								(*nc) = (*nc) + 1
								//*
								//*                          Save every datum before calling the
								//*                          subroutine.
								//*
								(*sides) = (*side)
								(*uplos) = (*uplo)
								(*tranas) = (*transa)
								(*diags) = (*diag)
								(*ms) = (*m)
								(*ns) = (*n)
								(*als) = (*alpha)
								for (*i) = 1; (*i) <= (*laa); (*i)++ {
									(*as)[(*i)-1] = (*aa)[(*i)-1]
								//Label30:
								}
								(*ldas) = (*lda)
								for (*i) = 1; (*i) <= (*lbb); (*i)++ {
									(*bs)[(*i)-1] = (*bb)[(*i)-1]
								//Label40:
								}
								(*ldbs) = (*ldb)
								//*
								//*                          Call the subroutine.
								//*
								if (*sname)[4     - 1] == "mm" {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d)               .\n"); return &y}(), (*nc), (*sname), (*side), (*uplo), (*transa), (*diag), (*m), (*n), (*alpha), (*lda), (*ldb))
									}
									if *rewi {
										REWIND(NTRA)
									}
									Ztrmm(side, uplo, transa, diag, m, n, alpha, aa, lda, bb, ldb)
								} else if (*sname)[4     - 1] == "sm" {
									if *trace {
										WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d)               .\n"); return &y}(), (*nc), (*sname), (*side), (*uplo), (*transa), (*diag), (*m), (*n), (*alpha), (*lda), (*ldb))
									}
									if *rewi {
										REWIND(NTRA)
									}
									Ztrsm(side, uplo, transa, diag, m, n, alpha, aa, lda, bb, ldb)
								}
								//*
								//*                          Check if error-exit was taken incorrectly.
								//*
								if !(*ok) {
									WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
									(*fatal) = true
									goto Label150
								}
								//*
								//*                          See what data changed inside subroutines.
								//*
								(*isame)[0] = (*sides) == (*side)
								(*isame)[1] = (*uplos) == (*uplo)
								(*isame)[2] = (*tranas) == (*transa)
								(*isame)[3] = (*diags) == (*diag)
								(*isame)[4] = (*ms) == (*m)
								(*isame)[5] = (*ns) == (*n)
								(*isame)[6] = (*als) == (*alpha)
								(*isame)[7] = (*LZE(as, aa, laa))
								(*isame)[8] = (*ldas) == (*lda)
								if *null {
									(*isame)[9] = (*LZE(bs, bb, lbb))
								} else {
									(*isame)[9] = (*LZERES(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), m, n, bs, bb, ldb))
								}
								(*isame)[10] = (*ldbs) == (*ldb)
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
								//Label50:
								}
								if !(*same) {
									(*fatal) = true
									goto Label150
								}
								//*
								if !(*null) {
									if (*sname)[4     - 1] == "mm" {
										//*
										//*                                Check the result.
										//*
										if *left {
											zmmch ((*transa), "n", (*m), (*n), (*m), (*alpha), (*a), (*nmax), (*b), (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*bb), (*ldb), (*eps), (*err), (*fatal), (*nout), true)
										} else {
											zmmch ( "n", (*transa), (*m), (*n), (*n), (*alpha), (*b), (*nmax), (*a), (*nmax), (*zero), (*c), (*nmax), (*ct), (*g), (*bb), (*ldb), (*eps), (*err), (*fatal), (*nout), true)
										}
									} else if (*sname)[4     - 1] == "sm" {
										//*
										//*                                Compute approximation to original
										//*                                matrix.
										//*
										for (*j) = 1; (*j) <= (*n); (*j)++ {
											for (*i) = 1; (*i) <= (*m); (*i)++ {
												(*c)[(*i)-1][(*j)-1] = (*bb)[(*i)+((*j)-1)*(*ldb)-1]
												(*bb)[(*i)+((*j)-1)*(*ldb)-1] = (*alpha) * (*b)[(*i)-1][(*j)-1]
											//Label60:
											}
										//Label70:
										}
										//*
										if *left {
											zmmch ((*transa), "n", (*m), (*n), (*m), (*one), (*a), (*nmax), (*c), (*nmax), (*zero), (*b), (*nmax), (*ct), (*g), (*bb), (*ldb), (*eps), (*err), (*fatal), (*nout), false)
										} else {
											zmmch ( "n", (*transa), (*m), (*n), (*n), (*one), (*c), (*nmax), (*a), (*nmax), (*zero), (*b), (*nmax), (*ct), (*g), (*bb), (*ldb), (*eps), (*err), (*fatal), (*nout), false)
										}
									}
									(*errmax) = (MAX((*errmax), (*err)))
									//*                             If got really bad answer, report and
									//*                             return.
									if *fatal {
										goto Label150
									}
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
				//Label110:
				}
				//*
			//Label120:
			}
			//*
		Label130:
		}
		//*
	//Label140:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label160
	//*
Label150:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d)               .\n"); return &y}(), (*nc), (*sname), (*side), (*uplo), (*transa), (*diag), (*m), (*n), (*alpha), (*lda), (*ldb))
	//*
Label160:
	;
	return
	//*
	//*
	//*     End of ZCHK3.
	//*
}

func Zchk4(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, nmax *int, a *[][]complex128, aa *[]complex128, as *[]complex128, b *[][]complex128, bb *[]complex128, bs *[]complex128, c *[][]complex128, cc *[]complex128, cs *[]complex128, ct *[]complex128, g *[]float64) {
	zero := new(complex128)
	rone := new(float64)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bets := new(complex128)
	err := new(float64)
	errmax := new(float64)
	ralpha := new(float64)
	rals := new(float64)
	rbeta := new(float64)
	rbets := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ict := new(int)
	icu := new(int)
	ik := new(int)
	in := new(int)
	j := new(int)
	jc := new(int)
	jj := new(int)
	k := new(int)
	ks := new(int)
	laa := new(int)
	lcc := new(int)
	lda := new(int)
	ldas := new(int)
	ldc := new(int)
	ldcs := new(int)
	lj := new(int)
	ma := new(int)
	n := new(int)
	na := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	conj := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	tran := new(bool)
	upper := new(bool)
	trans := new(byte)
	transs := new(byte)
	transt := new(byte)
	uplo := new(byte)
	uplos := new(byte)
	icht := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	ichu := func() *[]byte {
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
	//*  Tests Zherk and Zsyrk.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*rone) = 1.0
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
	(*icht)[0], (*icht)[1], (*ichu)[0], (*ichu)[1] = 'n', 'c', 'u', 'l'
	//*     .. Executable Statements ..
	(*conj) = (*sname)[2     - 1] == "he"
	//*
	(*nargs) = 10
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*        Set LDC to 1 more than minimum value if room.
		(*ldc) = (*n)
		if (*ldc) < (*nmax) {
			(*ldc) = (*ldc) + 1
		}
		//*        Skip tests if not enough room.
		if (*ldc) > (*nmax) {
			goto Label100
		}
		(*lcc) = (*ldc) * (*n)
		//*
		for (*ik) = 1; (*ik) <= (*nidim); (*ik)++ {
			(*k) = (*idim)[(*ik)-1]
			//*
			for (*ict) = 1; (*ict) <= 2; (*ict)++ {
				(*trans) = (*icht)[(*ict)-1]
				(*tran) = (*trans) == "c"
				if (*tran) && !(*conj) {
					(*trans) = 't'
				}
				if *tran {
					(*ma) = (*k)
					(*na) = (*n)
				} else {
					(*ma) = (*n)
					(*na) = (*k)
				}
				//*              Set LDA to 1 more than minimum value if room.
				(*lda) = (*ma)
				if (*lda) < (*nmax) {
					(*lda) = (*lda) + 1
				}
				//*              Skip tests if not enough room.
				if (*lda) > (*nmax) {
					goto Label80
				}
				(*laa) = (*lda) * (*na)
				//*
				//*              Generate the matrix A.
				//*
				ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, a, nmax, aa, lda, reset, zero)
				//*
				for (*icu) = 1; (*icu) <= 2; (*icu)++ {
					(*uplo) = (*ichu)[(*icu)-1]
					(*upper) = (*uplo) == "u"
					//*
					for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
						(*alpha) = (*alf)[(*ia)-1]
						if *conj {
							(*ralpha) = (DBLE((*alpha)))
							(*alpha) = (*DCMPLX(ralpha, rzero))
						}
						//*
						for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
							(*beta) = (*bet)[(*ib)-1]
							if *conj {
								(*rbeta) = (DBLE((*beta)))
								(*beta) = (*DCMPLX(rbeta, rzero))
							}
							(*null) = (*n) <= 0
							if *conj {
								(*null) = (*null) || ( ((*k) <= 0|| (*ralpha) == (*rzero)) && (*rbeta) == (*rone))
							}
							//*
							//*                       Generate the matrix C.
							//*
							ZMAKE(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), n, n, c, nmax, cc, ldc, reset, zero)
							//*
							(*nc) = (*nc) + 1
							//*
							//*                       Save every datum before calling the subroutine.
							//*
							(*uplos) = (*uplo)
							(*transs) = (*trans)
							(*ns) = (*n)
							(*ks) = (*k)
							if *conj {
								(*rals) = (*ralpha)
							} else {
								(*als) = (*alpha)
							}
							for (*i) = 1; (*i) <= (*laa); (*i)++ {
								(*as)[(*i)-1] = (*aa)[(*i)-1]
							//Label10:
							}
							(*ldas) = (*lda)
							if *conj {
								(*rbets) = (*rbeta)
							} else {
								(*bets) = (*beta)
							}
							for (*i) = 1; (*i) <= (*lcc); (*i)++ {
								(*cs)[(*i)-1] = (*cc)[(*i)-1]
							//Label20:
							}
							(*ldcs) = (*ldc)
							//*
							//*                       Call the subroutine.
							//*
							if *conj {
								if *trace {
									WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, a,%3d,%4.1f, c,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*ralpha), (*lda), (*rbeta), (*ldc))
								}
								if *rewi {
									REWIND(NTRA)
								}
								Zherk(uplo, trans, n, k, ralpha, aa, lda, rbeta, cc, ldc)
							} else {
								if *trace {
									WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d,(%4.1f,%4.1f), c,%3d)          .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*beta), (*ldc))
								}
								if *rewi {
									REWIND(NTRA)
								}
								Zsyrk(uplo, trans, n, k, alpha, aa, lda, beta, cc, ldc)
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
							(*isame)[0] = (*uplos) == (*uplo)
							(*isame)[1] = (*transs) == (*trans)
							(*isame)[2] = (*ns) == (*n)
							(*isame)[3] = (*ks) == (*k)
							if *conj {
								(*isame)[4] = (*rals) == (*ralpha)
							} else {
								(*isame)[4] = (*als) == (*alpha)
							}
							(*isame)[5] = (*LZE(as, aa, laa))
							(*isame)[6] = (*ldas) == (*lda)
							if *conj {
								(*isame)[7] = (*rbets) == (*rbeta)
							} else {
								(*isame)[7] = (*bets) == (*beta)
							}
							if *null {
								(*isame)[8] = (*LZE(cs, cc, lcc))
							} else {
								(*isame)[8] = (*LZERES(&((*sname)[1]), uplo, n, n, cs, cc, ldc))
							}
							(*isame)[9] = (*ldcs) == (*ldc)
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
							//Label30:
							}
							if !(*same) {
								(*fatal) = true
								goto Label120
							}
							//*
							if !(*null) {
								//*
								//*                          Check the result column by column.
								//*
								if *conj {
									(*transt) = 'c'
								} else {
									(*transt) = 't'
								}
								(*jc) = 1
								for (*j) = 1; (*j) <= (*n); (*j)++ {
									if *upper {
										(*jj) = 1
										(*lj) = (*j)
									} else {
										(*jj) = (*j)
										(*lj) = (*n) - (*j) + 1
									}
									if *tran {
										zmmch ((*transt), "n", (*lj), 1, (*k), (*alpha), (*a)[0][(*jj) - 1], (*nmax), (*a)[0][(*j) - 1], (*nmax), (*beta), (*c)[(*jj) - 1][(*j) - 1], (*nmax), (*ct), (*g), (*cc)[(*jc) - 1], (*ldc), (*eps), (*err), (*fatal), (*nout), true)
									} else {
										zmmch ( "n", (*transt), (*lj), 1, (*k), (*alpha), (*a)[(*jj) - 1][0], (*nmax), (*a)[(*j) - 1][0], (*nmax), (*beta), (*c)[(*jj) - 1][(*j) - 1], (*nmax), (*ct), (*g), (*cc)[(*jc) - 1], (*ldc), (*eps), (*err), (*fatal), (*nout), true)
									}
									if *upper {
										(*jc) = (*jc) + (*ldc)
									} else {
										(*jc) = (*jc) + (*ldc) + 1
									}
									(*errmax) = (MAX((*errmax), (*err)))
									//*                             If got really bad answer, report and
									//*                             return.
									if *fatal {
										goto Label110
									}
								//Label40:
								}
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
			Label80:
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
	if (*n) > 1 {
		WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	}
	//*
Label120:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *conj {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,%4.1f, a,%3d,%4.1f, c,%3d)                         .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*ralpha), (*lda), (*rbeta), (*ldc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d,(%4.1f,%4.1f), c,%3d)          .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*beta), (*ldc))
	}
	//*
Label130:
	;
	return
	//*
	//*
	//*     End of ZCHK4.
	//*
}

func Zchk5(sname *[]byte, eps *float64, thresh *float64, nout *int, ntra *int, trace *bool, rewi *bool, fatal *bool, nidim *int, idim *[]int, nalf *int, alf *[]complex128, nbet *int, bet *[]complex128, nmax *int, ab *[]complex128, aa *[]complex128, as *[]complex128, bb *[]complex128, bs *[]complex128, c *[][]complex128, cc *[]complex128, cs *[]complex128, ct *[]complex128, g *[]float64, w *[]complex128) {
	zero := new(complex128)
	one := new(complex128)
	rone := new(float64)
	rzero := new(float64)
	alpha := new(complex128)
	als := new(complex128)
	beta := new(complex128)
	bets := new(complex128)
	err := new(float64)
	errmax := new(float64)
	rbeta := new(float64)
	rbets := new(float64)
	i := new(int)
	ia := new(int)
	ib := new(int)
	ict := new(int)
	icu := new(int)
	ik := new(int)
	in := new(int)
	j := new(int)
	jc := new(int)
	jj := new(int)
	jjab := new(int)
	k := new(int)
	ks := new(int)
	laa := new(int)
	lbb := new(int)
	lcc := new(int)
	lda := new(int)
	ldas := new(int)
	ldb := new(int)
	ldbs := new(int)
	ldc := new(int)
	ldcs := new(int)
	lj := new(int)
	ma := new(int)
	n := new(int)
	na := new(int)
	nargs := new(int)
	nc := new(int)
	ns := new(int)
	conj := new(bool)
	null := new(bool)
	reset := new(bool)
	same := new(bool)
	tran := new(bool)
	upper := new(bool)
	trans := new(byte)
	transs := new(byte)
	transt := new(byte)
	uplo := new(byte)
	uplos := new(byte)
	icht := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	ichu := func() *[]byte {
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
	//*  Tests Zher2k and Zsyr2k.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Parameters ..
	(*zero) = (0.0 + (0.0)*1i)
	(*one) = (1.0 + (0.0)*1i)
	(*rone) = 1.0
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
	(*icht)[0], (*icht)[1], (*ichu)[0], (*ichu)[1] = 'n', 'c', 'u', 'l'
	//*     .. Executable Statements ..
	(*conj) = (*sname)[2     - 1] == "he"
	//*
	(*nargs) = 12
	(*nc) = 0
	(*reset) = true
	(*errmax) = (*rzero)
	//*
	for (*in) = 1; (*in) <= (*nidim); (*in)++ {
		(*n) = (*idim)[(*in)-1]
		//*        Set LDC to 1 more than minimum value if room.
		(*ldc) = (*n)
		if (*ldc) < (*nmax) {
			(*ldc) = (*ldc) + 1
		}
		//*        Skip tests if not enough room.
		if (*ldc) > (*nmax) {
			goto Label130
		}
		(*lcc) = (*ldc) * (*n)
		//*
		for (*ik) = 1; (*ik) <= (*nidim); (*ik)++ {
			(*k) = (*idim)[(*ik)-1]
			//*
			for (*ict) = 1; (*ict) <= 2; (*ict)++ {
				(*trans) = (*icht)[(*ict)-1]
				(*tran) = (*trans) == "c"
				if (*tran) && !(*conj) {
					(*trans) = 't'
				}
				if *tran {
					(*ma) = (*k)
					(*na) = (*n)
				} else {
					(*ma) = (*n)
					(*na) = (*k)
				}
				//*              Set LDA to 1 more than minimum value if room.
				(*lda) = (*ma)
				if (*lda) < (*nmax) {
					(*lda) = (*lda) + 1
				}
				//*              Skip tests if not enough room.
				if (*lda) > (*nmax) {
					goto Label110
				}
				(*laa) = (*lda) * (*na)
				//*
				//*              Generate the matrix A.
				//*
				if *tran {
					ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, ab, 2*(*nmax), aa, lda, reset, zero)
				} else {
					ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, ab, nmax, aa, lda, reset, zero)
				}
				//*
				//*              Generate the matrix B.
				//*
				(*ldb) = (*lda)
				(*lbb) = (*laa)
				if *tran {
					ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, &((*ab)[(*k)+0]), 2*(*nmax), bb, ldb, reset, zero)
				} else {
					ZMAKE(func() *[]byte{y := []byte("ge"); return &y}(), func() *byte{y := byte(' '); return &y}(), func() *byte{y := byte(' '); return &y}(), ma, na, &((*ab)[(*k)*(*nmax)+0]), nmax, bb, ldb, reset, zero)
				}
				//*
				for (*icu) = 1; (*icu) <= 2; (*icu)++ {
					(*uplo) = (*ichu)[(*icu)-1]
					(*upper) = (*uplo) == "u"
					//*
					for (*ia) = 1; (*ia) <= (*nalf); (*ia)++ {
						(*alpha) = (*alf)[(*ia)-1]
						//*
						for (*ib) = 1; (*ib) <= (*nbet); (*ib)++ {
							(*beta) = (*bet)[(*ib)-1]
							if *conj {
								(*rbeta) = (DBLE((*beta)))
								(*beta) = (*DCMPLX(rbeta, rzero))
							}
							(*null) = (*n) <= 0
							if *conj {
								(*null) = (*null) || ( ((*k) <= 0|| (*alpha) == (*zero)) && (*rbeta) == (*rone))
							}
							//*
							//*                       Generate the matrix C.
							//*
							ZMAKE(&((*sname)[1]), uplo, func() *byte{y := byte(' '); return &y}(), n, n, c, nmax, cc, ldc, reset, zero)
							//*
							(*nc) = (*nc) + 1
							//*
							//*                       Save every datum before calling the subroutine.
							//*
							(*uplos) = (*uplo)
							(*transs) = (*trans)
							(*ns) = (*n)
							(*ks) = (*k)
							(*als) = (*alpha)
							for (*i) = 1; (*i) <= (*laa); (*i)++ {
								(*as)[(*i)-1] = (*aa)[(*i)-1]
							//Label10:
							}
							(*ldas) = (*lda)
							for (*i) = 1; (*i) <= (*lbb); (*i)++ {
								(*bs)[(*i)-1] = (*bb)[(*i)-1]
							//Label20:
							}
							(*ldbs) = (*ldb)
							if *conj {
								(*rbets) = (*rbeta)
							} else {
								(*bets) = (*beta)
							}
							for (*i) = 1; (*i) <= (*lcc); (*i)++ {
								(*cs)[(*i)-1] = (*cc)[(*i)-1]
							//Label30:
							}
							(*ldcs) = (*ldc)
							//*
							//*                       Call the subroutine.
							//*
							if *conj {
								if *trace {
									WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,%4.1f, c,%3d)           .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*ldb), (*rbeta), (*ldc))
								}
								if *rewi {
									REWIND(NTRA)
								}
								Zher2k(uplo, trans, n, k, alpha, aa, lda, bb, ldb, rbeta, cc, ldc)
							} else {
								if *trace {
									WRITE((*ntra), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d)    .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
								}
								if *rewi {
									REWIND(NTRA)
								}
								Zsyr2k(uplo, trans, n, k, alpha, aa, lda, bb, ldb, beta, cc, ldc)
							}
							//*
							//*                       Check if error-exit was taken incorrectly.
							//*
							if !(*ok) {
								WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - error-exit taken on valid call *******\n"); return &y}())
								(*fatal) = true
								goto Label150
							}
							//*
							//*                       See what data changed inside subroutines.
							//*
							(*isame)[0] = (*uplos) == (*uplo)
							(*isame)[1] = (*transs) == (*trans)
							(*isame)[2] = (*ns) == (*n)
							(*isame)[3] = (*ks) == (*k)
							(*isame)[4] = (*als) == (*alpha)
							(*isame)[5] = (*LZE(as, aa, laa))
							(*isame)[6] = (*ldas) == (*lda)
							(*isame)[7] = (*LZE(bs, bb, lbb))
							(*isame)[8] = (*ldbs) == (*ldb)
							if *conj {
								(*isame)[9] = (*rbets) == (*rbeta)
							} else {
								(*isame)[9] = (*bets) == (*beta)
							}
							if *null {
								(*isame)[10] = (*LZE(cs, cc, lcc))
							} else {
								(*isame)[10] = (*LZERES(func() *[]byte{y := []byte("he"); return &y}(), uplo, n, n, cs, cc, ldc))
							}
							(*isame)[11] = (*ldcs) == (*ldc)
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
								goto Label150
							}
							//*
							if !(*null) {
								//*
								//*                          Check the result column by column.
								//*
								if *conj {
									(*transt) = 'c'
								} else {
									(*transt) = 't'
								}
								(*jjab) = 1
								(*jc) = 1
								for (*j) = 1; (*j) <= (*n); (*j)++ {
									if *upper {
										(*jj) = 1
										(*lj) = (*j)
									} else {
										(*jj) = (*j)
										(*lj) = (*n) - (*j) + 1
									}
									if *tran {
										for (*i) = 1; (*i) <= (*k); (*i)++ {
											(*w)[(*i)-1] = (*alpha) * (*ab)[((*j)-1)*2*(*nmax)+(*k)+(*i)-1]
											if *conj {
												(*w)[(*k)+(*i)-1] = DCONJG((*alpha)) * (*ab)[((*j)-1)*2*(*nmax)+(*i)-1]
											} else {
												(*w)[(*k)+(*i)-1] = (*alpha) * (*ab)[((*j)-1)*2*(*nmax)+(*i)-1]
											}
										//Label50:
										}
										zmmch ((*transt), "n", (*lj), 1, 2 * (*k), (*one), (*ab)[(*jjab) - 1], 2 * (*nmax), (*w), 2 * (*nmax), (*beta), (*c)[(*jj) - 1][(*j) - 1], (*nmax), (*ct), (*g), (*cc)[(*jc) - 1], (*ldc), (*eps), (*err), (*fatal), (*nout), true)
									} else {
										for (*i) = 1; (*i) <= (*k); (*i)++ {
											if *conj {
												(*w)[(*i)-1] = (*alpha) * DCONJG(((*ab)[((*k)+(*i)-1)*(*nmax)+(*j)-1]))
												(*w)[(*k)+(*i)-1] = (DCONJG((*alpha) * (*ab)[((*i)-1)*(*nmax)+(*j)-1]))
											} else {
												(*w)[(*i)-1] = (*alpha) * (*ab)[((*k)+(*i)-1)*(*nmax)+(*j)-1]
												(*w)[(*k)+(*i)-1] = (*alpha) * (*ab)[((*i)-1)*(*nmax)+(*j)-1]
											}
										//Label60:
										}
										zmmch ( "n", "n", (*lj), 1, 2 * (*k), (*one), (*ab)[(*jj) - 1], (*nmax), (*w), 2 * (*nmax), (*beta), (*c)[(*jj) - 1][(*j) - 1], (*nmax), (*ct), (*g), (*cc)[(*jc) - 1], (*ldc), (*eps), (*err), (*fatal), (*nout), true)
									}
									if *upper {
										(*jc) = (*jc) + (*ldc)
									} else {
										(*jc) = (*jc) + (*ldc) + 1
										if *tran {
											(*jjab) = (*jjab) + 2*(*nmax)
										}
									}
									(*errmax) = (MAX((*errmax), (*err)))
									//*                             If got really bad answer, report and
									//*                             return.
									if *fatal {
										goto Label140
									}
								//Label70:
								}
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
	Label130:
	}
	//*
	//*     Report result.
	//*
	if (*errmax) < (*thresh) {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the computational tests (%6d calls)\n"); return &y}(), (*sname), (*nc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s completed the computational tests (%6d calls)\n ******* but with maximum test ratio%8.2f - suspect *******\n"); return &y}(), (*sname), (*nc), (*errmax))
	}
	goto Label160
	//*
Label140:
	;
	if (*n) > 1 {
		WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	}
	//*
Label150:
	;
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed on call number:\n"); return &y}(), (*sname))
	if *conj {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,%4.1f, c,%3d)           .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*ldb), (*rbeta), (*ldc))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6d: %6s('%c',%3d,(%4.1f,%4.1f), a,%3d, b,%3d,(%4.1f,%4.1f), c,%3d)    .\n"); return &y}(), (*nc), (*sname), (*uplo), (*trans), (*n), (*k), (*alpha), (*lda), (*ldb), (*beta), (*ldc))
	}
	//*
Label160:
	;
	return
	//*
	//*
	//*     End of ZCHK5.
	//*
}

func Zchke(isnum *int, srnamt *[]byte, nout *int) {
	infot := new(int)
	noutc := new(int)
	lerr := new(bool)
	ok := new(bool)
	one := new(float64)
	two := new(float64)
	alpha := new(complex128)
	beta := new(complex128)
	ralpha := new(float64)
	rbeta := new(float64)
	a := func() *[][]complex128 {
		arr := make([][]complex128, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]complex128, 1)
		}
		return &arr
	}()
	b := func() *[][]complex128 {
		arr := make([][]complex128, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]complex128, 1)
		}
		return &arr
	}()
	c := func() *[][]complex128 {
		arr := make([][]complex128, 2)
		for u := 0; u < 2; u++ {
			arr[u] = make([]complex128, 1)
		}
		return &arr
	}()
	common.infoc.lerr := new(float64)
	common.infoc.ok := new(float64)
	common.infoc.noutc := new(float64)
	common.infoc.infot := new(int)
	//*
	//*  Tests the error exits from the Level 3 Blas.
	//*  Requires a special version of the error-handling routine Xerbla.
	//*  A, B and C should not need to be defined.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*  3-19-92:  Initialize ALPHA, BETA, RALPHA, and RBETA  (eca)
	//*  3-19-92:  Fix argument 12 in calls to Zsymm and Zhemm
	//*            with infot = 9  (eca)
	//*  10-9-00:  Declared INTRINSIC DCMPLX (susan)
	//*
	//*     .. Scalar Arguments ..
	//*     .. Scalars in Common ..
	//*     .. Parameters ..
	(*one) = 1.0
	(*two) = 2.0
	//*     .. Local Scalars ..
	//*     .. Local Arrays ..
	//*     .. External Subroutines ..
	//*     .. Intrinsic Functions ..
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
	//*
	//*     Initialize ALPHA, BETA, RALPHA, and RBETA.
	//*
	(*alpha) = (*DCMPLX(one, -(*one)))
	(*beta) = (*DCMPLX(two, -(*two)))
	(*ralpha) = (*one)
	(*rbeta) = (*two)
	//*
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
	}
Label10:
	;
	(*infot) = 1
	Zgemm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 1
	Zgemm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 1
	Zgemm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 8
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 13
	Zgemm(func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label20:
	;
	(*infot) = 1
	Zhemm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zhemm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zhemm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label30:
	;
	(*infot) = 1
	Zsymm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsymm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsymm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label40:
	;
	(*infot) = 1
	Ztrmm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrmm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label50:
	;
	(*infot) = 1
	Ztrsm(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('/'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 5
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 6
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 11
	Ztrsm(func() *byte{y := byte('r'); return &y}(), func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label60:
	;
	(*infot) = 1
	Zherk(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 2; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zherk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 2; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zherk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), ralpha, a, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label70:
	;
	(*infot) = 1
	Zsyrk(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zsyrk(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 10
	Zsyrk(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label80:
	;
	(*infot) = 1
	Zher2k(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zher2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zher2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), rbeta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	goto Label100
Label90:
	;
	(*infot) = 1
	Zsyr2k(func() *byte{y := byte('/'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 2
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('c'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 3
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), -1, func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 4
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), -1, alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 7
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 2; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 9
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 0; return &y}(), func() *int{y := 2; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsyr2k(func() *byte{y := byte('u'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('n'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 2; return &y}(), b, func() *int{y := 2; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	(*infot) = 12
	Zsyr2k(func() *byte{y := byte('l'); return &y}(), func() *byte{y := byte('t'); return &y}(), func() *int{y := 2; return &y}(), func() *int{y := 0; return &y}(), alpha, a, func() *int{y := 1; return &y}(), b, func() *int{y := 1; return &y}(), beta, c, func() *int{y := 1; return &y}())
	CHKXER(srnamt, infot, nout, lerr, ok)
	//*
Label100:
	;
	if *ok {
		WRITE((*nout), *func() *[]byte{y := []byte(" %6s passed the tests of error-exits\n"); return &y}(), (*srnamt))
	} else {
		WRITE((*nout), *func() *[]byte{y := []byte(" ******* %6s failed the tests of error-exits *******\n"); return &y}(), (*srnamt))
	}
	return
	//*
	//*
	//*     End of ZCHKE.
	//*
}

func Zmake(_type *int, uplo *byte, diag *byte, m *int, n *int, a *[][]complex128, nmax *int, aa *[]complex128, lda *int, reset *bool, transl *complex128) {
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
	ibeg := new(int)
	iend := new(int)
	j := new(int)
	jj := new(int)
	gen := new(bool)
	her := new(bool)
	lower := new(bool)
	sym := new(bool)
	tri := new(bool)
	unit := new(bool)
	upper := new(bool)
	//*
	//*  Generates values for an M by N matrix A.
	//*  Stores the values in the array AA in the data structure required
	//*  by the routine, with unwanted elements set to rogue value.
	//*
	//*  TYPE is 'GE', 'HE', 'SY' or 'TR'.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	(*gen) = (*type) == "ge"
	(*her) = (*type) == "he"
	(*sym) = (*type) == "sy"
	(*tri) = (*type) == "tr"
	(*upper) = ((*her) || (*sym) || (*tri)) && (*uplo) == "u"
	(*lower) = ((*her) || (*sym) || (*tri)) && (*uplo) == "l"
	(*unit) = (*tri) && (*diag) == "u"
	//*
	//*     Generate data in array A.
	//*
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		for (*i) = 1; (*i) <= (*m); (*i)++ {
			if (*gen) || ((*upper) && (*i) <= (*j)) || ((*lower) && (*i) >= (*j)) {
				(*a)[(*i)-1][(*j)-1] = ZBEG(reset) + (*transl)
				if (*i) != (*j) {
					//*                 Set some elements to zero
					if (*n) > 3&& (*j) == (*n) / 2 {
						(*a)[(*i)-1][(*j)-1] = (*zero)
					}
					if *her {
						(*a)[(*j)-1][(*i)-1] = (DCONJG(((*a)[(*i)-1][(*j)-1])))
					} else if *sym {
						(*a)[(*j)-1][(*i)-1] = (*a)[(*i)-1][(*j)-1]
					} else if *tri {
						(*a)[(*j)-1][(*i)-1] = (*zero)
					}
				}
			}
		//Label10:
		}
		if *her {
			(*a)[(*j)-1][(*j)-1] = (*DCMPLX(DBLE(((*a)[(*j)-1][(*j)-1])), rzero))
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
	} else if (*type) == "he" || (*type) == "sy" || (*type) == "tr" {
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
			//Label60:
			}
			for (*i) = (*ibeg); (*i) <= (*iend); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*a)[(*i)-1][(*j)-1]
			//Label70:
			}
			for (*i) = (*iend) + 1; (*i) <= (*lda); (*i)++ {
				(*aa)[(*i)+((*j)-1)*(*lda)-1] = (*rogue)
			//Label80:
			}
			if *her {
				(*jj) = (*j) + ((*j)-1)*(*lda)
				(*aa)[(*jj)-1] = (*DCMPLX(DBLE(((*aa)[(*jj)-1])), rrogue))
			}
		//Label90:
		}
	}
	return
	//*
	//*     End of ZMAKE.
	//*
}

func Zmmch(transa *byte, transb *byte, m *int, n *int, kk *int, alpha *complex128, a *[][]complex128, lda *int, b *[][]complex128, ldb *int, beta *complex128, c *[][]complex128, ldc *int, ct *[]complex128, g *[]float64, cc *[][]complex128, ldcc *int, eps *float64, err *float64, fatal *bool, nout *int, mv *bool) {
	zero := new(complex128)
	rzero := new(float64)
	rone := new(float64)
	cl := new(complex128)
	erri := new(float64)
	i := new(int)
	j := new(int)
	k := new(int)
	ctrana := new(bool)
	ctranb := new(bool)
	trana := new(bool)
	tranb := new(bool)
	abs1 := new(float64)
	//*
	//*  Checks the results of the computational tests.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	ABS1(cl) = ABS(DBLE((*cl))) + ABS(DIMAG(cl))
	//*     .. Executable Statements ..
	(*trana) = (*transa) == "t" || (*transa) == "c"
	(*tranb) = (*transb) == "t" || (*transb) == "c"
	(*ctrana) = (*transa) == "c"
	(*ctranb) = (*transb) == "c"
	//*
	//*     Compute expected result, one column at a time, in CT using data
	//*     in A, B and C.
	//*     Compute gauges in G.
	//*
	for (*j) = 1; (*j) <= (*n); (*j)++ {
		//*
		for (*i) = 1; (*i) <= (*m); (*i)++ {
			(*ct)[(*i)-1] = (*zero)
			(*g)[(*i)-1] = (*rzero)
		//Label10:
		}
		if !(*trana) && !(*tranb) {
			for (*k) = 1; (*k) <= (*kk); (*k)++ {
				for (*i) = 1; (*i) <= (*m); (*i)++ {
					(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*i)-1][(*k)-1]*(*b)[(*k)-1][(*j)-1]
					(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*i)-1][(*k)-1]))*ABS1(&((*b)[(*k)-1][(*j)-1]))
				//Label20:
				}
			//Label30:
			}
		} else if (*trana) && !(*tranb) {
			if *ctrana {
				for (*k) = 1; (*k) <= (*kk); (*k)++ {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*ct)[(*i)-1] = (*ct)[(*i)-1] + DCONJG(((*a)[(*k)-1][(*i)-1]))*(*b)[(*k)-1][(*j)-1]
						(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*k)-1][(*j)-1]))
					//Label40:
					}
				//Label50:
				}
			} else {
				for (*k) = 1; (*k) <= (*kk); (*k)++ {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*k)-1][(*i)-1]*(*b)[(*k)-1][(*j)-1]
						(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*k)-1][(*j)-1]))
					//Label60:
					}
				//Label70:
				}
			}
		} else if !(*trana) && (*tranb) {
			if *ctranb {
				for (*k) = 1; (*k) <= (*kk); (*k)++ {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*i)-1][(*k)-1]*DCONJG(((*b)[(*j)-1][(*k)-1]))
						(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*i)-1][(*k)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
					//Label80:
					}
				//Label90:
				}
			} else {
				for (*k) = 1; (*k) <= (*kk); (*k)++ {
					for (*i) = 1; (*i) <= (*m); (*i)++ {
						(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*i)-1][(*k)-1]*(*b)[(*j)-1][(*k)-1]
						(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*i)-1][(*k)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
					//Label100:
					}
				//Label110:
				}
			}
		} else if (*trana) && (*tranb) {
			if *ctrana {
				if *ctranb {
					for (*k) = 1; (*k) <= (*kk); (*k)++ {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*ct)[(*i)-1] = (*ct)[(*i)-1] + DCONJG(((*a)[(*k)-1][(*i)-1]))*DCONJG(((*b)[(*j)-1][(*k)-1]))
							(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
						//Label120:
						}
					//Label130:
					}
				} else {
					for (*k) = 1; (*k) <= (*kk); (*k)++ {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*ct)[(*i)-1] = (*ct)[(*i)-1] + DCONJG(((*a)[(*k)-1][(*i)-1]))*(*b)[(*j)-1][(*k)-1]
							(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
						//Label140:
						}
					//Label150:
					}
				}
			} else {
				if *ctranb {
					for (*k) = 1; (*k) <= (*kk); (*k)++ {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*k)-1][(*i)-1]*DCONJG(((*b)[(*j)-1][(*k)-1]))
							(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
						//Label160:
						}
					//Label170:
					}
				} else {
					for (*k) = 1; (*k) <= (*kk); (*k)++ {
						for (*i) = 1; (*i) <= (*m); (*i)++ {
							(*ct)[(*i)-1] = (*ct)[(*i)-1] + (*a)[(*k)-1][(*i)-1]*(*b)[(*j)-1][(*k)-1]
							(*g)[(*i)-1] = (*g)[(*i)-1] + ABS1(&((*a)[(*k)-1][(*i)-1]))*ABS1(&((*b)[(*j)-1][(*k)-1]))
						//Label180:
						}
					//Label190:
					}
				}
			}
		}
		for (*i) = 1; (*i) <= (*m); (*i)++ {
			(*ct)[(*i)-1] = (*alpha)*(*ct)[(*i)-1] + (*beta)*(*c)[(*i)-1][(*j)-1]
			(*g)[(*i)-1] = ABS1(alpha)*(*g)[(*i)-1] + ABS1(beta)*ABS1(&((*c)[(*i)-1][(*j)-1]))
		//Label200:
		}
		//*
		//*        Compute the error ratio for this result.
		//*
		(*err) = (*zero)
		for (*i) = 1; (*i) <= (*m); (*i)++ {
			(*erri) = ABS1((*ct)[(*i)-1]-(*cc)[(*i)-1][(*j)-1]) / (*eps)
			if (*g)[(*i)-1] != (*rzero) {
				(*erri) = (*erri) / (*g)[(*i)-1]
			}
			(*err) = (MAX((*err), (*erri)))
			if (*err) * SQRT((*eps)) >= (*rone) {
				goto Label230
			}
		//Label210:
		}
		//*
	//Label220:
	}
	//*
	//*     If the loop completes, all results are at least half accurate.
	goto Label250
	//*
	//*     Report fatal error.
	//*
Label230:
	;
	(*fatal) = true
	WRITE((*nout), *func() *[]byte{y := []byte(" ******* fatal error - computed result is less than half accurate *******\n                       expected result                    computed result\n"); return &y}())
	for (*i) = 1; (*i) <= (*m); (*i)++ {
		if *mv {
			WRITE((*nout), *func() *[]byte{y := []byte(" %7d  (%15.6f,%15.6f)\n"); return &y}(), (*i), (*ct)[(*i)-1], (*cc)[(*i)-1][(*j)-1])
		} else {
			WRITE((*nout), *func() *[]byte{y := []byte(" %7d  (%15.6f,%15.6f)\n"); return &y}(), (*i), (*cc)[(*i)-1][(*j)-1], (*ct)[(*i)-1])
		}
	//Label240:
	}
	if (*n) > 1 {
		WRITE((*nout), *func() *[]byte{y := []byte("      these are the results for column %3d\n"); return &y}(), (*j))
	}
	//*
Label250:
	;
	return
	//*
	//*
	//*     End of ZMMCH.
	//*
}

func Lze(ri *[]complex128, rj *[]complex128, lr *int) (lzeReturn *bool) {
	lzereturn := new(bool)
	i := new(int)
	//*
	//*  Tests if two arrays are identical.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	(*lze) = true
	goto Label30
Label20:
	;
	(*lze) = false
Label30:
	;
	return
	//*
	//*     End of LZE.
	//*
}

func Lzeres(_type *int, uplo *byte, m *int, n *int, aa *[][]complex128, as *[][]complex128, lda *int) (lzeresReturn *bool) {
	lzeresreturn := new(bool)
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
	//*  TYPE is 'GE' or 'HE' or 'SY'.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	} else if (*type) == "he" || (*type) == "sy" {
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
	(*lzeres) = true
	goto Label80
Label70:
	;
	(*lzeres) = false
Label80:
	;
	return
	//*
	//*     End of LZERES.
	//*
}

func Zbeg(reset *bool) (zbegReturn *complex128) {
	zbegreturn := new(complex128)
	i := new(int)
	ic := new(int)
	j := new(int)
	mi := new(int)
	mj := new(int)
	//*
	//*  Generates complex numbers as pairs of random numbers uniformly
	//*  distributed between -0.5 and 0.5.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	(*zbeg) = (*DCMPLX(((*i)-500)/1001.0, ((*j)-500)/1001.0))
	return
	//*
	//*     End of ZBEG.
	//*
}

func Ddiff(x *float64, y *float64) (ddiffReturn *float64) {
	ddiffreturn := new(float64)
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
	//*
	//*     .. Scalar Arguments ..
	//*     .. Executable Statements ..
	(*ddiff) = (*x) - (*y)
	return
	//*
	//*     End of DDIFF.
	//*
}

func Chkxer(srnamt *[]byte, infot *int, nout *int, lerr *bool, ok *bool) {
	//*
	//*  Tests whether Xerbla has detected an error when it should.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
	//*  the test program for testing error exits from the Level 3 BLAS
	//*  routines.
	//*
	//*  Xerbla  is an error handler for the Level 3 BLAS routines.
	//*
	//*  It is called by the Level 3 BLAS routines if an input parameter is
	//*  invalid.
	//*
	//*  Auxiliary routine for test program for Level 3 Blas.
	//*
	//*  -- Written on 8-February-1989.
	//*     Jack Dongarra, Argonne National Laboratory.
	//*     Iain Duff, AERE Harwell.
	//*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
	//*     Sven Hammarling, Numerical Algorithms Group Ltd.
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
