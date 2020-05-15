package goblas

import (
	"time"

	"github.com/whipstein/golapack/blas"
)

// Dchkaa is the main test program for the DOUBLE PRECISION lapACK
// linear equation routines
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       PROGRAM Dchkaa
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Dchkaa is the main test program for the DOUBLE PRECISION lapACK
// linear equation routines
//
// The program must be driven by a short data file. The first 15 records
// (not including the first comment  line) specify problem dimensions
// and program options using list-directed input. The remaining lines
// specify the lapACK test paths and the number of matrix types to use
// in testing.  An annotated example of a data file can be obtained by
// deleting the first 3 characters from the following 40 lines:
// Data file for testing DOUBLE PRECISION lapACK linear eqn. routines
// 7                      Number of values of M
// 0 1 2 3 5 10 16        Values of M (row dimension)
// 7                      Number of values of N
// 0 1 2 3 5 10 16        Values of N (column dimension)
// 1                      Number of values of nrhs
// 2                      Values of nrhs (number of right hand sides)
// 5                      Number of values of nb
// 1 3 3 3 20             Values of nb (the blocksize)
// 1 0 5 9 1              Values of nx (crossover point)
// 3                      Number of values of rank
// 30 50 90               Values of rank (as a % of N)
// 20.0                   Threshold value of test ratio
// T                      Put T to test the lapACK routines
// T                      Put T to test the driver routines
// T                      Put T to test the error exits
// DGE   11               list types on next line if 0 < ntypes < 11
// DGB    8               list types on next line if 0 < ntypes <  8
// DGT   12               list types on next line if 0 < ntypes < 12
// DPO    9               list types on next line if 0 < ntypes <  9
// DPS    9               list types on next line if 0 < ntypes <  9
// DPP    9               list types on next line if 0 < ntypes <  9
// DPB    8               list types on next line if 0 < ntypes <  8
// DPT   12               list types on next line if 0 < ntypes < 12
// DSY   10               list types on next line if 0 < ntypes < 10
// DSR   10               list types on next line if 0 < ntypes < 10
// DSK   10               list types on next line if 0 < ntypes < 10
// DSA   10               list types on next line if 0 < ntypes < 10
// Ds2   10               list types on next line if 0 < ntypes < 10
// DSP   10               list types on next line if 0 < ntypes < 10
// DTR   18               list types on next line if 0 < ntypes < 18
// DTP   18               list types on next line if 0 < ntypes < 18
// DTB   17               list types on next line if 0 < ntypes < 17
// DQR    8               list types on next line if 0 < ntypes <  8
// DRQ    8               list types on next line if 0 < ntypes <  8
// DLQ    8               list types on next line if 0 < ntypes <  8
// DQL    8               list types on next line if 0 < ntypes <  8
// DQP    6               list types on next line if 0 < ntypes <  6
// DTZ    3               list types on next line if 0 < ntypes <  3
// DLS    6               list types on next line if 0 < ntypes <  6
// DEQ
// DQT
// DQX
// DTQ
// DXQ
// Dts
// DHH
// \endverbatim
//
//  Parameters:
//  ==========
//
// \verbatim
//  nmax    intEGER
//          The maximum allowable value for M and N.
//
//  maxin   intEGER
//          The number of different values that can be used for each of
//          m, n, nrhs, nb, nx and rank
//
//  maxrhs  intEGER
//          The maximum number of right hand sides
//
//  matmax  intEGER
//          The maximum number of matrix types to use for testing
//
//  nin     intEGER
//          The unit number for input
//
//  nout    intEGER
//          The unit number for output
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
// \date November 2019
//
// \ingroup double_lin
//
//  =====================================================================
func Dchkaa() {
	start := time.Now()
	nmax := new(int)
	maxin := new(int)
	maxrhs := new(int)
	matmax := new(int)
	nin := new(int)
	nout := new(int)
	kdmax := new(int)
	fatal := new(bool)
	tstchk := new(bool)
	tstdrv := new(bool)
	tsterr := new(bool)
	c1 := new(byte)
	c2 := func() *[]byte {
		arr := make([]byte, 2)
		return &arr
	}()
	path := func() *[]byte {
		arr := make([]byte, 3)
		return &arr
	}()
	intstr := func() *[]byte {
		arr := make([]byte, 10)
		return &arr
	}()
	aline := func() *[]byte {
		arr := make([]byte, 72)
		return &arr
	}()
	i := new(int)
	ic := new(int)
	j := new(int)
	k := new(int)
	la := new(int)
	lafac := new(int)
	lda := new(int)
	nb := new(int)
	nm := new(int)
	nmats := new(int)
	nn := new(int)
	nnb := new(int)
	nnb2 := new(int)
	nns := new(int)
	nrhs := new(int)
	ntypes := new(int)
	nrank := new(int)
	versMajor := new(int)
	versMinor := new(int)
	versPatch := new(int)
	eps := new(float64)
	s1 := new(float64)
	s2 := new(float64)
	threq := new(float64)
	thresh := new(float64)
	dotype := func() *[]bool {
		arr := make([]bool, 30)
		return &arr
	}()
	iwork := func() *[]int {
		arr := make([]int, -1)
		return &arr
	}()
	mval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	nbval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	nbval2 := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	nsval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	nval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	nxval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	rankval := func() *[]int {
		arr := make([]int, 12)
		return &arr
	}()
	piv := func() *[]int {
		arr := make([]int, 132)
		return &arr
	}()
	a := func() *[][]float64 {
		arr := make([][]float64, -1)
		for u := 0; u < -1; u++ {
			arr[u] = make([]float64, 7)
		}
		return &arr
	}()
	b := func() *[][]float64 {
		arr := make([][]float64, -1)
		for u := 0; u < -1; u++ {
			arr[u] = make([]float64, 4)
		}
		return &arr
	}()
	e := func() *[]float64 {
		arr := make([]float64, 132)
		return &arr
	}()
	rwork := func() *[]float64 {
		arr := make([]float64, -1)
		return &arr
	}()
	s := func() *[]float64 {
		arr := make([]float64, -1)
		return &arr
	}()
	work := func() *[][]float64 {
		arr := make([][]float64, 132)
		for u := 0; u < 132; u++ {
			arr[u] = make([]float64, -1)
		}
		return &arr
	}()
	lerr := new(bool)
	ok := new(bool)
	srnamt := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	infot := new(int)
	nunit := new(int)
	iparms := func() *[]int {
		arr := make([]int, 100)
		return &arr
	}()
	common.infoc.lerr = new(bool)
	common.infoc.ok = new(bool)
	common.infoc.nunit = new(int)
	common.infoc.infot = new(int)
	common.srnamc.srnamt = new([]byte)
	common.claenv.iparms = new([]int)
	//
	//  -- lapACK test routine (version 3.9.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     Novemebr 2019
	//
	//  =====================================================================
	//
	//     .. Parameters ..
	(*nmax) = 132
	(*maxin) = 12
	(*maxrhs) = 16
	(*matmax) = 30
	(*nin) = 5
	(*nout) = 6
	(*kdmax) = (*nmax) + ((*nmax)+1)/4
	//     ..
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. External Subroutines ..
	//     ..
	//     .. Scalars in common ..
	//     ..
	//     .. Arrays in common ..
	//     ..
	//     .. common blocks ..
	infot = common.infoc.infot
	nunit = common.infoc.nunit
	ok = common.infoc.ok
	lerr = common.infoc.lerr
	srnamt = common.srnamc.srnamt
	iparms = common.claenv.iparms
	//     ..
	//     .. Data statements ..
	(*threq), (*intstr)[0], (*intstr)[1], (*intstr)[2], (*intstr)[3], (*intstr)[4], (*intstr)[5], (*intstr)[6], (*intstr)[7], (*intstr)[8], (*intstr)[9] = 2.0, '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'
	//     ..
	//     .. Executable Statements ..
	//
	(*s1) = (*Dsecnd(&start))
	(*lda) = (*nmax)
	(*fatal) = false
	//
	//     Read a dummy line.
	//
	REAd(nin, []byte(" %v\n"))
	//
	//     Report values of parameters.
	//
	Ilaver(versMajor, versMinor, versPatch)
	WRITE((*nout), *func() *[]byte {
		y := []byte(" Tests of the DOUBLE PRECISION lapACK routines \n lapACK VERSION %1d.%1d.%1d\n\n The following parameter values will be used:\n")
		return &y
	}(), (*versMajor), (*versMinor), (*versPatch))
	//
	//     Read the values of M
	//
	REAd(nin, []byte(" %v\n"), nm)
	if (*nm) < 1 {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nm "); return &y }(), (*nm), 1)
		(*nm) = 0
		(*fatal) = true
	} else if (*nm) > (*maxin) {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nm "); return &y }(), (*nm), (*maxin))
		(*nm) = 0
		(*fatal) = true
	}
	WRITE((*nin), []byte(" %v\n"), (*mval)[(*i)-1], (*i) == 1, (*nm))
	for (*i) = 1; (*i) <= (*nm); (*i)++ {
		if (*mval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" M  "); return &y }(), (*mval)[(*i)-1], 0)
			(*fatal) = true
		} else if (*mval)[(*i)-1] > (*nmax) {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" M  "); return &y }(), (*mval)[(*i)-1], (*nmax))
			(*fatal) = true
		}
		//Label10:
	}
	if (*nm) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "M   ", (*mval)[(*i)-1], (*i) == 1, (*nm))
	}
	//
	//     Read the values of N
	//
	REAd(nin, []byte(" %v\n"), nn)
	if (*nn) < 1 {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nn "); return &y }(), (*nn), 1)
		(*nn) = 0
		(*fatal) = true
	} else if (*nn) > (*maxin) {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nn "); return &y }(), (*nn), (*maxin))
		(*nn) = 0
		(*fatal) = true
	}
	WRITE((*nin), []byte(" %v\n"), (*nval)[(*i)-1], (*i) == 1, (*nn))
	for (*i) = 1; (*i) <= (*nn); (*i)++ {
		if (*nval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" N  "); return &y }(), (*nval)[(*i)-1], 0)
			(*fatal) = true
		} else if (*nval)[(*i)-1] > (*nmax) {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" N  "); return &y }(), (*nval)[(*i)-1], (*nmax))
			(*fatal) = true
		}
		//Label20:
	}
	if (*nn) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "N   ", (*nval)[(*i)-1], (*i) == 1, (*nn))
	}
	//
	//     Read the values of nrhs
	//
	REAd(nin, []byte(" %v\n"), nns)
	if (*nns) < 1 {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nns"); return &y }(), (*nns), 1)
		(*nns) = 0
		(*fatal) = true
	} else if (*nns) > (*maxin) {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nns"); return &y }(), (*nns), (*maxin))
		(*nns) = 0
		(*fatal) = true
	}
	WRITE((*nin), []byte(" %v\n"), (*nsval)[(*i)-1], (*i) == 1, (*nns))
	for (*i) = 1; (*i) <= (*nns); (*i)++ {
		if (*nsval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte("nrhs"); return &y }(), (*nsval)[(*i)-1], 0)
			(*fatal) = true
		} else if (*nsval)[(*i)-1] > (*maxrhs) {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte("nrhs"); return &y }(), (*nsval)[(*i)-1], (*maxrhs))
			(*fatal) = true
		}
		//Label30:
	}
	if (*nns) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "nrhs", (*nsval)[(*i)-1], (*i) == 1, (*nns))
	}
	//
	//     Read the values of nb
	//
	REAd(nin, []byte(" %v\n"), nnb)
	if (*nnb) < 1 {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte("nnb "); return &y }(), (*nnb), 1)
		(*nnb) = 0
		(*fatal) = true
	} else if (*nnb) > (*maxin) {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte("nnb "); return &y }(), (*nnb), (*maxin))
		(*nnb) = 0
		(*fatal) = true
	}
	WRITE((*nin), []byte(" %v\n"), (*nbval)[(*i)-1], (*i) == 1, (*nnb))
	for (*i) = 1; (*i) <= (*nnb); (*i)++ {
		if (*nbval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nb "); return &y }(), (*nbval)[(*i)-1], 0)
			(*fatal) = true
		}
		//Label40:
	}
	if (*nnb) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "nb  ", (*nbval)[(*i)-1], (*i) == 1, (*nnb))
	}
	//
	//     Set nbval2 to be the set of unique values of nb
	//
	(*nnb2) = 0
	for (*i) = 1; (*i) <= (*nnb); (*i)++ {
		(*nb) = (*nbval)[(*i)-1]
		for (*j) = 1; (*j) <= (*nnb2); (*j)++ {
			if (*nb) == (*nbval2)[(*j)-1] {
				goto Label60
			}
			//Label50:
		}
		(*nnb2) = (*nnb2) + 1
		(*nbval2)[(*nnb2)-1] = (*nb)
	Label60:
	}
	//
	//     Read the values of nx
	//
	WRITE((*nin), []byte(" %v\n"), (*nxval)[(*i)-1], (*i) == 1, (*nnb))
	for (*i) = 1; (*i) <= (*nnb); (*i)++ {
		if (*nxval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nx "); return &y }(), (*nxval)[(*i)-1], 0)
			(*fatal) = true
		}
		//Label70:
	}
	if (*nnb) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "nx  ", (*nxval)[(*i)-1], (*i) == 1, (*nnb))
	}
	//
	//     Read the values of rankval
	//
	REAd(nin, []byte(" %v\n"), nrank)
	if (*nn) < 1 {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nrank "); return &y }(), (*nrank), 1)
		(*nrank) = 0
		(*fatal) = true
	} else if (*nn) > (*maxin) {
		WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" nrank "); return &y }(), (*nrank), (*maxin))
		(*nrank) = 0
		(*fatal) = true
	}
	WRITE((*nin), []byte(" %v\n"), (*rankval)[(*i)-1], (*i) == 1, (*nrank))
	for (*i) = 1; (*i) <= (*nrank); (*i)++ {
		if (*rankval)[(*i)-1] < 0 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be >=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" rank  "); return &y }(), (*rankval)[(*i)-1], 0)
			(*fatal) = true
		} else if (*rankval)[(*i)-1] > 100 {
			WRITE((*nout), *func() *[]byte {y := []byte(" Invalid input value: %4s=%6d; must be <=%6d\n"); return &y }(), *func() *[]byte {y := []byte(" rank  "); return &y }(), (*rankval)[(*i)-1], 100)
			(*fatal) = true
		}
	}
	if (*nrank) > 0 {
		WRITE((*nout), []byte("    %4s:  %6d\n           %6d\n"), "rank % OF N", (*rankval)[(*i)-1], (*i) == 1, (*nrank))
	}
	//
	//     Read the threshold value for the test ratios.
	//
	REAd(nin, []byte(" %v\n"), thresh)
	WRITE((*nout), *func() *[]byte {
		y := []byte("\n Routines pass computational tests if test ratio is less than%8.2f\n\n")
		return &y
	}(), (*thresh))
	//
	//     Read the flag that indicates whether to test the lapACK routines.
	//
	REAd(nin, []byte(" %v\n"), tstchk)
	//
	//     Read the flag that indicates whether to test the driver routines.
	//
	REAd(nin, []byte(" %v\n"), tstdrv)
	//
	//     Read the flag that indicates whether to test the error exits.
	//
	REAd(nin, []byte(" %v\n"), tsterr)
	//
	if *fatal {
		WRITE((*nout), *func() *[]byte {y := []byte("\n Execution not attempted due to input errors\n"); return &y }())
		panic("")
	}
	//
	//     Calculate and print the machine dependent _constants.
	//
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("Underflow threshold"); return &y }()))
	WRITE((*nout), *func() *[]byte {y := []byte(" Relative machine %s is taken to be%16.6E\n"); return &y }(), *func() *[]byte {y := []byte("underflow"); return &y }(), (*eps))
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("Overflow threshold"); return &y }()))
	WRITE((*nout), *func() *[]byte {y := []byte(" Relative machine %s is taken to be%16.6E\n"); return &y }(), *func() *[]byte {y := []byte("overflow "); return &y }(), (*eps))
	(*eps) = (*Dlamch(func() *[]byte {y := []byte("epsilon"); return &y }()))
	WRITE((*nout), *func() *[]byte {y := []byte(" Relative machine %s is taken to be%16.6E\n"); return &y }(), *func() *[]byte {y := []byte("precision"); return &y }(), (*eps))
	WRITE((*nout), *func() *[]byte {y := []byte(" %v\n"); return &y }())
	//
Label80:
	;
	//
	//     Read a test path and the number of matrix types to use.
	//
	REAd(nin, []byte("%72s\n"), aline)
	(*path)[0] = (*aline)[0]
	(*nmats) = (*matmax)
	(*i) = 3
Label90:
	;
	(*i) = (*i) + 1
	if (*i) > 72 {
		(*nmats) = (*matmax)
		goto Label130
	}
	if (*aline)[(*i)-1] == ' ' {
		goto Label90
	}
	(*nmats) = 0
Label100:
	;
	(*c1) = (*aline)[(*i)-1]
	for (*k) = 1; (*k) <= 10; (*k)++ {
		if (*c1) == (*intstr)[(*k)-1] {
			(*ic) = (*k) - 1
			goto Label120
		}
		//Label110:
	}
	goto Label130
Label120:
	;
	(*nmats) = (*nmats)*10 + (*ic)
	(*i) = (*i) + 1
	if (*i) > 72 {
		goto Label130
	}
	goto Label100
Label130:
	;
	(*c1) = (*path)[0]
	(*c2)[0] = (*path)[1]
	(*nrhs) = (*nsval)[0]
	//
	//     Check first character for correct precision.
	//
	if !blas.Lsame(c1, func() *[]byte {y := []byte("Double precision"); return &y }()) {
		WRITE((*nout), *func() *[]byte {y := []byte("\n %3s:  Unrecognized path name\n"); return &y }(), (*path))
		//
	} else if (*nmats) <= 0 {
		//
		//        Check for a positive number of tests requested.
		//
		WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("GE"); return &y }()) {
		//
		//        GE:  general matrices
		//
		(*ntypes) = 11
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkge(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvge(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), s, work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("GB"); return &y }()) {
		//
		//        GB:  general banded matrices
		//
		(*la) = (2*(*kdmax) + 1) * (*nmax)
		(*lafac) = (3*(*kdmax) + 1) * (*nmax)
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkgb(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, &((*a)[0][0]), la, &((*a)[0][2]), lafac, &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvgb(dotype, nn, nval, nrhs, thresh, tsterr, &((*a)[0][0]), la, &((*a)[0][2]), lafac, &((*a)[0][5]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), s, work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("GT"); return &y }()) {
		//
		//        GT:  general tridiagonal matrices
		//
		(*ntypes) = 12
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkgt(dotype, nn, nval, nns, nsval, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvgt(dotype, nn, nval, nrhs, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PO"); return &y }()) {
		//
		//        PO:  positive definite matrices
		//
		(*ntypes) = 9
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkpo(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvpo(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), s, work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PS"); return &y }()) {
		//
		//        PS:  positive semi-definite matrices
		//
		(*ntypes) = 9
		//
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkps(dotype, nn, nval, nnb2, nbval2, nrank, rankval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), piv, work, rwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PP"); return &y }()) {
		//
		//        PP:  positive definite packed matrices
		//
		(*ntypes) = 9
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkpp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvpp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), s, work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PB"); return &y }()) {
		//
		//        PB:  positive definite banded matrices
		//
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkpb(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvpb(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), s, work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("PT"); return &y }()) {
		//
		//        PT:  positive definite tridiagonal matrices
		//
		(*ntypes) = 12
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkpt(dotype, nn, nval, nns, nsval, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvpt(dotype, nn, nval, nrhs, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SY"); return &y }()) {
		//
		//        SY:  symmetric indefinite matrices,
		//             with partial (Bunch-Kaufman) pivoting algorithm
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchksy(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvsy(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SR"); return &y }()) {
		//
		//        SR:  symmetric indefinite matrices,
		//             with bounded Bunch-Kaufman (rook) pivoting algorithm
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			DchksyRook(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			DdrvsyRook(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SK"); return &y }()) {
		//
		//        SK:  symmetric indefinite matrices,
		//             with bounded Bunch-Kaufman (rook) pivoting algorithm,
		//             differnet matrix storage format than SR path version.
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			DchksyRk(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), e, &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			DdrvsyRk(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), e, &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SA"); return &y }()) {
		//
		//        SA:  symmetric indefinite matrices,
		//             with partial (Aasen's) pivoting algorithm
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			DchksyAa(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			DdrvsyAa(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("s2"); return &y }()) {
		//
		//        SA:  symmetric indefinite matrices,
		//             with partial (Aasen's) pivoting algorithm
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			DchksyAa2stage(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			DdrvsyAa2stage(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("SP"); return &y }()) {
		//
		//        SP:  symmetric indefinite packed matrices,
		//             with partial (Bunch-Kaufman) pivoting algorithm
		//
		(*ntypes) = 10
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchksp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
		if *tstdrv {
			Ddrvsp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TR"); return &y }()) {
		//
		//        TR:  triangular matrices
		//
		(*ntypes) = 18
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchktr(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TP"); return &y }()) {
		//
		//        TP:  triangular packed matrices
		//
		(*ntypes) = 18
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchktp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TB"); return &y }()) {
		//
		//        TB:  triangular banded matrices
		//
		(*ntypes) = 17
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchktb(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QR"); return &y }()) {
		//
		//        QR:  QR factorization
		//
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkqr(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*a)[0][3]), &((*a)[0][4]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("LQ"); return &y }()) {
		//
		//        LQ:  LQ factorization
		//
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchklq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*a)[0][3]), &((*a)[0][4]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), work, rwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QL"); return &y }()) {
		//
		//        QL:  QL factorization
		//
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkql(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*a)[0][3]), &((*a)[0][4]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), work, rwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("RQ"); return &y }()) {
		//
		//        RQ:  RQ factorization
		//
		(*ntypes) = 8
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkrq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &((*a)[0][0]), &((*a)[0][1]), &((*a)[0][2]), &((*a)[0][3]), &((*a)[0][4]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), &((*b)[0][3]), work, rwork, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QP"); return &y }()) {
		//
		//        QP:  QR factorization with pivoting
		//
		(*ntypes) = 6
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchkq3(dotype, nm, mval, nn, nval, nnb, nbval, nxval, thresh, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][2]), work, iwork, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TZ"); return &y }()) {
		//
		//        TZ:  Trapezoidal matrix
		//
		(*ntypes) = 3
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstchk {
			Dchktz(dotype, nm, mval, nn, nval, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][2]), work, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("LS"); return &y }()) {
		//
		//        LS:  Least squares drivers
		//
		(*ntypes) = 6
		Alareq(path, nmats, dotype, ntypes, nin, nout)
		//
		if *tstdrv {
			Ddrvls(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, tsterr, &((*a)[0][0]), &((*a)[0][1]), &((*b)[0][0]), &((*b)[0][1]), &((*b)[0][2]), rwork, &((*rwork)[(*nmax)+0]), nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s driver routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("EQ"); return &y }()) {
		//
		//        EQ:  Equilibration routines for general and positive definite
		//             matrices (threq should be between 2 and 10)
		//
		if *tstchk {
			Dchkeq(threq, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QT"); return &y }()) {
		//
		//        QT:  QRT routines for general matrices
		//
		if *tstchk {
			Dchkqrt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("QX"); return &y }()) {
		//
		//        QX:  QRT routines for triangular-pentagonal matrices
		//
		if *tstchk {
			Dchkqrtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("TQ"); return &y }()) {
		//
		//        TQ:  LQT routines for general matrices
		//
		if *tstchk {
			Dchklqt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("XQ"); return &y }()) {
		//
		//        XQ:  LQT routines for triangular-pentagonal matrices
		//
		if *tstchk {
			Dchklqtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("ts"); return &y }()) {
		//
		//        ts:  QR routines for tall-skinny matrices
		//
		if *tstchk {
			Dchktsqr(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else if Lsamen(func() *int {y := 2; return &y }(), c2, func() *[]byte {y := []byte("HH"); return &y }()) {
		//
		//        HH:  Householder re_construction for tall-skinny matrices
		//
		if *tstchk {
			DchkorhrCol(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout)
		} else {
			WRITE((*nout), *func() *[]byte {y := []byte("\n %3s routines were not tested\n"); return &y }(), (*path))
		}
		//
	} else {
		//
		WRITE((*nout), *func() *[]byte {y := []byte("\n %3s:  Unrecognized path name\n"); return &y }(), (*path))
	}
	//
	//     Go back to get another input line.
	//
	goto Label80
	//
	//     Branch to this line when the last record is read.
	//
	//Label140:

	CLOSE(nin)
	(*s2) = (*Dsecnd())
	WRITE((*nout), *func() *[]byte {y := []byte("\n End of tests\n"); return &y }())
	WRITE((*nout), *func() *[]byte {y := []byte(" Total time used = %12.2f seconds\n\n"); return &y }(), (*s2)-(*s1))
	//
	//
	//     End of Dchkaa
	//
}
