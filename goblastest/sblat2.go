package main

import (
	"fmt"
	"strings"
)

type memory struct {
	combla struct {
		icase int
		n     int
		incx  int
		incy  int
		pass  bool
	}
	infoc struct {
		infot int
		nout  int
		ok    bool
		lerr  bool
		noutc int
	}
	srnamc struct {
		srnamt string
	}
}

var common memory

// \brief \b SBLAT2
//
// =========== DOCUMENTATION ===========
//
//Online html documentation available at
//           http://www.netlib.org/lapack/explore-html/
//
// Definition:
// ===========
//
//      PROGRAM SBLAT2
//
//
// \par Purpose:
// =============
//
// \verbatim
//
// Test program for the REAL Level 2 Blas.
//
// The program must be driven by a short data file. The first 18 records
// of the file are read using list-directed input, the last 16 records
// are read using the format ( A6, L2 ). An annotated example of a data
// file can be obtained by deleting the first 3 characters from the
// following 34 lines:
// 'sblat2.out'      NAME OF SUMMARY OUTPUT FILE
// 6                 UNIT NUMBER OF SUMMARY FILE
// 'SBLAT2.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
// -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
// F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
// F        LOGICAL FLAG, T TO STOP ON FAILURES.
// T        LOGICAL FLAG, T TO TEST ERROR EXITS.
// 16.0     THRESHOLD VALUE OF TEST RATIO
// 6                 NUMBER OF VALUES OF N
// 0 1 2 3 5 9       VALUES OF N
// 4                 NUMBER OF VALUES OF K
// 0 1 2 4           VALUES OF K
// 4                 NUMBER OF VALUES OF INCX AND INCY
// 1 2 -1 -2         VALUES OF INCX AND INCY
// 3                 NUMBER OF VALUES OF ALPHA
// 0.0 1.0 0.7       VALUES OF ALPHA
// 3                 NUMBER OF VALUES OF BETA
// 0.0 1.0 0.9       VALUES OF BETA
// SGEMV  T PUT F FOR NO TEST. SAME COLUMNS.
// SGBMV  T PUT F FOR NO TEST. SAME COLUMNS.
// SSYMV  T PUT F FOR NO TEST. SAME COLUMNS.
// SSBMV  T PUT F FOR NO TEST. SAME COLUMNS.
// SSPMV  T PUT F FOR NO TEST. SAME COLUMNS.
// STRMV  T PUT F FOR NO TEST. SAME COLUMNS.
// STBMV  T PUT F FOR NO TEST. SAME COLUMNS.
// STPMV  T PUT F FOR NO TEST. SAME COLUMNS.
// STRSV  T PUT F FOR NO TEST. SAME COLUMNS.
// STBSV  T PUT F FOR NO TEST. SAME COLUMNS.
// STPSV  T PUT F FOR NO TEST. SAME COLUMNS.
// SGER   T PUT F FOR NO TEST. SAME COLUMNS.
// SSYR   T PUT F FOR NO TEST. SAME COLUMNS.
// SSPR   T PUT F FOR NO TEST. SAME COLUMNS.
// SSYR2  T PUT F FOR NO TEST. SAME COLUMNS.
// SSPR2  T PUT F FOR NO TEST. SAME COLUMNS.
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
// Authors:
// ========
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
// =====================================================================
// VERIFIED
func main() {
	var _err error
	var nin, ntra, nsubs, nmax, incmax, ninmax, nidmax, nkbmax, nalmax, nbemax int
	var eps, err, thresh float32
	var i, isnum, j, n, nalf, nbet, nidim, ninc, nkb int
	var fatal, ltestt, rewi, same, sfatal, trace, tsterr bool
	var trans byte
	var line, snamet, snaps, summry string
	//
	// -- Reference BLAS test routine (version 3.7.0) --
	// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//    April 2012
	//
	// =====================================================================
	//
	nin = 5
	nsubs = 16
	nmax = 65
	incmax = 2
	ninmax = 7
	nidmax = 9
	nkbmax = 7
	nalmax = 7
	nbemax = 7
	snames := []string{
		"SGEMV",
		"SGBMV",
		"SSYMV",
		"SSBMV",
		"SSPMV",
		"STRMV",
		"STBMV",
		"STPMV",
		"STRSV",
		"STBSV",
		"STPSV",
		"SGER",
		"SSYR",
		"SSPR",
		"SSYR2",
		"SSPR2",
	}
	a := func() [][]float32 {
		arr := make([][]float32, nmax)
		for u := 0; u < nmax; u++ {
			arr[u] = make([]float32, nmax)
		}
		return arr
	}()
	aa := make([]float32, nmax*nmax)
	alf := make([]float32, nalmax)
	as := make([]float32, nmax*nmax)
	bet := make([]float32, nbemax)
	g := make([]float32, nmax)
	x := make([]float32, nmax)
	xs := make([]float32, nmax*incmax)
	xx := make([]float32, nmax*incmax)
	y := make([]float32, nmax)
	ys := make([]float32, nmax*incmax)
	yt := make([]float32, nmax)
	yy := make([]float32, nmax*incmax)
	// z := make([]float32, 2*nmax)
	idim := make([]int, nidmax)
	inc := make([]int, ninmax)
	kb := make([]int, nkbmax)
	ltest := make([]bool, nsubs)
	//    .. Executable Statements ..
	//
	//    Read name and unit number for summary output file and open file.
	//
	openRead(nin, "sblat2.in")
	defer units[nin].Close()

	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &summry)
	chkerr(_err)
	summry = strings.ReplaceAll(string(summry), "'", "")
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", common.infoc.nout)
	chkerr(_err)

	openWrite(common.infoc.nout, summry)
	defer units[common.infoc.nout].Close()
	//
	//    Read name and unit number for snapshot output file and open file.
	//
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &snaps)
	chkerr(_err)
	snaps = strings.ReplaceAll(snaps, "'", "")
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &ntra)
	chkerr(_err)

	trace = ntra >= 0
	if trace {
		openWrite(ntra, snaps)
		defer units[ntra].Close()
	}
	//    Read the flag that directs rewinding of the snapshot file.
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &rewi)
	chkerr(_err)
	rewi = rewi && trace
	//    Read the flag that directs stopping on any failure.
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &sfatal)
	chkerr(_err)
	//    Read the flag that indicates whether error exits are to be tested.
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &tsterr)
	chkerr(_err)
	//    Read the threshold value of the test ratio
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &thresh)
	chkerr(_err)
	//
	//    Read and check the parameter values for the tests.
	//
	//    Values of N
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &nidim)
	chkerr(_err)
	if nidim < 1 || nidim > nidmax {
		writeString(common.infoc.nout, " NUMBER OF VALUES OF %s IS LESS THAN 1 OR GREATER THAN %2d\n", 'N', nidmax)
		goto Label230
	}
	for i = 1; i <= nidim; i++ {
		_, _err = fmt.Sscanf(readval(readers[nin]), "%v", &idim[i-1])
		chkerr(_err)
	}
	_ = readln(readers[nin])
	for i = 1; i <= nidim; i++ {
		if idim[i-1] < 0 || idim[i-1] > nmax {
			writeString(common.infoc.nout, " VALUE OF N IS LESS THAN 0 OR GREATER THAN %2d\n", nmax)
			goto Label230
		}
		//Label10:
	}
	//    Values of K
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &nkb)
	chkerr(_err)
	if nkb < 1 || nkb > nkbmax {
		writeString(common.infoc.nout, " NUMBER OF VALUES OF %s IS LESS THAN 1 OR GREATER THAN %2d\n", 'K', nkbmax)
		goto Label230
	}
	for i = 1; i <= nkb; i++ {
		_, _err = fmt.Sscanf(readval(readers[nin]), "%v", &kb[i-1])
		chkerr(_err)
	}
	_ = readln(readers[nin])
	for i = 1; i <= nkb; i++ {
		if kb[i-1] < 0 {
			writeString(common.infoc.nout, " VALUE OF K IS LESS THAN 0\n")
			goto Label230
		}
		//Label20:
	}
	//    Values of INCX and INCY
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &ninc)
	chkerr(_err)
	if ninc < 1 || ninc > ninmax {
		writeString(common.infoc.nout, " NUMBER OF VALUES OF %s IS LESS THAN 1 OR GREATER THAN %2d\n", "INCX AND INCY", ninmax)
		goto Label230
	}
	for i = 1; i <= ninc; i++ {
		_, _err = fmt.Sscanf(readval(readers[nin]), "%v", &inc[i-1])
		chkerr(_err)
	}
	_ = readln(readers[nin])
	for i = 1; i <= ninc; i++ {
		if inc[i-1] == 0 || absint(inc[i-1]) > incmax {
			writeString(common.infoc.nout, " ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN %2d\n", incmax)
			goto Label230
		}
		//Label30:
	}
	//    Values of ALPHA
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &nalf)
	chkerr(_err)
	if nalf < 1 || nalf > nalmax {
		writeString(common.infoc.nout, " NUMBER OF VALUES OF %s IS LESS THAN 1 OR GREATER THAN %2d\n", "ALPHA", nalmax)
		goto Label230
	}
	for i = 1; i <= nalf; i++ {
		_, _err = fmt.Sscanf(readval(readers[nin]), "%v", &alf[i-1])
		chkerr(_err)
	}
	_ = readln(readers[nin])
	//    Values of BETA
	_, _err = fmt.Sscanf(readln(readers[nin]), "%v", &nbet)
	chkerr(_err)
	if nbet < 1 || nbet > nbemax {
		writeString(common.infoc.nout, " NUMBER OF VALUES OF %s IS LESS THAN 1 OR GREATER THAN %2d\n", "BETA", nbemax)
		goto Label230
	}
	for i = 1; i <= nbet; i++ {
		_, _err = fmt.Sscanf(readval(readers[nin]), "%v", &bet[i-1])
		chkerr(_err)
	}
	_ = readln(readers[nin])
	//
	//    Report values of parameters.
	//
	writeString(common.infoc.nout, " TESTS OF THE REAL             LEVEL 2 BLAS\n THE FOLLOWING PARAMETER VALUES WILL BE USED:\n")
	for i = 1; i <= nidim; i++ {
		writeString(common.infoc.nout, "   FOR N              %6d\n", idim[i-1])
	}
	for i = 1; i <= nkb; i++ {
		writeString(common.infoc.nout, "   FOR K              %6d\n", kb[i-1])
	}
	for i = 1; i <= ninc; i++ {
		writeString(common.infoc.nout, "   FOR INCX AND INCY  %6d\n", inc[i-1])
	}
	for i = 1; i <= nalf; i++ {
		writeString(common.infoc.nout, "   FOR ALPHA          %6.1f\n", alf[i-1])
	}
	for i = 1; i <= nbet; i++ {
		writeString(common.infoc.nout, "   FOR BETA           %6.1f\n", bet[i-1])
	}
	if !tsterr {
		writeString(common.infoc.nout, "\n")
		writeString(common.infoc.nout, " ERROR-EXITS WILL NOT BE TESTED\n")
	}
	writeString(common.infoc.nout, "\n")
	writeString(common.infoc.nout, " ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LESS THAN %8.2f\n", thresh)
	writeString(common.infoc.nout, "\n")
	//
	//    Read names of subroutines and flags which indicate
	//    whether they are to be tested.
	//
	// for i = 1; i <= nsubs; i++ {
	// 	ltest[i-1] = false
	// 	//Label40:
	// }
Label50:
	;
	line = readln(readers[nin])
	if line == "" {
		goto Label80
	}
	_, _err = fmt.Sscanf(line, "%6s %t", &snamet, &ltestt)
	chkerr(_err)
	snamet = strings.ToUpper(snamet)
	for i = 1; i <= nsubs; i++ {
		if snamet == snames[i-1] {
			goto Label70
		}
	}
	writeString(common.infoc.nout, " SUBPROGRAM NAME %6s NOT RECOGNIZED\n ******* TESTS ABANDONED *******\n", snamet)
	panic("")
Label70:
	;
	ltest[i-1] = ltestt
	goto Label50
	//
Label80:
	;
	_close(nin)
	//
	//    Compute EPS (the machine precision).
	//
	eps = epsilonf32()
	writeString(common.infoc.nout, " RELATIVE MACHINE PRECISION IS TAKEN TO BE %9.1e\n", eps)
	//
	//    Check the reliability of SMVCH using exact data.
	//
	n = min(32, nmax)
	for j = 1; j <= n; j++ {
		for i = 1; i <= n; i++ {
			a[j-1][i-1] = float32(max((i+1)-(j+1)+1, 0))
			//Label110:
		}
		x[j-1] = float32(j)
		y[j-1] = 0
		//Label120:
	}
	for j = 1; j <= n; j++ {
		yy[j-1] = float32(j*((j+1)*j))/2 - float32(((j+1)*j*(j-1)))/3
		//Label130:
	}
	//    YY holds the exact result. On exit from SMVCH YT holds
	//    the result computed by SMVCH.
	trans = 'N'
	smvch(trans, n, n, 1, &a, nmax, &x, 1, 0, &y, 1, &yt, &g, &yy, eps, err, fatal, common.infoc.nout, true)
	same = lse(&yy, &yt, n)
	if !same || err != 0 {
		writeString(common.infoc.nout, " ERROR IN SMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n SMVCH WAS CALLED WITH TRANS = %c AND RETURNED SAME =  %t AND ERR = %12.3f.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******\n", trans, same, err)
		panic("")
	}
	trans = 'T'
	smvch(trans, n, n, 1, &a, nmax, &x, -1, 0, &y, -1, &yt, &g, &yy, eps, err, fatal, common.infoc.nout, true)
	same = lse(&yy, &yt, n)
	if !same || err != 0 {
		writeString(common.infoc.nout, " ERROR IN SMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n SMVCH WAS CALLED WITH TRANS = %c AND RETURNED SAME =  %t AND ERR = %12.3f.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******\n", trans, same, err)
		panic("")
	}
	//
	//    Test each subroutine in turn.
	//
	for isnum = 1; isnum <= nsubs; isnum++ {
		writeString(common.infoc.nout, "\n")
		if !ltest[isnum-1] {
			//          Subprogram is not to be tested.
			writeString(common.infoc.nout, " %6s WAS NOT TESTED\n", snames[isnum-1])
		} else {
			common.srnamc.srnamt = snames[isnum-1]
			//          Test error exits.
			if tsterr {
				// schke(isnum, snames[isnum-1], common.infoc.nout)
				writeString(common.infoc.nout, "\n")
			}
			//          Test computations.
			common.infoc.infot = 0
			common.infoc.ok = true
			fatal = false
			switch isnum {
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
				goto Label180
			case 14:
				goto Label180
			case 15:
				goto Label190
			case 16:
				goto Label190
			}
			//          Test SGEMV, 01, and SGBMV, 02.
		Label140:
			;
			schk1(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, &idim, nkb, &kb, nalf, &alf, nbet, &bet, ninc, &inc, nmax, incmax, &a, &aa, &as, &x, &xx, &xs, &y, &yy, &ys, &yt, &g)
			goto Label200
			//          Test SSYMV, 03, SSBMV, 04, and SSPMV, 05.
		Label150:
			;
			// schk2(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, idim, nkb, kb, nalf, alf, nbet, bet, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g)
			goto Label200
			//          Test STRMV, 06, STBMV, 07, STPMV, 08,
			//          STRSV, 09, STBSV, 10, and STPSV, 11.
		Label160:
			;
			// schk3(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, idim, nkb, kb, ninc, inc, nmax, incmax, a, aa, as, y, yy, ys, yt, g, z)
			goto Label200
			//          Test SGER, 12.
		Label170:
			;
			// schk4(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			goto Label200
			//          Test SSYR, 13, and SSPR, 14.
		Label180:
			;
			// schk5(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			goto Label200
			//          Test SSYR2, 15, and SSPR2, 16.
		Label190:
			;
			// schk6(snames[isnum-1], eps, thresh, common.infoc.nout, ntra, trace, rewi, fatal, nidim, idim, nalf, alf, ninc, inc, nmax, incmax, a, aa, as, x, xx, xs, y, yy, ys, yt, g, z)
			//
		Label200:
			;
			if fatal && sfatal {
				goto Label220
			}
		}
		//Label210:
	}
	writeString(common.infoc.nout, "\n END OF TESTS\n")
	goto Label240
	//
Label220:
	;
	writeString(common.infoc.nout, "\n ******* FATAL ERROR - TESTS ABANDONED *******\n")
	goto Label240
	//
Label230:
	;
	writeString(common.infoc.nout, " AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******\n")
	//
Label240:
	;
	if trace {
		_close(ntra)
	}
	_close(common.infoc.nout)
	//
	//
	//    End of SBLAT2.
	//
}

// VERIFIED
func schk1(sname string, eps float32, thresh float32, nout int, ntra int, trace bool, rewi bool, fatal bool, nidim int, idim *[]int, nkb int, kb *[]int, nalf int, alf *[]float32, nbet int, bet *[]float32, ninc int, inc *[]int, nmax int, incmax int, a *[][]float32, aa *[]float32, as *[]float32, x *[]float32, xx *[]float32, xs *[]float32, y *[]float32, yy *[]float32, ys *[]float32, yt *[]float32, g *[]float32) {
	var alpha, als, beta, bls, err, errmax, transl float32
	var i, ia, ib, ic, iku, im, in, incx, incxs, incy, incys, ix, iy, kl, kls, ku, kus, laa, lda, ldas, lx, ly, m, ml, ms, n, nargs, nc, nd, nk, nl, ns int
	var banded, full, null, reset, same, tran bool
	var trans, transs byte
	ich := []byte{'N', 'T', 'C'}
	isame := make([]bool, 13)
	//
	// Tests SGEMV and SGBMV.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	full = sname[2] == 'E'
	banded = sname[2] == 'B'
	//    Define the number of arguments.
	if full {
		nargs = 11
	} else if banded {
		nargs = 13
	}
	//
	common.combla.n = 0
	nc = 0
	reset = true
	errmax = 0
	//
	for in = 1; in <= nidim; in++ {
		n = (*idim)[in-1]
		nd = n/2 + 1
		//
		for im = 1; im <= 2; im++ {
			if im == 1 {
				m = max(n-nd, 0)
			}
			if im == 2 {
				m = min(n+nd, nmax)
			}
			//
			if banded {
				nk = nkb
			} else {
				nk = 1
			}
			for iku = 1; iku <= nk; iku++ {
				if banded {
					ku = (*kb)[iku-1]
					kl = max(ku-1, 0)
				} else {
					ku = n - 1
					kl = m - 1
				}
				//             Set LDA to 1 more than minimum value if room.
				if banded {
					lda = kl + ku + 1
				} else {
					lda = m
				}
				if lda < nmax {
					lda++
				}
				//             Skip tests if not enough room.
				if lda > nmax {
					goto Label100
				}
				laa = lda * n
				null = n <= 0 || m <= 0
				//
				//             Generate the matrix A.
				//
				transl = 0
				smake(sname[1:3], ' ', ' ', m, n, a, nmax, aa, lda, kl, ku, reset, transl)
				//
				for ic = 1; ic <= 3; ic++ {
					trans = ich[ic-1]
					tran = trans == 'T' || trans == 'C'
					//
					if tran {
						ml = n
						nl = m
					} else {
						ml = m
						nl = n
					}
					//
					for ix = 1; ix <= ninc; ix++ {
						incx = (*inc)[ix-1]
						lx = absint(incx) * nl
						//
						//                   Generate the vector X.
						//
						transl = 0.5
						_x := make([][]float32, 1)
						for _j := range _x {
							_x[_j] = make([]float32, nl)
							for _i := 0; _i < nl; _i++ {
								_x[_j][_i] = (*x)[(nl-1)*_j+_i]
							}
						}
						smake("GE", ' ', ' ', 1, nl, &_x, 1, xx, absint(incx), 0, nl-1, reset, transl)
						for _i := range _x[0] {
							(*x)[_i] = _x[0][_i]
						}
						if nl > 1 {
							(*x)[nl/2-1] = 0
							(*xx)[1+absint(incx)*(nl/2-1)-1] = 0
						}
						//
						for iy = 1; iy <= ninc; iy++ {
							incy = (*inc)[iy-1]
							ly = absint(incy) * ml
							//
							for ia = 1; ia <= nalf; ia++ {
								alpha = (*alf)[ia-1]
								//
								for ib = 1; ib <= nbet; ib++ {
									beta = (*bet)[ib-1]
									//
									//                            Generate the vector Y.
									//
									transl = 0
									_y := make([][]float32, 1)
									for _j := range _y {
										_y[_j] = make([]float32, ml)
										for _i := 0; _i < ml; _i++ {
											_y[_j][_i] = (*y)[(ml-1)*_j+_i]
										}
									}
									smake("GE", ' ', ' ', 1, ml, &_y, 1, yy, absint(incy), 0, ml-1, reset, transl)
									for _i := range _y[0] {
										(*y)[_i] = _y[0][_i]
									}
									//
									nc++
									//
									//                            Save every datum before calling the
									//                            subroutine.
									//
									transs = trans
									ms = m
									ns = n
									kls = kl
									kus = ku
									als = alpha
									for i = 1; i <= laa; i++ {
										(*as)[i-1] = (*aa)[i-1]
										//Label10:
									}
									ldas = lda
									for i = 1; i <= lx; i++ {
										(*xs)[i-1] = (*xx)[i-1]
										//Label20:
									}
									incxs = incx
									bls = beta
									for i = 1; i <= ly; i++ {
										(*ys)[i-1] = (*yy)[i-1]
										//Label30:
									}
									incys = incy
									//
									//                            Call the subroutine.
									//
									common.combla.n++
									if full {
										if trace {
											writeString(ntra, " %6d: %6s('%c', %3d, %3d, %4.1f, a, %3d, x, %2d, %4.1f, y, %2d)\n", nc, sname, trans, m, n, alpha, lda, incx, beta, incy)
										}
										if rewi {
											rewind(ntra)
										}
										_aa := make([][]float32, m)
										for _j := range _aa {
											_aa[_j] = make([]float32, n)
											for _i := 0; _i < n; _i++ {
												_aa[_j][_i] = (*aa)[(n-1)*_j+_i]
											}
										}
										Sgemv(trans, m, n, alpha, &_aa, lda, xx, incx, beta, yy, incy)
									} else if banded {
										if trace {
											writeString(ntra, " %6d: %6s('%c', %3d, %3d, %3d, %3d, %4.1f, a, %3d, x, %2d, %4.1f, y, %2d)\n", nc, sname, trans, m, n, kl, ku, alpha, lda, incx, beta, incy)
										}
										if rewi {
											rewind(ntra)
										}
										_aa := make([][]float32, m)
										for _j := range _aa {
											_aa[_j] = make([]float32, n)
											for _i := 0; _i < n; _i++ {
												_aa[_j][_i] = (*aa)[(n-1)*_j+_i]
											}
										}
										Sgbmv(trans, m, n, kl, ku, alpha, &_aa, lda, xx, incx, beta, yy, incy)
									}
									//
									//                            Check if error-exit was taken incorrectly.
									//
									if !common.infoc.ok {
										writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n")
										fatal = true
										goto Label130
									}
									//
									//                            See what data changed inside subroutines.
									//
									isame[0] = trans == transs
									isame[1] = ms == m
									isame[2] = ns == n
									if full {
										isame[3] = als == alpha
										isame[4] = lse(as, aa, laa)
										isame[5] = ldas == lda
										isame[6] = lse(xs, xx, lx)
										isame[7] = incxs == incx
										isame[8] = bls == beta
										if null {
											isame[9] = lse(ys, yy, ly)
										} else {
											// Generate m x n matrix from vector
											_ys := make([][]float32, 1)
											_yy := make([][]float32, 1)
											for _j := range _ys {
												_ys[_j] = make([]float32, ml*incy)
												_yy[_j] = make([]float32, ml*incy)
												for _i := 0; _i < ml; _i++ {
													_ys[_j][_i] = (*ys)[(ml*incy-1)*_j+_i]
													_yy[_j][_i] = (*yy)[(ml*incy-1)*_j+_i]
												}
											}
											isame[9] = lseres("GE", ' ', 1, ml, &_ys, &_yy, absint(incy))
										}
										isame[10] = incys == incy
									} else if banded {
										isame[3] = kls == kl
										isame[4] = kus == ku
										isame[5] = als == alpha
										isame[6] = lse(as, aa, laa)
										isame[7] = ldas == lda
										isame[8] = lse(xs, xx, lx)
										isame[9] = incxs == incx
										isame[10] = bls == beta
										if null {
											isame[11] = lse(ys, yy, ly)
										} else {
											_ys := make([][]float32, 1)
											_yy := make([][]float32, 1)
											for _j := range _ys {
												_ys[_j] = make([]float32, ml)
												_yy[_j] = make([]float32, ml)
												for _i := 0; _i < ml; _i++ {
													_ys[_j][_i] = (*ys)[(ml-1)*_j+_i]
													_yy[_j][_i] = (*yy)[(ml-1)*_j+_i]
												}
											}
											isame[11] = lseres("GE", ' ', 1, ml, &_ys, &_yy, absint(incy))
										}
										isame[12] = incys == incy
									}
									//
									//                            If data was incorrectly changed, report
									//                            and return.
									//
									same = true
									for i = 1; i <= nargs; i++ {
										same = same && isame[i-1]
										if !isame[i-1] {
											writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n", i)
										}
										//Label40:
									}
									if !same {
										fatal = true
										goto Label130
									}
									//
									if !null {
										//
										//                               Check the result.
										//
										smvch(trans, m, n, alpha, a, nmax, x, incx, beta, y, incy, yt, g, yy, eps, err, fatal, nout, true)
										errmax = maxf32(errmax, err)
										//                               If got really bad answer, report and
										//                               return.
										if fatal {
											goto Label130
										}
									} else {
										//                               Avoid repeating tests with M.le.0 or
										//                               N.le.0.
										goto Label110
									}
									//
									//Label50:
								}
								//
								//Label60:
							}
							//
							//Label70:
						}
						//
						//Label80:
					}
					//
					//Label90:
				}
				//
			Label100:
				switch i {
				case -1:
					fmt.Println("")
				}
			}
			//
		Label110:
			switch i {
			case -1:
				fmt.Println("")
			}
		}
		//
		//Label120:
	}
	//
	//    Report result.
	//
	if errmax < thresh {
		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n", sname, nc)
	} else {
		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n", sname, nc, errmax)
	}
	goto Label140
	//
Label130:
	;
	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n", sname)
	if full {
		writeString(common.infoc.nout, " %6d: %6s('%c', %3d, %3d, %4.1f, a, %3d, x, %2d, %4.1f, y, %2d)\n", nc, sname, trans, m, n, alpha, lda, incx, beta, incy)
	} else if banded {
		writeString(common.infoc.nout, " %6d: %6s('%c', %3d, %3d, %3d, %3d, %4.1f, a, %3d, x, %2d, %4.1f, y, %2d)\n", nc, sname, trans, m, n, kl, ku, alpha, lda, incx, beta, incy)
	}
	//
Label140:
	;
	return
	//
	//
	//    End of SCHK1.
	//
}

// func schk2(SNAME []byte, EPS float32, THRESH float32, NOUT int, NTRA int, TRACE bool, REWI bool, FATAL bool, NIDIM int, IDIM *[]int, NKB int, KB *[]int, NALF int, ALF *[]float32, NBET int, BET *[]float32, NINC int, INC *[]int, NMAX int, INCMAX int, A *[][]float32, AA *[]float32, AS *[]float32, X *[]float32, XX *[]float32, XS *[]float32, Y *[]float32, YY *[]float32, YS *[]float32, YT *[]float32, G *[]float32) {
// 	ZERO := new(float32)
// 	HALF := new(float32)
// 	ALPHA := new(float32)
// 	ALS := new(float32)
// 	BETA := new(float32)
// 	BLS := new(float32)
// 	ERR := new(float32)
// 	ERRMAX := new(float32)
// 	TRANSL := new(float32)
// 	I := new(int)
// 	IA := new(int)
// 	IB := new(int)
// 	IC := new(int)
// 	IK := new(int)
// 	IN := new(int)
// 	INCX := new(int)
// 	INCXS := new(int)
// 	INCY := new(int)
// 	INCYS := new(int)
// 	IX := new(int)
// 	IY := new(int)
// 	K := new(int)
// 	KS := new(int)
// 	LAA := new(int)
// 	LDA := new(int)
// 	LDAS := new(int)
// 	LX := new(int)
// 	LY := new(int)
// 	N := new(int)
// 	NARGS := new(int)
// 	NC := new(int)
// 	NK := new(int)
// 	NS := new(int)
// 	BANDED := new(bool)
// 	FULL := new(bool)
// 	NULL := new(bool)
// 	PACKED := new(bool)
// 	RESET := new(bool)
// 	SAME := new(bool)
// 	UPLO := new(byte)
// 	UPLOS := new(byte)
// 	ICH := func() *[]byte {
// 		arr := make([]byte, 2)
// 		return &arr
// 	}()
// 	ISAME := func() *[]bool {
// 		arr := make([]bool, 13)
// 		return &arr
// 	}()
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests SSYMV, SSBMV and SSPMV.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Parameters ..
// 	0 = 0.0
// 	0.5 = 0.5
// 	//    .. Scalar Arguments ..
// 	//    .. Array Arguments ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Functions ..
// 	//    .. External Subroutines ..
// 	//    .. Intrinsic Functions ..
// 	//    .. Scalars in Common ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Data statements ..
// 	ich[0], ich[1] = 'U', 'L'
// 	//    .. Executable Statements ..
// 	full = sname[2] == 'Y'
// 	banded = sname[2] == 'B'
// 	(*PACKED) = sname[2] == 'P'
// 	//    Define the number of arguments.
// 	if full {
// 		nargs = 10
// 	} else if banded {
// 		nargs = 11
// 	} else if *PACKED {
// 		nargs = 9
// 	}
// 	//
// 	nc = 0
// 	reset = true
// 	errmax = 0
// 	//
// 	for in = 1; in <= nidim; in++ {
// 		n = (*idim)[in-1]
// 		//
// 		if banded {
// 			nk = nkb
// 		} else {
// 			nk = 1
// 		}
// 		for (*IK) = 1; (*IK) <= nk; (*IK)++ {
// 			if banded {
// 				(*K) = (*kb)[(*IK)-1]
// 			} else {
// 				(*K) = n - 1
// 			}
// 			//          Set LDA to 1 more than minimum value if room.
// 			if banded {
// 				lda = (*K) + 1
// 			} else {
// 				lda = n
// 			}
// 			if lda < nmax {
// 				lda = lda + 1
// 			}
// 			//          Skip tests if not enough room.
// 			if lda > nmax {
// 				goto Label100
// 			}
// 			if *PACKED {
// 				laa = (n * (n + 1)) / 2
// 			} else {
// 				laa = lda * n
// 			}
// 			null = n <= 0
// 			//
// 			for ic = 1; ic <= 2; ic++ {
// 				(*UPLO) = ich[ic-1]
// 				//
// 				//             Generate the matrix A.
// 				//
// 				transl = 0
// 				smake(sname[1], UPLO, ' '), n, n, a, nmax, aa, lda, K, K, reset, transl)
// 				//
// 				for ix = 1; ix <= ninc; ix++ {
// 					incx = (*inc)[ix-1]
// 					lx = absint(incx) * n
// 					//
// 					//                Generate the vector X.
// 					//
// 					transl = 0.5
// 					smake("GE"), ' '), ' '), 1, n, x, 1, xx, absint(incx), 0, n-1, reset, transl)
// 					if n > 1 {
// 						(*x)[n/1] = 0
// 						(*xx)[1+absint(incx)*(n/2-1)-1] = 0
// 					}
// 					//
// 					for iy = 1; iy <= ninc; iy++ {
// 						incy = (*inc)[iy-1]
// 						ly = absint(incy) * n
// 						//
// 						for ia = 1; ia <= nalf; ia++ {
// 							alpha = (*alf)[ia-1]
// 							//
// 							for ib = 1; ib <= nbet; ib++ {
// 								beta = (*bet)[ib-1]
// 								//
// 								//                         Generate the vector Y.
// 								//
// 								transl = 0
// 								smake("GE"), ' '), ' '), 1, n, y, 1, yy, absint(incy), 0, n-1, reset, transl)
// 								//
// 								nc = nc + 1
// 								//
// 								//                         Save every datum before calling the
// 								//                         subroutine.
// 								//
// 								(*UPLOS) = (*UPLO)
// 								ns = n
// 								(*KS) = (*K)
// 								als = alpha
// 								for i = 1; i <= laa; i++ {
// 									(*as)[i-1] = (*aa)[i-1]
// 								//Label10:
// 								}
// 								ldas = lda
// 								for i = 1; i <= lx; i++ {
// 									(*xs)[i-1] = (*xx)[i-1]
// 								//Label20:
// 								}
// 								incxs = incx
// 								bls = beta
// 								for i = 1; i <= ly; i++ {
// 									(*ys)[i-1] = (*yy)[i-1]
// 								//Label30:
// 								}
// 								incys = incy
// 								//
// 								//                         Call the subroutine.
// 								//
// 								if full {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, a,%3d, x,%2d,%4.1f, y,%2d)             .\n"), nc, sname, (*UPLO), n, alpha, lda, incx, beta, incy)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									SSYMV(UPLO, n, alpha, aa, lda, xx, incx, beta, yy, incy)
// 								} else if banded {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, a,%3d, x,%2d,%4.1f, y,%2d)         .\n"), nc, sname, (*UPLO), n, (*K), alpha, lda, incx, beta, incy)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									SSBMV(UPLO, n, K, alpha, aa, lda, xx, incx, beta, yy, incy)
// 								} else if *PACKED {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, AP, x,%2d,%4.1f, y,%2d)                .\n"), nc, sname, (*UPLO), n, alpha, incx, beta, incy)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									SSPMV(UPLO, n, alpha, aa, xx, incx, beta, yy, incy)
// 								}
// 								//
// 								//                         Check if error-exit was taken incorrectly.
// 								//
// 								if !common.infoc.ok {
// 									writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n"))
// 									fatal = true
// 									goto Label120
// 								}
// 								//
// 								//                         See what data changed inside subroutines.
// 								//
// 								isame[0] = (*UPLO) == (*UPLOS)
// 								isame[1] = ns == n
// 								if full {
// 									isame[2] = als == alpha
// 									isame[3] = lse(as, aa, laa))
// 									isame[4] = ldas == lda
// 									isame[5] = lse(xs, xx, lx))
// 									isame[6] = incxs == incx
// 									isame[7] = bls == beta
// 									if null {
// 										isame[8] = lse(ys, yy, ly))
// 									} else {
// 										isame[8] = (*lseres("GE"), ' '), 1, n, ys, yy, absint(incy)))
// 									}
// 									isame[9] = incys == incy
// 								} else if banded {
// 									isame[2] = (*KS) == (*K)
// 									isame[3] = als == alpha
// 									isame[4] = lse(as, aa, laa))
// 									isame[5] = ldas == lda
// 									isame[6] = lse(xs, xx, lx))
// 									isame[7] = incxs == incx
// 									isame[8] = bls == beta
// 									if null {
// 										isame[9] = lse(ys, yy, ly))
// 									} else {
// 										isame[9] = (*lseres("GE"), ' '), 1, n, ys, yy, absint(incy)))
// 									}
// 									isame[10] = incys == incy
// 								} else if *PACKED {
// 									isame[2] = als == alpha
// 									isame[3] = lse(as, aa, laa))
// 									isame[4] = lse(xs, xx, lx))
// 									isame[5] = incxs == incx
// 									isame[6] = bls == beta
// 									if null {
// 										isame[7] = lse(ys, yy, ly))
// 									} else {
// 										isame[7] = (*lseres("GE"), ' '), 1, n, ys, yy, absint(incy)))
// 									}
// 									isame[8] = incys == incy
// 								}
// 								//
// 								//                         If data was incorrectly changed, report and
// 								//                         return.
// 								//
// 								same = true
// 								for i = 1; i <= nargs; i++ {
// 									same = same && isame[i-1]
// 									if !isame[i-1] {
// 										writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n"), i)
// 									}
// 								//Label40:
// 								}
// 								if !same {
// 									fatal = true
// 									goto Label120
// 								}
// 								//
// 								if !null {
// 									//
// 									//                            Check the result.
// 									//
// 									smvch('N'), n, n, alpha, a, nmax, x, incx, beta, y, incy, yt, g, yy, eps, err, fatal, nout, true)
// 									errmax = (max(errmax, err))
// 									//                            If got really bad answer, report and
// 									//                            return.
// 									if fatal {
// 										goto Label120
// 									}
// 								} else {
// 									//                            Avoid repeating tests with N.le.0
// 									goto Label110
// 								}
// 								//
// 							//Label50:
// 							}
// 							//
// 						//Label60:
// 						}
// 						//
// 					//Label70:
// 					}
// 					//
// 				//Label80:
// 				}
// 				//
// 			//Label90:
// 			}
// 			//
// 		Label100:
// 		}
// 		//
// 	Label110:
// 	}
// 	//
// 	//    Report result.
// 	//
// 	if errmax < thresh {
// 		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n"), sname, nc)
// 	} else {
// 		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n"), sname, nc, errmax)
// 	}
// 	goto Label130
// 	//
// Label120:
// 	;
// 	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n"), sname)
// 	if full {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, a,%3d, x,%2d,%4.1f, y,%2d)             .\n"), nc, sname, (*UPLO), n, alpha, lda, incx, beta, incy)
// 	} else if banded {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, a,%3d, x,%2d,%4.1f, y,%2d)         .\n"), nc, sname, (*UPLO), n, (*K), alpha, lda, incx, beta, incy)
// 	} else if *PACKED {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, AP, x,%2d,%4.1f, y,%2d)                .\n"), nc, sname, (*UPLO), n, alpha, incx, beta, incy)
// 	}
// 	//
// Label130:
// 	;
// 	return
// 	//
// 	//
// 	//    End of SCHK2.
// 	//
// }

// func schk3(SNAME []byte, EPS float32, THRESH float32, NOUT int, NTRA int, TRACE bool, REWI bool, FATAL bool, NIDIM int, IDIM *[]int, NKB int, KB *[]int, NINC int, INC *[]int, NMAX int, INCMAX int, A *[][]float32, AA *[]float32, AS *[]float32, X *[]float32, XX *[]float32, XS *[]float32, XT *[]float32, G *[]float32, Z *[]float32) {
// 	ZERO := new(float32)
// 	HALF := new(float32)
// 	ONE := new(float32)
// 	ERR := new(float32)
// 	ERRMAX := new(float32)
// 	TRANSL := new(float32)
// 	I := new(int)
// 	ICD := new(int)
// 	ICT := new(int)
// 	ICU := new(int)
// 	IK := new(int)
// 	IN := new(int)
// 	INCX := new(int)
// 	INCXS := new(int)
// 	IX := new(int)
// 	K := new(int)
// 	KS := new(int)
// 	LAA := new(int)
// 	LDA := new(int)
// 	LDAS := new(int)
// 	LX := new(int)
// 	N := new(int)
// 	NARGS := new(int)
// 	NC := new(int)
// 	NK := new(int)
// 	NS := new(int)
// 	BANDED := new(bool)
// 	FULL := new(bool)
// 	NULL := new(bool)
// 	PACKED := new(bool)
// 	RESET := new(bool)
// 	SAME := new(bool)
// 	DIAG := new(byte)
// 	DIAGS := new(byte)
// 	TRANS := new(byte)
// 	TRANSS := new(byte)
// 	UPLO := new(byte)
// 	UPLOS := new(byte)
// 	ICHD := func() *[]byte {
// 		arr := make([]byte, 2)
// 		return &arr
// 	}()
// 	ICHU := func() *[]byte {
// 		arr := make([]byte, 2)
// 		return &arr
// 	}()
// 	ICHT := func() *[]byte {
// 		arr := make([]byte, 3)
// 		return &arr
// 	}()
// 	ISAME := func() *[]bool {
// 		arr := make([]bool, 13)
// 		return &arr
// 	}()
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests STRMV, STBMV, STPMV, STRSV, STBSV and STPSV.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Parameters ..
// 	0 = 0.0
// 	0.5 = 0.5
// 	1 = 1.0
// 	//    .. Scalar Arguments ..
// 	//    .. Array Arguments ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Functions ..
// 	//    .. External Subroutines ..
// 	//    .. Intrinsic Functions ..
// 	//    .. Scalars in Common ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Data statements ..
// 	(*ICHU)[0], (*ICHU)[1], (*ICHT)[0], (*ICHT)[1], (*ICHT)[2], (*ICHD)[0], (*ICHD)[1] = 'U', 'L', 'N', 'T', 'C', 'U', 'N'
// 	//    .. Executable Statements ..
// 	full = sname[2] == 'R'
// 	banded = sname[2] == 'B'
// 	(*PACKED) = sname[2] == 'P'
// 	//    Define the number of arguments.
// 	if full {
// 		nargs = 8
// 	} else if banded {
// 		nargs = 9
// 	} else if *PACKED {
// 		nargs = 7
// 	}
// 	//
// 	nc = 0
// 	reset = true
// 	errmax = 0
// 	//    Set up zero vector for SMVCH.
// 	for i = 1; i <= nmax; i++ {
// 		(*(Z))[i-1] = 0
// 	//Label10:
// 	}
// 	//
// 	for in = 1; in <= nidim; in++ {
// 		n = (*idim)[in-1]
// 		//
// 		if banded {
// 			nk = nkb
// 		} else {
// 			nk = 1
// 		}
// 		for (*IK) = 1; (*IK) <= nk; (*IK)++ {
// 			if banded {
// 				(*K) = (*kb)[(*IK)-1]
// 			} else {
// 				(*K) = n - 1
// 			}
// 			//          Set LDA to 1 more than minimum value if room.
// 			if banded {
// 				lda = (*K) + 1
// 			} else {
// 				lda = n
// 			}
// 			if lda < nmax {
// 				lda = lda + 1
// 			}
// 			//          Skip tests if not enough room.
// 			if lda > nmax {
// 				goto Label100
// 			}
// 			if *PACKED {
// 				laa = (n * (n + 1)) / 2
// 			} else {
// 				laa = lda * n
// 			}
// 			null = n <= 0
// 			//
// 			for (*ICU) = 1; (*ICU) <= 2; (*ICU)++ {
// 				(*UPLO) = (*ICHU)[(*ICU)-1]
// 				//
// 				for (*ICT) = 1; (*ICT) <= 3; (*ICT)++ {
// 					trans = (*ICHT)[(*ICT)-1]
// 					//
// 					for (*ICD) = 1; (*ICD) <= 2; (*ICD)++ {
// 						(*DIAG) = (*ICHD)[(*ICD)-1]
// 						//
// 						//                   Generate the matrix A.
// 						//
// 						transl = 0
// 						smake(sname[1], UPLO, DIAG, n, n, a, nmax, aa, lda, K, K, reset, transl)
// 						//
// 						for ix = 1; ix <= ninc; ix++ {
// 							incx = (*inc)[ix-1]
// 							lx = absint(incx) * n
// 							//
// 							//                      Generate the vector X.
// 							//
// 							transl = 0.5
// 							smake("GE"), ' '), ' '), 1, n, x, 1, xx, absint(incx), 0, n-1, reset, transl)
// 							if n > 1 {
// 								(*x)[n/1] = 0
// 								(*xx)[1+absint(incx)*(n/2-1)-1] = 0
// 							}
// 							//
// 							nc = nc + 1
// 							//
// 							//                      Save every datum before calling the subroutine.
// 							//
// 							(*UPLOS) = (*UPLO)
// 							transs = trans
// 							(*DIAGS) = (*DIAG)
// 							ns = n
// 							(*KS) = (*K)
// 							for i = 1; i <= laa; i++ {
// 								(*as)[i-1] = (*aa)[i-1]
// 							//Label20:
// 							}
// 							ldas = lda
// 							for i = 1; i <= lx; i++ {
// 								(*xs)[i-1] = (*xx)[i-1]
// 							//Label30:
// 							}
// 							incxs = incx
// 							//
// 							//                      Call the subroutine.
// 							//
// 							if sname[3] == "MV") {
// 								if full {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                     .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, lda, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STRMV(UPLO, TRANS, DIAG, n, aa, lda, xx, INCX)
// 								} else if banded {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                 .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, (*K), lda, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STBMV(UPLO, TRANS, DIAG, n, K, aa, lda, xx, INCX)
// 								} else if *PACKED {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, AP, x,%2d)                        .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STPMV(UPLO, TRANS, DIAG, n, aa, xx, INCX)
// 								}
// 							} else if sname[3] == "SV") {
// 								if full {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                     .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, lda, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STRSV(UPLO, TRANS, DIAG, n, aa, lda, xx, INCX)
// 								} else if banded {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                 .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, (*K), lda, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STBSV(UPLO, TRANS, DIAG, n, K, aa, lda, xx, INCX)
// 								} else if *PACKED {
// 									if trace {
// 										writeString(ntra, " %6d: %6s('%c',%3d, AP, x,%2d)                        .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, incx)
// 									}
// 									if rewi {
// 										rewind((ntra))
// 									}
// 									STPSV(UPLO, TRANS, DIAG, n, aa, xx, INCX)
// 								}
// 							}
// 							//
// 							//                      Check if error-exit was taken incorrectly.
// 							//
// 							if !common.infoc.ok {
// 								writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n"))
// 								fatal = true
// 								goto Label120
// 							}
// 							//
// 							//                      See what data changed inside subroutines.
// 							//
// 							isame[0] = (*UPLO) == (*UPLOS)
// 							isame[1] = trans == transs
// 							isame[2] = (*DIAG) == (*DIAGS)
// 							isame[3] = ns == n
// 							if full {
// 								isame[4] = lse(as, aa, laa))
// 								isame[5] = ldas == lda
// 								if null {
// 									isame[6] = lse(xs, xx, lx))
// 								} else {
// 									isame[6] = (*lseres("GE"), ' '), 1, n, xs, xx, absint(incx)))
// 								}
// 								isame[7] = incxs == incx
// 							} else if banded {
// 								isame[4] = (*KS) == (*K)
// 								isame[5] = lse(as, aa, laa))
// 								isame[6] = ldas == lda
// 								if null {
// 									isame[7] = lse(xs, xx, lx))
// 								} else {
// 									isame[7] = (*lseres("GE"), ' '), 1, n, xs, xx, absint(incx)))
// 								}
// 								isame[8] = incxs == incx
// 							} else if *PACKED {
// 								isame[4] = lse(as, aa, laa))
// 								if null {
// 									isame[5] = lse(xs, xx, lx))
// 								} else {
// 									isame[5] = (*lseres("GE"), ' '), 1, n, xs, xx, absint(incx)))
// 								}
// 								isame[6] = incxs == incx
// 							}
// 							//
// 							//                      If data was incorrectly changed, report and
// 							//                      return.
// 							//
// 							same = true
// 							for i = 1; i <= nargs; i++ {
// 								same = same && isame[i-1]
// 								if !isame[i-1] {
// 									writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n"), i)
// 								}
// 							//Label40:
// 							}
// 							if !same {
// 								fatal = true
// 								goto Label120
// 							}
// 							//
// 							if !null {
// 								if sname[3] == "MV") {
// 									//
// 									//                            Check the result.
// 									//
// 									smvch(trans, n, n, 1, a, nmax, x, incx, 0, (Z), incx, (XT), g, xx, eps, err, fatal, nout, true)
// 								} else if sname[3] == "SV") {
// 									//
// 									//                            Compute approximation to original vector.
// 									//
// 									for i = 1; i <= n; i++ {
// 										(*(Z))[i-1] = (*xx)[1+(i-1)absint(incx)-1]
// 										(*xx)[1+(i-1)absint(incx)-1] = (*x)[i-1]
// 									//Label50:
// 									}
// 									smvch(trans, n, n, 1, a, nmax, (Z), incx, 0, x, incx, (XT), g, xx, eps, err, fatal, nout, false)
// 								}
// 								errmax = (max(errmax, err))
// 								//                         If got really bad answer, report and return.
// 								if fatal {
// 									goto Label120
// 								}
// 							} else {
// 								//                         Avoid repeating tests with N.le.0.
// 								goto Label110
// 							}
// 							//
// 						//Label60:
// 						}
// 						//
// 					//Label70:
// 					}
// 					//
// 				//Label80:
// 				}
// 				//
// 			//Label90:
// 			}
// 			//
// 		Label100:
// 		}
// 		//
// 	Label110:
// 	}
// 	//
// 	//    Report result.
// 	//
// 	if errmax < thresh {
// 		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n"), sname, nc)
// 	} else {
// 		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n"), sname, nc, errmax)
// 	}
// 	goto Label130
// 	//
// Label120:
// 	;
// 	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n"), sname)
// 	if full {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                     .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, lda, incx)
// 	} else if banded {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d, a,%3d, x,%2d)                 .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, (*K), lda, incx)
// 	} else if *PACKED {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d, AP, x,%2d)                        .\n"), nc, sname, (*UPLO), trans, (*DIAG), n, incx)
// 	}
// 	//
// Label130:
// 	;
// 	return
// 	//
// 	//
// 	//    End of SCHK3.
// 	//
// }

// func schk4(SNAME []byte, EPS float32, THRESH float32, NOUT int, NTRA int, TRACE bool, REWI bool, FATAL bool, NIDIM int, IDIM *[]int, NALF int, ALF *[]float32, NINC int, INC *[]int, NMAX int, INCMAX int, A *[][]float32, AA *[]float32, AS *[]float32, X *[]float32, XX *[]float32, XS *[]float32, Y *[]float32, YY *[]float32, YS *[]float32, YT *[]float32, G *[]float32, Z *[]float32) {
// 	ZERO := new(float32)
// 	HALF := new(float32)
// 	ONE := new(float32)
// 	ALPHA := new(float32)
// 	ALS := new(float32)
// 	ERR := new(float32)
// 	ERRMAX := new(float32)
// 	TRANSL := new(float32)
// 	I := new(int)
// 	IA := new(int)
// 	IM := new(int)
// 	IN := new(int)
// 	INCX := new(int)
// 	INCXS := new(int)
// 	INCY := new(int)
// 	INCYS := new(int)
// 	IX := new(int)
// 	IY := new(int)
// 	J := new(int)
// 	LAA := new(int)
// 	LDA := new(int)
// 	LDAS := new(int)
// 	LX := new(int)
// 	LY := new(int)
// 	M := new(int)
// 	MS := new(int)
// 	N := new(int)
// 	NARGS := new(int)
// 	NC := new(int)
// 	ND := new(int)
// 	NS := new(int)
// 	NULL := new(bool)
// 	RESET := new(bool)
// 	SAME := new(bool)
// 	W := func() *[]float32 {
// 		arr := make([]float32, 1)
// 		return &arr
// 	}()
// 	ISAME := func() *[]bool {
// 		arr := make([]bool, 13)
// 		return &arr
// 	}()
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests SGER.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Parameters ..
// 	0 = 0.0
// 	0.5 = 0.5
// 	1 = 1.0
// 	//    .. Scalar Arguments ..
// 	//    .. Array Arguments ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Functions ..
// 	//    .. External Subroutines ..
// 	//    .. Intrinsic Functions ..
// 	//    .. Scalars in Common ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Executable Statements ..
// 	//    Define the number of arguments.
// 	nargs = 9
// 	//
// 	nc = 0
// 	reset = true
// 	errmax = 0
// 	//
// 	for in = 1; in <= nidim; in++ {
// 		n = (*idim)[in-1]
// 		nd = n/2 + 1
// 		//
// 		for im = 1; im <= 2; im++ {
// 			if im == 1 {
// 				m = (max(n-nd, 0))
// 			}
// 			if im == 2 {
// 				m = (min(n+nd, nmax))
// 			}
// 			//
// 			//          Set LDA to 1 more than minimum value if room.
// 			lda = m
// 			if lda < nmax {
// 				lda = lda + 1
// 			}
// 			//          Skip tests if not enough room.
// 			if lda > nmax {
// 				goto Label110
// 			}
// 			laa = lda * n
// 			null = n <= 0 || m <= 0
// 			//
// 			for ix = 1; ix <= ninc; ix++ {
// 				incx = (*inc)[ix-1]
// 				lx = absint(incx) * m
// 				//
// 				//             Generate the vector X.
// 				//
// 				transl = 0.5
// 				smake("GE"), ' '), ' '), 1, m, x, 1, xx, absint(incx), 0, m-1, reset, transl)
// 				if m > 1 {
// 					(*x)[m/1] = 0
// 					(*xx)[1+absint(incx)*(m/2-1)-1] = 0
// 				}
// 				//
// 				for iy = 1; iy <= ninc; iy++ {
// 					incy = (*inc)[iy-1]
// 					ly = absint(incy) * n
// 					//
// 					//                Generate the vector Y.
// 					//
// 					transl = 0
// 					smake("GE"), ' '), ' '), 1, n, y, 1, yy, absint(incy), 0, n-1, reset, transl)
// 					if n > 1 {
// 						(*y)[n/1] = 0
// 						(*yy)[1+absint(incy)*(n/2-1)-1] = 0
// 					}
// 					//
// 					for ia = 1; ia <= nalf; ia++ {
// 						alpha = (*alf)[ia-1]
// 						//
// 						//                   Generate the matrix A.
// 						//
// 						transl = 0
// 						smake(sname[1], ' '), ' '), m, n, a, nmax, aa, lda, m-1, n-1, reset, transl)
// 						//
// 						nc = nc + 1
// 						//
// 						//                   Save every datum before calling the subroutine.
// 						//
// 						ms = m
// 						ns = n
// 						als = alpha
// 						for i = 1; i <= laa; i++ {
// 							(*as)[i-1] = (*aa)[i-1]
// 						//Label10:
// 						}
// 						ldas = lda
// 						for i = 1; i <= lx; i++ {
// 							(*xs)[i-1] = (*xx)[i-1]
// 						//Label20:
// 						}
// 						incxs = incx
// 						for i = 1; i <= ly; i++ {
// 							(*ys)[i-1] = (*yy)[i-1]
// 						//Label30:
// 						}
// 						incys = incy
// 						//
// 						//                   Call the subroutine.
// 						//
// 						if trace {
// 							writeString(ntra, " %6d: %6s(%3d,%4.1f, x,%2d, y,%2d, a,%3d)                  .\n"), nc, sname, m, n, alpha, incx, incy, lda)
// 						}
// 						if rewi {
// 							rewind((ntra))
// 						}
// 						SGER(M, n, alpha, xx, incx, yy, incy, aa, LDA)
// 						//
// 						//                   Check if error-exit was taken incorrectly.
// 						//
// 						if !common.infoc.ok {
// 							writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n"))
// 							fatal = true
// 							goto Label140
// 						}
// 						//
// 						//                   See what data changed inside subroutine.
// 						//
// 						isame[0] = ms == m
// 						isame[1] = ns == n
// 						isame[2] = als == alpha
// 						isame[3] = lse(xs, xx, lx))
// 						isame[4] = incxs == incx
// 						isame[5] = lse(ys, yy, ly))
// 						isame[6] = incys == incy
// 						if null {
// 							isame[7] = lse(as, aa, laa))
// 						} else {
// 							isame[7] = (*lseres("GE"), ' '), m, n, as, aa, LDA))
// 						}
// 						isame[8] = ldas == lda
// 						//
// 						//                   If data was incorrectly changed, report and return.
// 						//
// 						same = true
// 						for i = 1; i <= nargs; i++ {
// 							same = same && isame[i-1]
// 							if !isame[i-1] {
// 								writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n"), i)
// 							}
// 						//Label40:
// 						}
// 						if !same {
// 							fatal = true
// 							goto Label140
// 						}
// 						//
// 						if !null {
// 							//
// 							//                      Check the result column by column.
// 							//
// 							if incx > 0 {
// 								for i = 1; i <= m; i++ {
// 									(*(Z))[i-1] = (*x)[i-1]
// 								//Label50:
// 								}
// 							} else {
// 								for i = 1; i <= m; i++ {
// 									(*(Z))[i-1] = (*x)[m-i+0]
// 								//Label60:
// 								}
// 							}
// 							for j = 1; j <= n; j++ {
// 								if incy > 0 {
// 									(*W)[0] = (*y)[j-1]
// 								} else {
// 									(*W)[0] = (*y)[n-j+0]
// 								}
// 								smvch('N'), m, 1, alpha, (Z), nmax, W, 1, 1, &((*a)[0][j-1]), 1, yt, g, &((*aa)[1+(j-1)*lda-1]), eps, err, fatal, nout, true)
// 								errmax = (max(errmax, err))
// 								//                         If got really bad answer, report and return.
// 								if fatal {
// 									goto Label130
// 								}
// 							//Label70:
// 							}
// 						} else {
// 							//                      Avoid repeating tests with M.le.0 or N.le.0.
// 							goto Label110
// 						}
// 						//
// 					//Label80:
// 					}
// 					//
// 				//Label90:
// 				}
// 				//
// 			//Label100:
// 			}
// 			//
// 		Label110:
// 		}
// 		//
// 	//Label120:
// 	}
// 	//
// 	//    Report result.
// 	//
// 	if errmax < thresh {
// 		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n"), sname, nc)
// 	} else {
// 		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n"), sname, nc, errmax)
// 	}
// 	goto Label150
// 	//
// Label130:
// 	;
// 	writeString(common.infoc.nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n"), j)
// 	//
// Label140:
// 	;
// 	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n"), sname)
// 	writeString(common.infoc.nout, " %6d: %6s(%3d,%4.1f, x,%2d, y,%2d, a,%3d)                  .\n"), nc, sname, m, n, alpha, incx, incy, lda)
// 	//
// Label150:
// 	;
// 	return
// 	//
// 	//
// 	//    End of SCHK4.
// 	//
// }

// func schk5(SNAME []byte, EPS float32, THRESH float32, NOUT int, NTRA int, TRACE bool, REWI bool, FATAL bool, NIDIM int, IDIM *[]int, NALF int, ALF *[]float32, NINC int, INC *[]int, NMAX int, INCMAX int, A *[][]float32, AA *[]float32, AS *[]float32, X *[]float32, XX *[]float32, XS *[]float32, Y *[]float32, YY *[]float32, YS *[]float32, YT *[]float32, G *[]float32, Z *[]float32) {
// 	ZERO := new(float32)
// 	HALF := new(float32)
// 	ONE := new(float32)
// 	ALPHA := new(float32)
// 	ALS := new(float32)
// 	ERR := new(float32)
// 	ERRMAX := new(float32)
// 	TRANSL := new(float32)
// 	I := new(int)
// 	IA := new(int)
// 	IC := new(int)
// 	IN := new(int)
// 	INCX := new(int)
// 	INCXS := new(int)
// 	IX := new(int)
// 	J := new(int)
// 	JA := new(int)
// 	JJ := new(int)
// 	LAA := new(int)
// 	LDA := new(int)
// 	LDAS := new(int)
// 	LJ := new(int)
// 	LX := new(int)
// 	N := new(int)
// 	NARGS := new(int)
// 	NC := new(int)
// 	NS := new(int)
// 	FULL := new(bool)
// 	NULL := new(bool)
// 	PACKED := new(bool)
// 	RESET := new(bool)
// 	SAME := new(bool)
// 	UPPER := new(bool)
// 	UPLO := new(byte)
// 	UPLOS := new(byte)
// 	ICH := func() *[]byte {
// 		arr := make([]byte, 2)
// 		return &arr
// 	}()
// 	W := func() *[]float32 {
// 		arr := make([]float32, 1)
// 		return &arr
// 	}()
// 	ISAME := func() *[]bool {
// 		arr := make([]bool, 13)
// 		return &arr
// 	}()
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests SSYR and SSPR.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Parameters ..
// 	0 = 0.0
// 	0.5 = 0.5
// 	1 = 1.0
// 	//    .. Scalar Arguments ..
// 	//    .. Array Arguments ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Functions ..
// 	//    .. External Subroutines ..
// 	//    .. Intrinsic Functions ..
// 	//    .. Scalars in Common ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Data statements ..
// 	ich[0], ich[1] = 'U', 'L'
// 	//    .. Executable Statements ..
// 	full = sname[2] == 'Y'
// 	(*PACKED) = sname[2] == 'P'
// 	//    Define the number of arguments.
// 	if full {
// 		nargs = 7
// 	} else if *PACKED {
// 		nargs = 6
// 	}
// 	//
// 	nc = 0
// 	reset = true
// 	errmax = 0
// 	//
// 	for in = 1; in <= nidim; in++ {
// 		n = (*idim)[in-1]
// 		//       Set LDA to 1 more than minimum value if room.
// 		lda = n
// 		if lda < nmax {
// 			lda = lda + 1
// 		}
// 		//       Skip tests if not enough room.
// 		if lda > nmax {
// 			goto Label100
// 		}
// 		if *PACKED {
// 			laa = (n * (n + 1)) / 2
// 		} else {
// 			laa = lda * n
// 		}
// 		//
// 		for ic = 1; ic <= 2; ic++ {
// 			(*UPLO) = ich[ic-1]
// 			upper = (*UPLO) == 'U'
// 			//
// 			for ix = 1; ix <= ninc; ix++ {
// 				incx = (*inc)[ix-1]
// 				lx = absint(incx) * n
// 				//
// 				//             Generate the vector X.
// 				//
// 				transl = 0.5
// 				smake("GE"), ' '), ' '), 1, n, x, 1, xx, absint(incx), 0, n-1, reset, transl)
// 				if n > 1 {
// 					(*x)[n/1] = 0
// 					(*xx)[1+absint(incx)*(n/2-1)-1] = 0
// 				}
// 				//
// 				for ia = 1; ia <= nalf; ia++ {
// 					alpha = (*alf)[ia-1]
// 					null = n <= 0 || alpha == 0
// 					//
// 					//                Generate the matrix A.
// 					//
// 					transl = 0
// 					smake(sname[1], UPLO, ' '), n, n, a, nmax, aa, lda, n-1, n-1, reset, transl)
// 					//
// 					nc = nc + 1
// 					//
// 					//                Save every datum before calling the subroutine.
// 					//
// 					(*UPLOS) = (*UPLO)
// 					ns = n
// 					als = alpha
// 					for i = 1; i <= laa; i++ {
// 						(*as)[i-1] = (*aa)[i-1]
// 					//Label10:
// 					}
// 					ldas = lda
// 					for i = 1; i <= lx; i++ {
// 						(*xs)[i-1] = (*xx)[i-1]
// 					//Label20:
// 					}
// 					incxs = incx
// 					//
// 					//                Call the subroutine.
// 					//
// 					if full {
// 						if trace {
// 							writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, x,%2d, a,%3d)                        .\n"), nc, sname, (*UPLO), n, alpha, incx, lda)
// 						}
// 						if rewi {
// 							rewind((ntra))
// 						}
// 						SSYR(UPLO, n, alpha, xx, incx, aa, LDA)
// 					} else if *PACKED {
// 						if trace {
// 							writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, x,%2d, AP)                           .\n"), nc, sname, (*UPLO), n, alpha, incx)
// 						}
// 						if rewi {
// 							rewind((ntra))
// 						}
// 						SSPR(UPLO, n, alpha, xx, incx, aa)
// 					}
// 					//
// 					//                Check if error-exit was taken incorrectly.
// 					//
// 					if !common.infoc.ok {
// 						writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n"))
// 						fatal = true
// 						goto Label120
// 					}
// 					//
// 					//                See what data changed inside subroutines.
// 					//
// 					isame[0] = (*UPLO) == (*UPLOS)
// 					isame[1] = ns == n
// 					isame[2] = als == alpha
// 					isame[3] = lse(xs, xx, lx))
// 					isame[4] = incxs == incx
// 					if null {
// 						isame[5] = lse(as, aa, laa))
// 					} else {
// 						isame[5] = (*lseres(sname[1], UPLO, n, n, as, aa, LDA))
// 					}
// 					if !(*PACKED) {
// 						isame[6] = ldas == lda
// 					}
// 					//
// 					//                If data was incorrectly changed, report and return.
// 					//
// 					same = true
// 					for i = 1; i <= nargs; i++ {
// 						same = same && isame[i-1]
// 						if !isame[i-1] {
// 							writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n"), i)
// 						}
// 					//Label30:
// 					}
// 					if !same {
// 						fatal = true
// 						goto Label120
// 					}
// 					//
// 					if !null {
// 						//
// 						//                   Check the result column by column.
// 						//
// 						if incx > 0 {
// 							for i = 1; i <= n; i++ {
// 								(*(Z))[i-1] = (*x)[i-1]
// 							//Label40:
// 							}
// 						} else {
// 							for i = 1; i <= n; i++ {
// 								(*(Z))[i-1] = (*x)[n-i+0]
// 							//Label50:
// 							}
// 						}
// 						(*JA) = 1
// 						for j = 1; j <= n; j++ {
// 							(*W)[0] = (*(Z))[j-1]
// 							if upper {
// 								(*JJ) = 1
// 								(*LJ) = j
// 							} else {
// 								(*JJ) = j
// 								(*LJ) = n - j + 1
// 							}
// 							smvch('N'), LJ, 1, alpha, &((*(Z))[(*JJ)-1]), LJ, W, 1, 1, &((*a)[(*JJ)-1][j-1]), 1, yt, g, &((*aa)[(*JA)-1]), eps, err, fatal, nout, true)
// 							if full {
// 								if upper {
// 									(*JA) = (*JA) + lda
// 								} else {
// 									(*JA) = (*JA) + lda + 1
// 								}
// 							} else {
// 								(*JA) = (*JA) + (*LJ)
// 							}
// 							errmax = (max(errmax, err))
// 							//                      If got really bad answer, report and return.
// 							if fatal {
// 								goto Label110
// 							}
// 						//Label60:
// 						}
// 					} else {
// 						//                   Avoid repeating tests if N.le.0.
// 						if n <= 0 {
// 							goto Label100
// 						}
// 					}
// 					//
// 				//Label70:
// 				}
// 				//
// 			//Label80:
// 			}
// 			//
// 		//Label90:
// 		}
// 		//
// 	Label100:
// 	}
// 	//
// 	//    Report result.
// 	//
// 	if errmax < thresh {
// 		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n"), sname, nc)
// 	} else {
// 		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n"), sname, nc, errmax)
// 	}
// 	goto Label130
// 	//
// Label110:
// 	;
// 	writeString(common.infoc.nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n"), j)
// 	//
// Label120:
// 	;
// 	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n"), sname)
// 	if full {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, x,%2d, a,%3d)                        .\n"), nc, sname, (*UPLO), n, alpha, incx, lda)
// 	} else if *PACKED {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, x,%2d, AP)                           .\n"), nc, sname, (*UPLO), n, alpha, incx)
// 	}
// 	//
// Label130:
// 	;
// 	return
// 	//
// 	//
// 	//    End of SCHK5.
// 	//
// }

// func schk6(SNAME []byte, EPS float32, THRESH float32, NOUT int, NTRA int, TRACE bool, REWI bool, FATAL bool, NIDIM int, IDIM *[]int, NALF int, ALF *[]float32, NINC int, INC *[]int, NMAX int, INCMAX int, A *[][]float32, AA *[]float32, AS *[]float32, X *[]float32, XX *[]float32, XS *[]float32, Y *[]float32, YY *[]float32, YS *[]float32, YT *[]float32, G *[]float32, Z *[ ][]float32) {
// 	ZERO := new(float32)
// 	HALF := new(float32)
// 	ONE := new(float32)
// 	ALPHA := new(float32)
// 	ALS := new(float32)
// 	ERR := new(float32)
// 	ERRMAX := new(float32)
// 	TRANSL := new(float32)
// 	I := new(int)
// 	IA := new(int)
// 	IC := new(int)
// 	IN := new(int)
// 	INCX := new(int)
// 	INCXS := new(int)
// 	INCY := new(int)
// 	INCYS := new(int)
// 	IX := new(int)
// 	IY := new(int)
// 	J := new(int)
// 	JA := new(int)
// 	JJ := new(int)
// 	LAA := new(int)
// 	LDA := new(int)
// 	LDAS := new(int)
// 	LJ := new(int)
// 	LX := new(int)
// 	LY := new(int)
// 	N := new(int)
// 	NARGS := new(int)
// 	NC := new(int)
// 	NS := new(int)
// 	FULL := new(bool)
// 	NULL := new(bool)
// 	PACKED := new(bool)
// 	RESET := new(bool)
// 	SAME := new(bool)
// 	UPPER := new(bool)
// 	UPLO := new(byte)
// 	UPLOS := new(byte)
// 	ICH := func() *[]byte {
// 		arr := make([]byte, 2)
// 		return &arr
// 	}()
// 	W := func() *[]float32 {
// 		arr := make([]float32, 2)
// 		return &arr
// 	}()
// 	ISAME := func() *[]bool {
// 		arr := make([]bool, 13)
// 		return &arr
// 	}()
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests SSYR2 and SSPR2.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Parameters ..
// 	0 = 0.0
// 	0.5 = 0.5
// 	1 = 1.0
// 	//    .. Scalar Arguments ..
// 	//    .. Array Arguments ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Functions ..
// 	//    .. External Subroutines ..
// 	//    .. Intrinsic Functions ..
// 	//    .. Scalars in Common ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Data statements ..
// 	ich[0], ich[1] = 'U', 'L'
// 	//    .. Executable Statements ..
// 	full = sname[2] == 'Y'
// 	(*PACKED) = sname[2] == 'P'
// 	//    Define the number of arguments.
// 	if full {
// 		nargs = 9
// 	} else if *PACKED {
// 		nargs = 8
// 	}
// 	//
// 	nc = 0
// 	reset = true
// 	errmax = 0
// 	//
// 	for in = 1; in <= nidim; in++ {
// 		n = (*idim)[in-1]
// 		//       Set LDA to 1 more than minimum value if room.
// 		lda = n
// 		if lda < nmax {
// 			lda = lda + 1
// 		}
// 		//       Skip tests if not enough room.
// 		if lda > nmax {
// 			goto Label140
// 		}
// 		if *PACKED {
// 			laa = (n * (n + 1)) / 2
// 		} else {
// 			laa = lda * n
// 		}
// 		//
// 		for ic = 1; ic <= 2; ic++ {
// 			(*UPLO) = ich[ic-1]
// 			upper = (*UPLO) == 'U'
// 			//
// 			for ix = 1; ix <= ninc; ix++ {
// 				incx = (*inc)[ix-1]
// 				lx = absint(incx) * n
// 				//
// 				//             Generate the vector X.
// 				//
// 				transl = 0.5
// 				smake("GE"), ' '), ' '), 1, n, x, 1, xx, absint(incx), 0, n-1, reset, transl)
// 				if n > 1 {
// 					(*x)[n/1] = 0
// 					(*xx)[1+absint(incx)*(n/2-1)-1] = 0
// 				}
// 				//
// 				for iy = 1; iy <= ninc; iy++ {
// 					incy = (*inc)[iy-1]
// 					ly = absint(incy) * n
// 					//
// 					//                Generate the vector Y.
// 					//
// 					transl = 0
// 					smake("GE"), ' '), ' '), 1, n, y, 1, yy, absint(incy), 0, n-1, reset, transl)
// 					if n > 1 {
// 						(*y)[n/1] = 0
// 						(*yy)[1+absint(incy)*(n/2-1)-1] = 0
// 					}
// 					//
// 					for ia = 1; ia <= nalf; ia++ {
// 						alpha = (*alf)[ia-1]
// 						null = n <= 0 || alpha == 0
// 						//
// 						//                   Generate the matrix A.
// 						//
// 						transl = 0
// 						smake(sname[1], UPLO, ' '), n, n, a, nmax, aa, lda, n-1, n-1, reset, transl)
// 						//
// 						nc = nc + 1
// 						//
// 						//                   Save every datum before calling the subroutine.
// 						//
// 						(*UPLOS) = (*UPLO)
// 						ns = n
// 						als = alpha
// 						for i = 1; i <= laa; i++ {
// 							(*as)[i-1] = (*aa)[i-1]
// 						//Label10:
// 						}
// 						ldas = lda
// 						for i = 1; i <= lx; i++ {
// 							(*xs)[i-1] = (*xx)[i-1]
// 						//Label20:
// 						}
// 						incxs = incx
// 						for i = 1; i <= ly; i++ {
// 							(*ys)[i-1] = (*yy)[i-1]
// 						//Label30:
// 						}
// 						incys = incy
// 						//
// 						//                   Call the subroutine.
// 						//
// 						if full {
// 							if trace {
// 								writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, x,%2d, y,%2d, a,%3d)                  .\n"), nc, sname, (*UPLO), n, alpha, incx, incy, lda)
// 							}
// 							if rewi {
// 								rewind((ntra))
// 							}
// 							SSYR2(UPLO, n, alpha, xx, incx, yy, incy, aa, LDA)
// 						} else if *PACKED {
// 							if trace {
// 								writeString(ntra, " %6d: %6s('%c',%3d,%4.1f, x,%2d, y,%2d, AP)                     .\n"), nc, sname, (*UPLO), n, alpha, incx, incy)
// 							}
// 							if rewi {
// 								rewind((ntra))
// 							}
// 							SSPR2(UPLO, n, alpha, xx, incx, yy, incy, aa)
// 						}
// 						//
// 						//                   Check if error-exit was taken incorrectly.
// 						//
// 						if !common.infoc.ok {
// 							writeString(common.infoc.nout, " ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******\n"))
// 							fatal = true
// 							goto Label160
// 						}
// 						//
// 						//                   See what data changed inside subroutines.
// 						//
// 						isame[0] = (*UPLO) == (*UPLOS)
// 						isame[1] = ns == n
// 						isame[2] = als == alpha
// 						isame[3] = lse(xs, xx, lx))
// 						isame[4] = incxs == incx
// 						isame[5] = lse(ys, yy, ly))
// 						isame[6] = incys == incy
// 						if null {
// 							isame[7] = lse(as, aa, laa))
// 						} else {
// 							isame[7] = (*lseres(sname[1], UPLO, n, n, as, aa, LDA))
// 						}
// 						if !(*PACKED) {
// 							isame[8] = ldas == lda
// 						}
// 						//
// 						//                   If data was incorrectly changed, report and return.
// 						//
// 						same = true
// 						for i = 1; i <= nargs; i++ {
// 							same = same && isame[i-1]
// 							if !isame[i-1] {
// 								writeString(common.infoc.nout, " ******* FATAL ERROR - PARAMETER NUMBER %2d WAS CHANGED INCORRECTLY *******\n"), i)
// 							}
// 						//Label40:
// 						}
// 						if !same {
// 							fatal = true
// 							goto Label160
// 						}
// 						//
// 						if !null {
// 							//
// 							//                      Check the result column by column.
// 							//
// 							if incx > 0 {
// 								for i = 1; i <= n; i++ {
// 									(*(Z))[i-1][0] = (*x)[i-1]
// 								//Label50:
// 								}
// 							} else {
// 								for i = 1; i <= n; i++ {
// 									(*(Z))[i-1][0] = (*x)[n-i+0]
// 								//Label60:
// 								}
// 							}
// 							if incy > 0 {
// 								for i = 1; i <= n; i++ {
// 									(*(Z))[i-1][1] = (*y)[i-1]
// 								//Label70:
// 								}
// 							} else {
// 								for i = 1; i <= n; i++ {
// 									(*(Z))[i-1][1] = (*y)[n-i+0]
// 								//Label80:
// 								}
// 							}
// 							(*JA) = 1
// 							for j = 1; j <= n; j++ {
// 								(*W)[0] = (*(Z))[j-1][1]
// 								(*W)[1] = (*(Z))[j-1][0]
// 								if upper {
// 									(*JJ) = 1
// 									(*LJ) = j
// 								} else {
// 									(*JJ) = j
// 									(*LJ) = n - j + 1
// 								}
// 								smvch('N'), LJ, 2, alpha, &((*(Z))[(*JJ)-1][0]), nmax, W, 1, 1, &((*a)[(*JJ)-1][j-1]), 1, yt, g, &((*aa)[(*JA)-1]), eps, err, fatal, nout, true)
// 								if full {
// 									if upper {
// 										(*JA) = (*JA) + lda
// 									} else {
// 										(*JA) = (*JA) + lda + 1
// 									}
// 								} else {
// 									(*JA) = (*JA) + (*LJ)
// 								}
// 								errmax = (max(errmax, err))
// 								//                         If got really bad answer, report and return.
// 								if fatal {
// 									goto Label150
// 								}
// 							//Label90:
// 							}
// 						} else {
// 							//                      Avoid repeating tests with N.le.0.
// 							if n <= 0 {
// 								goto Label140
// 							}
// 						}
// 						//
// 					//Label100:
// 					}
// 					//
// 				//Label110:
// 				}
// 				//
// 			//Label120:
// 			}
// 			//
// 		//Label130:
// 		}
// 		//
// 	Label140:
// 	}
// 	//
// 	//    Report result.
// 	//
// 	if errmax < thresh {
// 		writeString(common.infoc.nout, " %6s PASSED THE COMPUTATIONAL TESTS (%6d CALLS)\n"), sname, nc)
// 	} else {
// 		writeString(common.infoc.nout, " %6s COMPLETED THE COMPUTATIONAL TESTS (%6d CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO %8.2f - SUSPECT *******\n"), sname, nc, errmax)
// 	}
// 	goto Label170
// 	//
// Label150:
// 	;
// 	writeString(common.infoc.nout, "      THESE ARE THE RESULTS FOR COLUMN %3d\n"), j)
// 	//
// Label160:
// 	;
// 	writeString(common.infoc.nout, " ******* %6s FAILED ON CALL NUMBER:\n"), sname)
// 	if full {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, x,%2d, y,%2d, a,%3d)                  .\n"), nc, sname, (*UPLO), n, alpha, incx, incy, lda)
// 	} else if *PACKED {
// 		writeString(common.infoc.nout, " %6d: %6s('%c',%3d,%4.1f, x,%2d, y,%2d, AP)                     .\n"), nc, sname, (*UPLO), n, alpha, incx, incy)
// 	}
// 	//
// Label170:
// 	;
// 	return
// 	//
// 	//
// 	//    End of SCHK6.
// 	//
// }

// func schke(ISNUM int, SRNAMT []byte, NOUT int) {
// 	INFOT := new(int)
// 	NOUTC := new(int)
// 	LERR := new(bool)
// 	OK := new(bool)
// 	ALPHA := new(float32)
// 	BETA := new(float32)
// 	A := func() *[][]float32 {
// 		arr := make([][]float32, 1)
// 		for u := 0; u < 1; u++ {
// 			arr[u] = make([]float32, 1)
// 		}
// 		return &arr
// 	}()
// 	X := func() *[]float32 {
// 		arr := make([]float32, 1)
// 		return &arr
// 	}()
// 	Y := func() *[]float32 {
// 		arr := make([]float32, 1)
// 		return &arr
// 	}()
// 	COMMON.INFOC.LERR = new(float32)
// 	COMMON.INFOC.OK = new(float32)
// 	COMMON.INFOC.NOUTC = new(float32)
// 	COMMON.INFOC.INFOT = new(int)
// 	//
// 	// Tests the error exits from the Level 2 Blas.
// 	// Requires a special version of the error-handling routine Xerbla.
// 	// ALPHA, beta, a, X and Y should not need to be defined.
// 	//
// 	// Auxiliary routine for test program for Level 2 Blas.
// 	//
// 	// -- Written on 10-August-1987.
// 	//    Richard Hanson, Sandia National Labs.
// 	//    Jeremy Du Croz, NAG Central Office.
// 	//
// 	//    .. Scalar Arguments ..
// 	//    .. Scalars in Common ..
// 	//    .. Local Scalars ..
// 	//    .. Local Arrays ..
// 	//    .. External Subroutines ..
// 	//    .. Common blocks ..
// 	INFOT = COMMON.INFOC.INFOT
// 	NOUTC = COMMON.INFOC.NOUTC
// 	OK = COMMON.INFOC.OK
// 	LERR = COMMON.INFOC.LERR
// 	//    .. Executable Statements ..
// 	//    OK is set to .FALSE. by the special version of Xerbla or by CHKXER
// 	//    if anything is wrong.
// 	common.infoc.ok = true
// 	//    LERR is set to .TRUE. by the special version of Xerbla each time
// 	//    it is called, and is then tested and re-set by CHKXER.
// 	common.infoc.lerr = false
// 	switch *(ISNUM) {
// 	case 1:
// 		goto Label10
// 	case 2:
// 		goto Label20
// 	case 3:
// 		goto Label30
// 	case 4:
// 		goto Label40
// 	case 5:
// 		goto Label50
// 	case 6:
// 		goto Label60
// 	case 7:
// 		goto Label70
// 	case 8:
// 		goto Label80
// 	case 9:
// 		goto Label90
// 	case 10:
// 		goto Label100
// 	case 11:
// 		goto Label110
// 	case 12:
// 		goto Label120
// 	case 13:
// 		goto Label130
// 	case 14:
// 		goto Label140
// 	case 15:
// 		goto Label150
// 	case 16:
// 		goto Label160
// 	}
// Label10:
// 	;
// 	common.infoc.infot = 1
// 	Sgemv('/'), 0, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	Sgemv('N'), -1, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	Sgemv('N'), 0, -1, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 6
// 	Sgemv('N'), 2, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 8
// 	Sgemv('N'), 0, 0, alpha, a, 1, x, 0, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 11
// 	Sgemv('N'), 0, 0, alpha, a, 1, x, 1, beta, y, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label20:
// 	;
// 	common.infoc.infot = 1
// 	Sgbmv('/'), 0, 0, 0, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	Sgbmv('N'), -1, 0, 0, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	Sgbmv('N'), 0, -1, 0, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	Sgbmv('N'), 0, 0, -1, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	Sgbmv('N'), 2, 0, 0, -1, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 8
// 	Sgbmv('N'), 0, 0, 1, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 10
// 	Sgbmv('N'), 0, 0, 0, 0, alpha, a, 1, x, 0, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 13
// 	Sgbmv('N'), 0, 0, 0, 0, alpha, a, 1, x, 1, beta, y, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label30:
// 	;
// 	common.infoc.infot = 1
// 	SSYMV('/'), 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSYMV('U'), -1, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SSYMV('U'), 2, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	SSYMV('U'), 0, alpha, a, 1, x, 0, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 10
// 	SSYMV('U'), 0, alpha, a, 1, x, 1, beta, y, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label40:
// 	;
// 	common.infoc.infot = 1
// 	SSBMV('/'), 0, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSBMV('U'), -1, 0, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	SSBMV('U'), 0, -1, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 6
// 	SSBMV('U'), 0, 1, alpha, a, 1, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 8
// 	SSBMV('U'), 0, 0, alpha, a, 1, x, 0, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 11
// 	SSBMV('U'), 0, 0, alpha, a, 1, x, 1, beta, y, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label50:
// 	;
// 	common.infoc.infot = 1
// 	SSPMV('/'), 0, alpha, a, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSPMV('U'), -1, alpha, a, x, 1, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 6
// 	SSPMV('U'), 0, alpha, a, x, 0, beta, y, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 9
// 	SSPMV('U'), 0, alpha, a, x, 1, beta, y, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label60:
// 	;
// 	common.infoc.infot = 1
// 	STRMV('/'), 'N'), 'N'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STRMV('U'), '/'), 'N'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STRMV('U'), 'N'), '/'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STRMV('U'), 'N'), 'N'), -1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 6
// 	STRMV('U'), 'N'), 'N'), 2, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 8
// 	STRMV('U'), 'N'), 'N'), 0, a, 1, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label70:
// 	;
// 	common.infoc.infot = 1
// 	STBMV('/'), 'N'), 'N'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STBMV('U'), '/'), 'N'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STBMV('U'), 'N'), '/'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STBMV('U'), 'N'), 'N'), -1, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	STBMV('U'), 'N'), 'N'), 0, -1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	STBMV('U'), 'N'), 'N'), 0, 1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 9
// 	STBMV('U'), 'N'), 'N'), 0, 0, a, 1, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label80:
// 	;
// 	common.infoc.infot = 1
// 	STPMV('/'), 'N'), 'N'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STPMV('U'), '/'), 'N'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STPMV('U'), 'N'), '/'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STPMV('U'), 'N'), 'N'), -1, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	STPMV('U'), 'N'), 'N'), 0, a, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label90:
// 	;
// 	common.infoc.infot = 1
// 	STRSV('/'), 'N'), 'N'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STRSV('U'), '/'), 'N'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STRSV('U'), 'N'), '/'), 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STRSV('U'), 'N'), 'N'), -1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 6
// 	STRSV('U'), 'N'), 'N'), 2, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 8
// 	STRSV('U'), 'N'), 'N'), 0, a, 1, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label100:
// 	;
// 	common.infoc.infot = 1
// 	STBSV('/'), 'N'), 'N'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STBSV('U'), '/'), 'N'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STBSV('U'), 'N'), '/'), 0, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STBSV('U'), 'N'), 'N'), -1, 0, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	STBSV('U'), 'N'), 'N'), 0, -1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	STBSV('U'), 'N'), 'N'), 0, 1, a, 1, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 9
// 	STBSV('U'), 'N'), 'N'), 0, 0, a, 1, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label110:
// 	;
// 	common.infoc.infot = 1
// 	STPSV('/'), 'N'), 'N'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	STPSV('U'), '/'), 'N'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 3
// 	STPSV('U'), 'N'), '/'), 0, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 4
// 	STPSV('U'), 'N'), 'N'), -1, a, x, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	STPSV('U'), 'N'), 'N'), 0, a, x, 0)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label120:
// 	;
// 	common.infoc.infot = 1
// 	SGER(-1, 0, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SGER(0, -1, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SGER(0, 0, alpha, x, 0, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	SGER(0, 0, alpha, x, 1, y, 0, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 9
// 	SGER(2, 0, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label130:
// 	;
// 	common.infoc.infot = 1
// 	SSYR('/'), 0, alpha, x, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSYR('U'), -1, alpha, x, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SSYR('U'), 0, alpha, x, 0, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	SSYR('U'), 2, alpha, x, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label140:
// 	;
// 	common.infoc.infot = 1
// 	SSPR('/'), 0, alpha, x, 1, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSPR('U'), -1, alpha, x, 1, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SSPR('U'), 0, alpha, x, 0, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label150:
// 	;
// 	common.infoc.infot = 1
// 	SSYR2('/'), 0, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSYR2('U'), -1, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SSYR2('U'), 0, alpha, x, 0, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	SSYR2('U'), 0, alpha, x, 1, y, 0, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 9
// 	SSYR2('U'), 2, alpha, x, 1, y, 1, a, 1)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	goto Label170
// Label160:
// 	;
// 	common.infoc.infot = 1
// 	SSPR2('/'), 0, alpha, x, 1, y, 1, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 2
// 	SSPR2('U'), -1, alpha, x, 1, y, 1, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 5
// 	SSPR2('U'), 0, alpha, x, 0, y, 1, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	common.infoc.infot = 7
// 	SSPR2('U'), 0, alpha, x, 1, y, 0, A)
// 	chkxer((SRNAMT), INFOT, nout, LERR, OK)
// 	//
// Label170:
// 	;
// 	if *OK {
// 		writeString(common.infoc.nout, " %6s PASSED THE TESTS OF ERROR-EXITS\n"), common.srnamc.srnamt)
// 	} else {
// 		writeString(common.infoc.nout, " ******* %6s FAILED THE TESTS OF ERROR-EXITS *******\n"), common.srnamc.srnamt)
// 	}
// 	return
// 	//
// 	//
// 	//    End of SCHKE.
// 	//
// }

// VERIFIED
func smake(_type string, uplo byte, diag byte, m int, n int, a *[][]float32, nmax int, aa *[]float32, lda int, kl int, ku int, reset bool, transl float32) {
	var rogue float32
	var i, i1, i2, i3, ibeg, iend, ioff, j, kk int
	var gen, lower, sym, tri, unit, upper bool
	//
	// Generates values for an M by N matrix A within the bandwidth
	// defined by KL and KU.
	// Stores the values in the array AA in the data structure required
	// by the routine, with unwanted elements set to rogue value.
	//
	// TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	rogue = -1.0e10
	gen = _type[0] == 'G'
	sym = _type[0] == 'S'
	tri = _type[0] == 'T'
	upper = (sym || tri) && uplo == 'U'
	lower = (sym || tri) && uplo == 'L'
	unit = tri && diag == 'U'
	//
	//    Generate data in array A.
	//
	for j = 1; j <= n; j++ {
		for i = 1; i <= m; i++ {
			if gen || (upper && i <= j) || (lower && i >= j) {
				if (i <= j && j-i <= ku) || (i >= j && i-j <= kl) {
					(*a)[i-1][j-1] = sbeg(reset) + transl
				} else {
					(*a)[i-1][j-1] = 0
				}
				if i != j {
					if sym {
						(*a)[j-1][i-1] = (*a)[i-1][j-1]
					} else if tri {
						(*a)[j-1][i-1] = 0
					}
				}
			}
			//Label10:
		}
		if tri {
			(*a)[j-1][j-1]++
		}
		if unit {
			(*a)[j-1][j-1] = 1
		}
		//Label20:
	}
	//
	//    Store elements in array AS in data structure required by routine.
	//
	if _type == "GE" {
		for j = 1; j <= n; j++ {
			for i = 1; i <= m; i++ {
				(*aa)[i+(j-1)*lda-1] = (*a)[i-1][j-1]
				//Label30:
			}
			for i = m + 1; i <= lda; i++ {
				(*aa)[i+(j-1)*lda-1] = rogue
				//Label40:
			}
			//Label50:
		}
	} else if _type == "GB" {
		for j = 1; j <= n; j++ {
			for i1 = 1; i1 <= ku+1-j; i1++ {
				(*aa)[i1+(j-1)*lda-1] = rogue
				//Label60:
			}
			for i2 = i1; i2 <= min(kl+ku+1, ku+1+m-j); i2++ {
				(*aa)[i2+(j-1)*lda-1] = (*a)[i2+j-ku-2][j-1]
				//Label70:
			}
			for i3 = i2; i3 <= lda; i3++ {
				(*aa)[i3+(j-1)*lda-1] = rogue
				//Label80:
			}
			//Label90:
		}
	} else if _type == "SY" || _type == "TR" {
		for j = 1; j <= n; j++ {
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
				iend = n
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*lda-1] = rogue
				//Label100:
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*lda-1] = (*a)[i-1][j-1]
				//Label110:
			}
			for i = iend + 1; i <= lda; i++ {
				(*aa)[i+(j-1)*lda-1] = rogue
				//Label120:
			}
			//Label130:
		}
	} else if _type == "SB" || _type == "TB" {
		for j = 1; j <= n; j++ {
			if upper {
				kk = kl + 1
				ibeg = max(1, kl+2-j)
				if unit {
					iend = kl
				} else {
					iend = kl + 1
				}
			} else {
				kk = 1
				if unit {
					ibeg = 2
				} else {
					ibeg = 1
				}
				iend = min(kl+1, 1+m-j)
			}
			for i = 1; i <= ibeg-1; i++ {
				(*aa)[i+(j-1)*lda-1] = rogue
				//Label140:
			}
			for i = ibeg; i <= iend; i++ {
				(*aa)[i+(j-1)*lda-1] = (*a)[i+j-kk-1][j-1]
				//Label150:
			}
			for i = iend + 1; i <= lda; i++ {
				(*aa)[i+(j-1)*lda-1] = rogue
				//Label160:
			}
			//Label170:
		}
	} else if _type == "SP" || _type == "TP" {
		ioff = 0
		for j = 1; j <= n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = n
			}
			for i = ibeg; i <= iend; i++ {
				ioff++
				(*aa)[ioff-1] = (*a)[i-1][j-1]
				if i == j {
					if unit {
						(*aa)[ioff-1] = rogue
					}
				}
				//Label180:
			}
			//Label190:
		}
	}
	return
	//
	//    End of SMAKE.
	//
}

// VERIFIED
func smvch(trans byte, m int, n int, alpha float32, a *[][]float32, nmax int, x *[]float32, incx int, beta float32, y *[]float32, incy int, yt *[]float32, g *[]float32, yy *[]float32, eps float32, err float32, fatal bool, nout int, mv bool) {
	var erri float32
	var i, incxl, incyl, iy, j, jx, kx, ky, ml, nl int
	var tran bool
	//
	// Checks the results of the computational tests.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	tran = trans == 'T' || trans == 'C'
	if tran {
		ml = n
		nl = m
	} else {
		ml = m
		nl = n
	}
	if incx < 0 {
		kx = nl
		incxl = -1
	} else {
		kx = 1
		incxl = 1
	}
	if incy < 0 {
		ky = ml
		incyl = -1
	} else {
		ky = 1
		incyl = 1
	}
	//
	//    Compute expected result in YT using data in A, X and Y.
	//    Compute gauges in G.
	//
	iy = ky
	for i = 1; i <= ml; i++ {
		(*yt)[iy-1] = 0
		(*g)[iy-1] = 0
		jx = kx
		if tran {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[j-1][i-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[j-1][i-1] * (*x)[jx-1])
				jx += incxl
				//Label10:
			}
		} else {
			for j = 1; j <= nl; j++ {
				(*yt)[iy-1] += (*a)[i-1][j-1] * (*x)[jx-1]
				(*g)[iy-1] += absf32((*a)[i-1][j-1] * (*x)[jx-1])
				jx += incxl
				//Label20:
			}
		}
		(*yt)[iy-1] = alpha*(*yt)[iy-1] + beta*(*y)[iy-1]
		(*g)[iy-1] = absf32(alpha)*(*g)[iy-1] + absf32(beta*(*y)[iy-1])
		iy += incyl
		//Label30:
	}
	//
	//    Compute the error ratio for this result.
	//
	err = 0
	for i = 1; i <= ml; i++ {
		erri = absf32((*yt)[i-1]-(*yy)[1+(i-1)*absint(incy)-1]) / eps
		if (*g)[i-1] != 0 {
			erri /= (*g)[i-1]
		}
		err = maxf32(err, erri)
		if err*sqrtf32(eps) >= 1 {
			goto Label50
		}
		//Label40:
	}
	//    If the loop completes, all results are at least half accurate.
	goto Label70
	//
	//    Report fatal error.
	//
Label50:
	;
	fatal = true
	writeString(nout, " ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n")
	writeString(nout, " %6d: smvch('%c', %3d, %3d, %4.1f, a, %3d, x, %2d, %4.1f, y, %2d, yt, g, yy, %9.1e, %4.1f, %t, %2d, %t)\n", common.combla.n, trans, m, n, alpha, nmax, incx, beta, incy, eps, err, fatal, nout, mv)
	writeString(nout, "            EXPECTED RESULT    COMPUTED RESULT\n")
	for i = 1; i <= ml; i++ {
		if mv {
			writeString(nout, " %7d %18.8f %18.8f\n", i, (*yt)[i-1], (*yy)[1+(i-1)*absint(incy)-1])
		} else {
			writeString(nout, " %7d %18.8f %18.8f\n", i, (*yy)[1+(i-1)*absint(incy)-1], (*yt)[i-1])
		}
		//Label60:
	}
	writeString(nout, "\n")
	//
Label70:
	;
	return
	//
	//
	//    End of SMVCH.
	//
}

// VERIFIED
func lse(ri, rj *[]float32, lr int) bool {
	var i int
	//
	// Tests if two arrays are identical.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	//    .. Scalar Arguments ..
	//    .. Array Arguments ..
	//    .. Local Scalars ..
	//    .. Executable Statements ..
	for i = 1; i <= lr; i++ {
		if (*ri)[i-1] != (*rj)[i-1] {
			return false
		}
		//Label10:
	}
	return true
	//
	//    End of LSE.
	//
}

// VERIFIED
func lseres(_type string, uplo byte, m int, n int, aa *[][]float32, as *[][]float32, lda int) bool {
	var i, ibeg, iend, j int
	var upper bool
	//
	// Tests if selected elements in two arrays are equal.
	//
	// TYPE is 'GE', 'SY' or 'SP'.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	upper = uplo == 'U'
	if _type == "GE" {
		for j = 1; j <= n; j++ {
			for i = m + 1; i <= lda; i++ {
				if (*aa)[j-1][i-1] != (*as)[j-1][i-1] {
					goto Label70
				}
				//Label10:
			}
			//Label20:
		}
	} else if _type == "SY" {
		for j = 1; j <= n; j++ {
			if upper {
				ibeg = 1
				iend = j
			} else {
				ibeg = j
				iend = n
			}
			for i = 1; i <= ibeg-1; i++ {
				if (*aa)[j-1][i-1] != (*as)[j-1][i-1] {
					goto Label70
				}
				//Label30:
			}
			for i = iend + 1; i <= lda; i++ {
				if (*aa)[j-1][i-1] != (*as)[j-1][i-1] {
					goto Label70
				}
				//Label40:
			}
			//Label50:
		}
	}
	//
	return true
Label70:
	;
	return false
	//
	//    End of LSERES.
	//
}

// VERIFIED
func sbeg(reset bool) float32 {
	var i, ic, mi int
	//
	// Generates random numbers uniformly distributed between -0.5 and 0.5.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	//    .. Scalar Arguments ..
	//    .. Local Scalars ..
	//    .. Save statement ..
	//    .. Intrinsic Functions ..
	//    .. Executable Statements ..
	if reset {
		//       Initialize local variables.
		mi = 891
		i = 7
		ic = 0
		reset = false
	}
	//
	//    The sequence of values of I is bounded between 1 and 999.
	//    If initial I = 1,2,3,6,7 or 9, the period will be 50.
	//    If initial I = 4 or 8, the period will be 25.
	//    If initial I = 5, the period will be 10.
	//    IC is used to break up the period by skipping 1 value of I in 6.
	//
	ic++
Label10:
	;
	i *= mi
	i -= 1000 * (i / 1000)
	if ic >= 5 {
		ic = 0
		goto Label10
	}
	return float32(i-500) / 1001.0
	//
	//    End of SBEG.
	//
}

// VERIFIED
func chkxer() {
	//
	// Tests whether Xerbla has detected an error when it should.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	//    .. Scalar Arguments ..
	//    .. Executable Statements ..
	if !common.infoc.lerr {
		writeString(common.infoc.nout, " ***** ILLEGAL VALUE OF PARAMETER NUMBER %2d NOT DETECTED BY %6s *****\n", common.infoc.infot, common.srnamc.srnamt)
		common.infoc.ok = false
	}
	common.infoc.lerr = false
	return
	//
	//
	//    End of CHKXER.
	//
}

// Xerbla ...
// VERIFIED
func Xerbla(srname string, info int) {
	//
	// This is a special version of Xerbla to be used only as part of
	// the test program for testing error exits from the Level 2 BLAS
	// routines.
	//
	// Xerbla  is an error handler for the Level 2 BLAS routines.
	//
	// It is called by the Level 2 BLAS routines if an input parameter is
	// invalid.
	//
	// Auxiliary routine for test program for Level 2 Blas.
	//
	// -- Written on 10-August-1987.
	//    Richard Hanson, Sandia National Labs.
	//    Jeremy Du Croz, NAG Central Office.
	//
	common.infoc.lerr = true
	if info != common.infoc.infot {
		if common.infoc.infot != 0 {
			writeString(common.infoc.nout, " ******* Xerbla WAS CALLED WITH INFO = %6d INSTEAD OF %2d *******\n", info, common.infoc.infot)
		} else {
			writeString(common.infoc.nout, " ******* Xerbla WAS CALLED WITH INFO = %6d *******\n", info)
		}
		common.infoc.ok = false
	}
	if srname != common.srnamc.srnamt {
		writeString(common.infoc.nout, " ******* Xerbla WAS CALLED WITH SRNAME = %6s INSTEAD OF %6s *******\n", srname, common.srnamc.srnamt)
		common.infoc.ok = false
	}
	return
	//
	//
	//    End of Xerbla
	//
}

func chkerr(err error) {
	if err != nil {
		_ = fmt.Errorf("%v", err)
	}
}
