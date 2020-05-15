package goblas

import 

// Alareq handles input for the lapACK test program.  It is called
// to evaluate the input line which requested nmats matrix types for
// path.  The flow of control is as follows:
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTinE Alareq( path, nmats, dotype, ntypes, nin, nout)
//
//       .. Scalar Arguments ..
//       CHARACTER*3        path
//       intEGER            nin, nmats, nout, ntypes
//       ..
//       .. Array Arguments ..
//       LOGICAL            dotype(*)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Alareq handles input for the lapACK test program.  It is called
// to evaluate the input line which requested nmats matrix types for
// path.  The flow of control is as follows:
//
// If nmats = ntypes then
//    dotype(1:ntypes) = .TRUE.
// else
//    Read the next input line for nmats matrix types
//    Set dotype(i) = .TRUE. for each valid type I
// endif
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] path
// \verbatim
//          path is CHARACTER*3
//          An lapACK path name for testing.
// \endverbatim
//
// \param[in] nmats
// \verbatim
//          nmats is intEGER
//          The number of matrix types to be used in testing this path.
// \endverbatim
//
// \param[out] dotype
// \verbatim
//          dotype is LOGICAL array, dimension (ntypes)
//          The vector of flags indicating if each type will be tested.
// \endverbatim
//
// \param[in] ntypes
// \verbatim
//          ntypes is intEGER
//          The maximum number of matrix types for this path.
// \endverbatim
//
// \param[in] nin
// \verbatim
//          nin is intEGER
//          The unit number for input.  nin >= 1.
// \endverbatim
//
// \param[in] nout
// \verbatim
//          nout is intEGER
//          The unit number for output.  nout >= 1.
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
// \date December 2016
//
// \ingroup aux_lin
//
//  =====================================================================
func Alareq(path *[]byte, nmats *int, dotype *[]bool, ntypes *int, nin *int, nout *int) {
	firstt := new(bool)
	c1 := new(byte)
	intstr := func() *[]byte {
		arr := make([]byte, 10)
		return &arr
	}()
	line := func() *[]byte {
		arr := make([]byte, 80)
		return &arr
	}()
	i := new(int)
	i1 := new(int)
	ic := new(int)
	j := new(int)
	k := new(int)
	lenp := new(int)
	nt := new(int)
	nreq := func() *[]int {
		arr := make([]int, 100)
		return &arr
	}()
	//
	//  -- lapACK test routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//     .. Array Arguments ..
	//     ..
	//
	//  =====================================================================
	//
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Data statements ..
	(*intstr)[0], (*intstr)[1], (*intstr)[2], (*intstr)[3], (*intstr)[4], (*intstr)[5], (*intstr)[6], (*intstr)[7], (*intstr)[8], (*intstr)[9] = '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'
	//     ..
	//     .. Executable Statements ..
	//
	if (*(nmats)) >= (*(ntypes)) {
		//
		//        Test everything if nmats >= ntypes.
		//
		for (*i) = 1; (*i) <= (*(ntypes)); (*i)++ {
			(*(dotype))[(*i)-(1)] = true
			//Label10:
		}
	} else {
		for (*i) = 1; (*i) <= (*(ntypes)); (*i)++ {
			(*(dotype))[(*i)-(1)] = false
			//Label20:
		}
		(*firstt) = true
		//
		//        Read a line of matrix types if 0 < nmats < ntypes.
		//
		if (*(nmats)) > 0 {
			REAd((nin), []byte("%80s\n"), line)
			(*lenp) = len((*line))
			(*i) = 0
			for (*j) = 1; (*j) <= (*(nmats)); (*j)++ {
				(*nreq)[(*j)-(1)] = 0
				(*i1) = 0
			Label30:
				;
				(*i) = (*i) + 1
				if (*i) > (*lenp) {
					if (*j) == (*(nmats)) && (*i1) > 0 {
						goto Label60
					} else {
						WRITE((*(nout)), *func() *[]byte {y :=[]byte("%v *** Not enough matrix types on input line\n%79s\n"); return &y }(), (*line))
						WRITE((*(nout)), *func() *[]byte {
							y :=[]byte(" ==> Specify %4d matrix types on this line or adjust ntypes on previous line\n")
							return &y
						}(), (*(nmats)))
						goto Label80
					}
				}
				if (*line)[(*i)-(1)] != ' ' && (*line)[(*i)-(1)] != ',' {
					(*i1) = (*i)
					(*c1) = (*line)[(*i1)-(1)]
					//
					//              Check that a valid integer was read
					//
					for (*k) = 1; (*k) <= 10; (*k)++ {
						if (*c1) == (*intstr)[(*k)-(1)] {
							(*ic) = (*k) - 1
							goto Label50
						}
						//Label40:
					}
					WRITE((*(nout)), *func() *[]byte {
						y :=[]byte("%v *** Invalid integer value in column %2d of input line:\n%79s\n")
						return &y
					}(), (*i), (*line))
					WRITE((*(nout)), *func() *[]byte {
						y :=[]byte(" ==> Specify %4d matrix types on this line or adjust ntypes on previous line\n")
						return &y
					}(), (*(nmats)))
					goto Label80
				Label50:
					;
					(*nreq)[(*j)-(1)] = 10*(*nreq)[(*j)-(1)] + (*ic)
					goto Label30
				} else if (*i1) > 0 {
					goto Label60
				} else {
					goto Label30
				}
			Label60:
			}
		}
		for (*i) = 1; (*i) <= (*(nmats)); (*i)++ {
			(*nt) = (*nreq)[(*i)-(1)]
			if (*nt) > 0 && (*nt) <= (*(ntypes)) {
				if (*(dotype))[(*nt)-(1)] {
					if *firstt {
						WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %v\n"); return &y }())
					}
					(*firstt) = false
					WRITE((*(nout)), *func() *[]byte {
						y :=[]byte(" *** warning:  duplicate request of matrix type %2d for %3s\n")
						return &y
					}(), (*nt), (*(path)))
				}
				(*(dotype))[(*nt)-(1)] = true
			} else {
				WRITE((*(nout)), *func() *[]byte {
					y :=[]byte(" *** Invalid type request for %3s, type  %4d: must satisfy  1 <= type <= %2d\n")
					return &y
				}(), (*(path)), (*nt), (*(ntypes)))
			}
			//Label70:
		}
	Label80:
	}
	return
	//
	//Label90:

	WRITE((*(nout)), *func() *[]byte {
		y :=[]byte("\n *** End of file reached when trying to read matrix types for %3s\n *** Check that you are requesting the right number of types for each path\n\n")
		return &y
	}(), (*(path)))
	WRITE((*(nout)), *func() *[]byte {y :=[]byte(" %v\n"); return &y }())
	panic("")
	//
	//     End of Alareq
	//
}
