package goblas

import "time"

// Dsecnd Using the inTERNAL function Etime.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//      DOUBLE PRECISION FUNCTION Dsecnd()
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//  Dsecnd returns the user time for a process in seconds.
//  This version gets the time from the intERNAL function Etime.
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
// \ingroup auxOTHERauxiliary
//
//  =====================================================================
func Dsecnd(start *time.time) (dsecndReturn *float64) {
	dsecndReturn = new(float64)
	t1 := new(float64)
	tarray := func() *[]float64 {
		arr := make([]float64, 2)
		return &arr
	}()
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	// =====================================================================
	//
	//     .. Local Scalars ..
	//     ..
	//     .. Local Arrays ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	*t1 = time.Since(*start).Seconds()
	*dsecndReturn = (*tarray)[0]
	return
	//
	//     End of Dsecnd
	//
}
