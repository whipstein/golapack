package goblas

// Ilaver returns the lapACK version.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//     SUBROUTinE ILAVEr( versMajor, versMinor, versPatch)
//
//     intEGER versMajor, versMinor, versPatch
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
//  This subroutine returns the lapACK version.
// \endverbatim
//
//  Arguments:
//  ==========
//
//  \param[out] versMajor
//      versMajor is intEGER
//      return the lapack major version
//
//  \param[out] versMinor
//      versMinor is intEGER
//      return the lapack minor version from the major version
//
//  \param[out] versPatch
//      versPatch is intEGER
//      return the lapack patch version from the minor version
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
// \ingroup auxOTHERauxiliary
//
//  =====================================================================
func Ilaver(versMajor *int, versMinor *int, versPatch *int) {
	//
	//  -- lapACK computational routine --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//
	//  =====================================================================
	//
	//  =====================================================================
	*(versMajor) = 3
	*(versMinor) = 9
	*(versPatch) = 0
	//  =====================================================================
	//
	return
}
