package goblas

import "github.com/whipstein/golapack/blas"

// Lsamen tests if the first n letters of ca are the same as the
// first n letters of cb, regardless of case.
// Lsamen returns trud if ca and cb are equivalent except for case
// and false otherwise.  Lsamen also returns false if len(ca)
// or len(cb) is less than n.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
// \htmlonly
// Download Lsamen + dependencies
// <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/Lsamen.f">
//[TGZ]</a>
// <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/Lsamen.f">
//[Zip]</a>
// <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/Lsamen.f">
//[Txt]</a>
// \endhtmlonly
//
//  Definition:
//  ===========
//
//       func Lsamen(n, ca, cb) bool
//
//       .. Scalar Arguments ..
//       *[]byte    ca, cb
//       int            n
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Lsamen  tests if the first n letters of ca are the same as the
// first n letters of cb, regardless of case.
// Lsamen returns trud if ca and cb are equivalent except for case
// and false otherwise.  Lsamen also returns false if len(ca)
// or len(cb) is less than n.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] n
// \verbatim
//          n is int
//          The number of characters in ca and cb to be compared.
// \endverbatim
//
// \param[in] ca
// \verbatim
//          ca is *[]byte
// \endverbatim
//
// \param[in] cb
// \verbatim
//          cb is *[]byte
//          ca and cb specify two character strings of length at least n.
//          Only the first n characters of each string will be accessed.
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
// \ingroup OTHERauxiliary
//
//  =====================================================================
func Lsamen(n *int, ca *[]byte, cb *[]byte) (LsamenReturn bool) {
	//
	//  -- lapACK auxiliary routine (version 3.7.0) --
	//  -- lapACK is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	// =====================================================================
	//
	//     .. Local Scalars ..
	//     ..
	//     .. External Functions ..
	//     ..
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Executable Statements ..
	//
	LsamenReturn = false
	if len(*ca) < *n || len(*cb) < *n {
		return
	}
	//
	//     Do for each character in the two strings.
	//
	for i := 0; i < *n; i++ {
		//
		//        Test if the characters are equal using blas.Lsame.
		//
		if !blas.Lsame(&(*ca)[i], &(*cb)[i]) {
			return
		}
		//
		//Label10:
	}
	LsamenReturn = true
	//
	return
	//
	//     End of Lsamen
	//
}
