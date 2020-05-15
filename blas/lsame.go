package blas

// Lsame returns true if ca is the same letter as cb regardless of case.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       func Lsame(ca,cb *byte) bool
//
//       .. Scalar Arguments ..
//       byte ca,cb
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// Lsame returns true if ca is the same letter as cb regardless of
// case.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] ca
// \verbatim
//          ca is *byte
// \endverbatim
//
// \param[in] cb
// \verbatim
//          cb is *byte
//          ca and cb specify the single bytes to be compared.
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
// \ingroup aux_blas
//
//  =====================================================================
func Lsame(ca, cb *byte) (lsameReturn bool) {
	// var inta, intb, zcode rune
	var inta, intb, zcode int
	//
	//  -- Reference BLAS level1 routine (version 3.1) --
	//  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//     December 2016
	//
	//     .. Scalar Arguments ..
	//     ..
	//
	// =====================================================================
	//
	//     .. Intrinsic Functions ..
	//     ..
	//     .. Local Scalars ..
	//     ..
	//
	//     Test if the characters are equal
	//
	lsameReturn = ca == cb
	if lsameReturn {
		return
	}
	//
	//     Now test for equivalence if both characters are alphabetic.
	//
	zcode = int(func() byte { var y byte; y = byte('Z'); return y }())
	//
	//     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
	//     machines, on which int returns a value with bit 8 set.
	//     int('A') on Prime machines returns 193 which is the same as
	//     int('A') on an EBCDIC machine.
	//
	inta = int(*ca)
	intb = int(*cb)
	//
	if zcode == 90 || zcode == 122 {
		//
		//        ASCII is assumed - zcode is the ASCII code of either lower or
		//        upper case 'Z'.
		//
		if inta >= 97 && inta <= 122 {
			inta -= 32
		}
		if intb >= 97 && intb <= 122 {
			intb -= 32
		}
		//
	} else if zcode == 233 || zcode == 169 {
		//
		//        EBCDIC is assumed - zcode is the EBCDIC code of either lower or
		//        upper case 'Z'.
		//
		if inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta >= 162 && inta <= 169 {
			inta += 64
		}
		if intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb >= 162 && intb <= 169 {
			intb += 64
		}
		//
	} else if zcode == 218 || zcode == 250 {
		//
		//        ASCII is assumed, on Prime machines - zcode is the ASCII code
		//        plus 128 of either lower or upper case 'Z'.
		//
		if inta >= 225 && inta <= 250 {
			inta -= 32
		}
		if intb >= 225 && intb <= 250 {
			intb -= 32
		}
	}
	lsameReturn = inta == intb
	//
	//     return
	//
	//     End of Lsame
	//
	return
}
