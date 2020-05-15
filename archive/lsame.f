*> \brief \b lsame
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       LOGicaL FUNCTION lsame(ca,cb)
*
*       .. Scalar Arguments ..
*       CHARACTER ca,cb
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> lsame returns .TRUE. if ca is the same letter as cb regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ca
*> \verbatim
*>          ca is CHARACTER*1
*> \endverbatim
*>
*> \param[in] cb
*> \verbatim
*>          cb is CHARACTER*1
*>          ca and cb specify the single characters to be compared.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup aux_blas
*
*  =====================================================================
      FUNCTION lsame(ca,cb) result(lsame) bind(c)
            use iso_c_binding, only: c_int
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER ca,cb
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      intRinSIC ICHAR
*     ..
*     .. Local Scalars ..
      intEGER intA,intB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      lsame = ca .EQ. cb
      IF (lsame) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAr('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAr('A') on Prime machines returns 193 which is the same as
*     ICHAr('A') on an EBCDIC machine.
*
      intA = ICHAr(ca)
      intB = ICHAr(cb)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (intA.GE.97 .AND. intA.LE.122) intA = intA - 32
          IF (intB.GE.97 .AND. intB.LE.122) intB = intB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (intA.GE.129 .AND. intA.LE.137 .OR.
     +        intA.GE.145 .AND. intA.LE.153 .OR.
     +        intA.GE.162 .AND. intA.LE.169) intA = intA + 64
          IF (intB.GE.129 .AND. intB.LE.137 .OR.
     +        intB.GE.145 .AND. intB.LE.153 .OR.
     +        intB.GE.162 .AND. intB.LE.169) intB = intB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (intA.GE.225 .AND. intA.LE.250) intA = intA - 32
          IF (intB.GE.225 .AND. intB.LE.250) intB = intB - 32
      END IF
      lsame = intA .EQ. intB
*
*     RETURN
*
*     End of lsame
*
      END
