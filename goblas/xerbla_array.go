package goblas

import 

// \brief \b XerblaArray
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//  Definition:
//  ===========
//
//       SUBROUTINE XerblaArray(SRNAME_ARRAY, SRNAME_LEN, INFO)
//
//       .. Scalar Arguments ..
//       INTEGER SRNAME_LEN, INFO
//       ..
//       .. Array Arguments ..
//       CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)
//       ..
//
//
// \par Purpose:
//  =============
//
// \verbatim
//
// XerblaArray assists other languages in calling Xerbla, the LAPACK
// and BLAS error handler.  Rather than taking a Fortran string argument
// as the function's name, XerblaArray takes an array of single
// characters along with the array's length.  XerblaArray then copies
// up to 32 characters of that array into a Fortran string and passes
// that to Xerbla.  If called with a non-positive SRNAME_LEN,
// XerblaArray will call Xerbla with a string of all blank characters.
//
// Say some macro or other device makes XerblaArray available to C99
// by a name lapack_xerbla and with a common Fortran calling convention.
// Then a C99 program could invoke Xerbla via:
//    {
//      int flen = strlen(__func__);
//      lapack_xerbla(__func__, &flen, &info);
//    }
//
// Providing XerblaArray is not necessary for intercepting LAPACK
// errors.  XerblaArray calls Xerbla.
// \endverbatim
//
//  Arguments:
//  ==========
//
// \param[in] SRNAME_ARRAY
// \verbatim
//          SRNAME_ARRAY is CHARACTER(1) array, dimension (SRNAME_LEN)
//          The name of the routine which called XerblaArray.
// \endverbatim
//
// \param[in] SRNAME_LEN
// \verbatim
//          SRNAME_LEN is INTEGER
//          The length of the name in SRNAME_ARRAY.
// \endverbatim
//
// \param[in] INFO
// \verbatim
//          INFO is INTEGER
//          The position of the invalid parameter in the parameter list
//          of the calling routine.
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
func Xerbla_array(srname_array *[]byte, srname_len *int, info *int) {
	i := new(int)
	srname := func() *[]byte {
		arr := make([]byte, 32)
		return &arr
	}()
	//*
	//*  -- Reference BLAS level1 routine (version 3.7.0) --
	//*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
	//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//*     December 2016
	//*
	//*     .. Scalar Arguments ..
	//*     ..
	//*     .. Array Arguments ..
	//*     ..
	//*
	//* =====================================================================
	//*
	//*     ..
	//*     .. Local Scalars ..
	//*     ..
	//*     .. Local Arrays ..
	//*     ..
	//*     .. Intrinsic Functions ..
	//*     ..
	//*     .. External Functions ..
	//*     ..
	//*     .. Executable Statements ..
	(*srname) = *func() *[]byte {y :=[]byte(""); return &y }()
	for (*i) = 1; (*i) <= (MIN((*srname_len), len((*srname)))); (*i)++ {
		(*srname)[(*i)-(1)] = (*srname_array)[(*i)-(1)]
	}
	Xerbla(srname, info)
	return
}
